#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) University of St Andrews 2020-2021
# (c) University of Strathclyde 2020-2021
# Author:
# Emma E. M. Hobbs
#
# Contact
# eemh1@st-andrews.ac.uk
#
# Emma E. M. Hobbs,
# Biomolecular Sciences Building,
# University of St Andrews,
# North Haugh Campus,
# St Andrews,
# KY16 9ST
# Scotland,
# UK
#
# The MIT License
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#
"""Parse dbCAN output and a local CAZyme database to compile an annotated CAZome"""


import json
import logging
import re
import sys
import yaml

from datetime import datetime
from pathlib import Path
from typing import List, Optional


from cazy_webscraper.sql.sql_orm import (
  CazyFamily,
  Genbank,
  Session,
  get_db_connection,
)
from saintBioutils.utilities.file_io import get_paths, make_output_directory
from saintBioutils.utilities.logger import config_logger

from tqdm import tqdm

from pyrewton.cazymes.annotate_cazome.tools.parse_tools import (
    parse_dbcan_output,
    parse_ecami_output, 
    parse_cupp_output
)
from pyrewton.sql.sql_orm import get_cazome_db_connection
from pyrewton.sql.sql_interface import add_data
from pyrewton.sql.sql_interface.add_data import add_prot_tax_genomes, add_cazy_data
from pyrewton.sql.sql_interface.load_data import load_genbank_data


def main(args: Optional[List[str]] = None, logger: Optional[logging.Logger] = None):
    """Coordinate the retrieval of protein annotations from GenBank (.gbff) files.
    Including building parser, logger and output directory.
    Return dataframe of protein data.
    """
    if logger is None:
        config_logger(args)
    logger = logging.getLogger(__package__)

    time_stamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    
    if (args.new_db is None) and (args.db is None):
        db_path = Path(f"cazome_db_{time_stamp}.db")
        logger.warning(
            "Neither --new_db or --db flags used.\n"
            f"Building a new database at {db_path}"
        )

    elif (args.new_db is not None) and (args.db is None):  # make a new database, where specified
        logger.warning(f"Building a new CAZome database at {args.new_db}")
        db_path = args.new_db

    elif (args.new_db is None) and (args.db is not None):  # make a new database, where specified
        if args.db.exists() is False:
            logger.error(f"Could not find database at {args.db}\nCheck the path is correct.\nTermianting program.")
            sys.exit(1)

        logger.warning(f"Adding data to an existing CAZome database at {args.new_db}")
        db_path = args.db

    else:  # both are not none
        logger.error(
            "Specified to create a new CAZome database (using the --db flag)\n"
            "AND to add data to an existing CAZome database (using the --db flag).\n"
            "One or the other can be done, not both\n"
            "Terminating program"
        )
        sys.exit(1)

    if str(db_path.parent) != '.':  # make output directory
        make_output_directory(db_path.parent, args.force, args.nodelete)

    connection = get_cazome_db_connection(db_path, args)

    # parse protein FASTA files and add new proteins, taxonomies and genomic
    # assemblies to the db
    if args.protein_dir is not None:
        add_prot_tax_genomes.add_protein_data(connection, args, logger)

    # add annotations from CAZy that are stored in a local CAZyme database
    # built using cazy_webscraper
    if args.cazy is not None:
        add_cazy_annotations(connection, args, logger)
    
    # add data from cazyme classifiers specified in the config file
    if args.config_file is not None:
        config_data = load_config_file(args, logger) 
        # {toolname: 
        # {'version': '0.0.0',
        # 'cazy_release': 'date',
        # 'dir': 'path to parent directory containing the output directories with the tool output'}}
    
        # add classifiers to the database
        add_prediction_classifiers(config_data, logger, connection)

        # add prediction data


def gather_domain_annotations(config_data, logger, connection):
    """Gather all predicted CAZome domain annotations. 
    One tool + one family = a new domain
    
    :param config_data: dict {toolname: {version: str, cazy_release: str, dir: str}}
    :param logger: logger object
    :param connection: open connection to sql db
    """
    # existing_db_proteins = {genbank_accession: db_id}
    existing_db_proteins = load_genbank_data.get_protein_db_ids(connection)

    # detect output dirs for all classifiers
    # {classifier name: [path to output dirs]}
    output_dirs = get_classifier_paths(config_data, logger)

    if len(list(output_dirs.key())) == 0:
        return
    
    # get all cazyme annotations
    # {prot acc: {toolname: {family : [domain ranges]}}
    # # only contains classifiers with ids in the local db
    tool_outputs = get_all_cazyme_annotations(config_data, output_dirs, logger, connection)

    # add all new CAZy families
    all_cazy_families = set()
    for protein_acc in tool_outputs:
        for classifier in tool_outputs[protein_acc]:
            for family in tool_outputs[protein_acc][classifier]:
                all_cazy_families.add(family)

    add_cazy_data.add_families(all_cazy_families, connection, logger)

    # {classifier : db id}
    db_classifiers = load_genbank_data.get_classifier_db_ids(connection)

    # load existing families = {family: db_id}
    family_table_dict = load_genbank_data.get_family_db_ids(connection)

    # {protein_id : [ (classifier_id, family_id) ]}
    existing_domains = load_genbank_data.get_protein_annotations_table(connection)

    new_domains = []  # [ (protein_id, classifier_id, family_id) ]
    for protein_acc in tqdm(tool_outputs, desc="Parsing proteins to identify new CAZyme domains"):
        try:
            protein_id = existing_db_proteins[protein_acc]
        except KeyError:
            logger.warning(
                f"No record for protein {protein_acc} in local db.\n"
                f"Not adding domains for {protein_acc} to db"
            )
            continue

        for classifier in tool_outputs[protein_acc]:
            classifier_id = db_classifiers[classifier]
        
            for cazy_family in tool_outputs[protein_acc][classifier]:
                try:
                    family_id = family_table_dict[cazy_family]
                except KeyError:
                    logger.warning(
                        f"No record for CAZy family {cazy_family} in local db.\n"
                        f"Not adding domains for {cazy_family} to db"
                    )
                    continue

                for domain_range in tool_outputs[protein_acc][classifier][family]:
                    # test if domain already represented 
                    try:
                        if (classifier_id, family_id) not in existing_domains[protein_id]:
                            new_domains.append( (protein_id, classifier_id, family_id, domain_range) )

                    except KeyError:
                        new_domains.append( (protein_id, classifier_id, family_id, domain_range) )

    if len(new_domains) > 0:
        logger.warning(f"Adding {len(new_domains)} CAZyme domains to the Domains table")
        add_data.insert_data(connection, "Domains", ['protein_id','classifier_id','family_id','domain_range'], new_domains)


def get_all_cazyme_annotations(config_data, output_dirs, logger, connection):
    """Retrieve all CAZyme annotations from all classifiers
    Return dict {prot acc: {toolname: {family : [domain ranges]}}
    """
    tool_outputs = {}  # {genomic acc: {prot acc: {toolname: {family : [domain ranges]}}}
    
    # new_domains_dict = {}  # {ncbi acc: {prot_id: local db prot id, classifier-name: {classifier_id, db id, fam: str, fam_id: db id, domain_range:str}}}

    for classifier in config_data:
        # check we have paths to the output dirs
        try:
            if len(output_dirs[classifier]) == 0:
                logger.error(
                    f"Could not find any output directories for {classifier}\n"
                    "Check path provided to parent output dir in the config file is correct\n"
                    f"Not adding CAZyme domains from {classifier}"
                )
                continue
        except KeyError:
            logger.error(
                f"Could not find any output directories for {classifier}\n"
                "Check path provided to parent output dir in the config file is correct\n"
                f"Not adding CAZyme domains from {classifier}"
            )
            continue

        # {classifier : db id}
        db_classifiers = load_genbank_data.get_classifier_db_ids(connection)

        try:
            db_classifiers[classifier]
        except KeyError:
            logger.error(
                f"Could not find record that exactly matched {classifier} in the local db\n"
                f"Not adding CAZyme domains from {classifier}"
            )
            continue

        if classifier.lower().startswith('dbcan'):
            try:
                dbcan_version = int(config_data[classifier]['version'].split(".")[0])
            except KeyError:
                dbcan_version = None

            for dbcan_output_dir in tqdm(output_dirs[classifier], desc="Parsing dbCAN output files"):
                overview_file_path = Path(dbcan_output_dir) / "overview.txt"
                dbcan_output = parse_dbcan_output(overview_file_path, logger, dbcan_version)
                # {toolname: {protein acc: {family: domain range}}}
                if dbcan_output is not None:
                    tool_outputs = add_tool_outputs(tool_outputs, dbcan_output)

        elif classifier.lower().startswith('cupp'):
            for cupp_output_dir in tqdm(output_dirs[classifier], desc="Parsing CUPP output files"):
                cupp_out_file = Path(cupp_output_dir) / "cupp_output.fasta.log"
                cupp_output = parse_cupp_output(cupp_out_file, logger)
                if cupp_output is not None:
                    tool_outputs = add_tool_outputs(tool_outputs, cupp_output)

        elif classifier.lower().startswith('ecami'):
            for ecami_output_dir in tqdm(output_dirs[classifier], desc="Parsing eCAMI output files"):
                ecami_out_file = Path(ecami_output_dir) / "ecami_output.txt"
                ecami_output = parse_ecami_output(ecami_out_file, logger)
                if ecami_output is not None:
                    tool_outputs = add_tool_outputs(tool_outputs, ecami_output)

        else:
            logger.error(
                f"Pyrewton is not currently setup to parse output from {classifier}\n"
                "Please raise this as an issue at GitHub:\n"
                "https://github.com/HobnobMancer/pyrewton/issues\n"
                f"Not adding CAZyme domains from {classifier}"
            )

    return tool_outputs



    # load data from local CAZome db

    # add news proteins

    # load data from local CAZome db


    # get paths to directories containing dbCAN output
    dbcan_output_dirs = get_paths.get_dir_paths(args.dbcan_dir)

    if len(dbcan_output_dirs) == 0:
        logger.error(
            f"No dirs retrieved from {args.dbcan_dir}\nCheck the correct dir was provided\n"
            "The output for each genome must be stored within a separate dir\n"
            "within a parent dir, which is provided to pyrewton.\n"
            "Terminating program"
        )
        sys.exit(1)
        



    # retrieve protein data from the FASTA file parsed by dbCAN
    # tax_dict = {genomic_accession: {'genus': str, 'species': str, 'tax_id': str}}
    # protein_dict = {genomic_accession: {protein_accession: str(sequence)}}
    protein_dict, tax_dict = get_protein_data(protein_fasta_files)

    # retrieve CAZymes predicted by dbCAN
    cazome_dict, domain_range_dict = get_dbcan_annotations(dbcan_output_dirs, protein_dict)

    if args.cazy is not None:
        cazome_dict = add_cazy_annotations(cazome_dict, cazy_db_connection)

    add_data_to_db(cazome_dict, tax_dict, domain_range_dict, connection, args)

    # cache the dict into the output dir
    cache_dict(cazome_dict, time_stamp, args)


def load_config_file(args, logger):
    try:
        with open(args.config_file) as fh:
            return yaml.full_load(fh)
    except FileNotFoundError:
        logger.error(
            f"Could not found config file at {args.config_file}\n"
            "Check the path is correct\n"
            "Terminating program"
        )
        sys.exit(1)


def add_prediction_classifiers(config_data, logger, connection):
    """Add classifiers that predict CAZyme annotations to db
    
    :param config_data: dict {toolname: {version: str, cazy_release: str, dir: str}}
    :param logger: logger object
    :param connection: open connection to SQL db
    """
    for classifier in config_data:
        try:
            c_version = config_data[classifier]['version']
        except KeyError:
            logger.warning(f"No version number found for {classifier}. Setting to null")
            c_version = None
        try:
            cazy_release = config_data[classifier]['cazy_release']
        except KeyError:
            logger.warning(f"No training dataset release or version number found for {classifier}. Setting to null")
            cazy_release = None

        add_cazy_data.add_classifiers(
            'classifier',
            connection, 
            logger,
            classifier_version=c_version,
            classifier_training_date=cazy_release,
        )

        if classifier.lower().startswith('dbcan'):
            try:
                c_version = config_data[classifier]['version']
                if str(c_version).startswith('2'):
                    kmer_classifier = 'dbCAN-sub'
                elif str(c_version).startswith('3'):
                    kmer_classifier = 'eCAMI'
                else:
                    kmer_classifier = 'dbCAN-sub'
            except KeyError:
                logger.warning(
                    f"No version number found for {classifier}\n"
                    "Presuming using dbCAN version 4, so adding classifiers HMMER, DIAMOND and dbCAN-sub to db"
                )
                kmer_classifier = 'dbCAN-sub'
                c_version=None

            add_cazy_data.add_classifiers(
                'HMMER',
                connection, 
                logger,
                classifier_version=f"dbCAN=={c_version}",
                classifier_training_date=cazy_release,
            )
            add_cazy_data.add_classifiers(
                'DIAMOND',
                connection, 
                logger,
                classifier_version=f"dbCAN=={c_version}",
                classifier_training_date=cazy_release,
            )
            add_cazy_data.add_classifiers(
                kmer_classifier,
                connection, 
                logger,
                classifier_version=f"dbCAN=={c_version}",
                classifier_training_date=cazy_release,
            )

            
def get_classifier_paths(config_data, logger):
    """Get paths to all output directories for all classifiers in config_data
    
    :param config_data: dict {toolname: {version: str, cazy_release: str, dir: str}}
    :param logger: logger object

    Return dict {classifier name: [path to output dirs]}
    """
    output_dirs = {}  # {toolname: [paths to output dirs]}
    for classifier in config_data:
        try:
            classifier_parent_dir = config_data[classifier]['dir']
        except KeyError:
            logger.error(
                f"Could not find dir specified for {classifier}."
                f"Not adding data for {classifier} to the local CAZome db"
            )
            continue

        classifier_paths = get_paths.get_dir_paths(classifier_parent_dir)

        if len(classifier_paths) == 0:
            logger.error(
                f"No dirs retrieved fpr {classifier} from {classifier_parent_dir}\n"
                "Check the correct path was provided.\n"
                f"Not adding data for {classifier} to the local CAZome db"
            )
            continue

        output_dirs[classifier] = classifier_paths

    return output_dirs


def add_tool_outputs(tool_outputs, tool_dict):
    """Add output from tool (tool_dict) to collection of all tool outputs (tool_outputs)"""
    for prot_acc in tool_dict:
        try:  # check if parsed protein before
            tool_outputs[prot_acc]
        except KeyError:
            tool_outputs[prot_acc] = {}

        for toolname in tool_dict[prot_acc]:
            try:  # check if parsed tool before
                tool_outputs[prot_acc][toolname]
            except KeyError:
                tool_outputs[prot_acc][toolname] = {}

            for family in tool_dict[prot_acc][toolname]:
                for domain_range in tool_dict[prot_acc][toolname][family]:

                    try:  # check if parsed this domain before
                        tool_outputs[prot_acc][toolname][family].append(domain_range)

                    except KeyError:
                        tool_outputs[prot_acc][toolname][family] = [domain_range]

    return tool_outputs


def add_cazy_annotations(connection, args, logger):
    """Add CAZy family annotations from CAZy to the data dict
    
    :param connection: open sqlalchemy connection to an SQLite3 db engine
    :param args: cli args parser
    :param logger: logger object
    """
    # existing_db_proteins = {genbank_accession: db_id}
    existing_db_proteins = load_genbank_data.get_protein_db_ids(connection)

    # get connection to local cazyme database
    cazy_db_connection = get_db_connection(args.cazy, args, False)

    # load cazy annotations into dict {protein acc: {families}}
    cazy_dict = {}
    cazy_families = set()
    with Session(bind=cazy_db_connection) as session:
        db_query = session.query(Genbank, CazyFamily).\
            join(CazyFamily, Genbank.families).\
            all()
        
    for record in tqdm(db_query, desc="Getting all CAZy family annotations from the CAZy database"):
        prot_accession = record[0].genbank_accession
        if prot_accession not in existing_db_proteins.keys():
            continue

        family = record[1].family

        try:
            cazy_dict[prot_accession].add(family)
        except KeyError:
            cazy_dict[prot_accession] = {family}

        cazy_families.add(family)

    # test if need to add new fams
    add_cazy_data.add_families(cazy_families, connection, logger)

    # add cazy as a classifier
    add_cazy_data.add_classifiers('CAZy', connection, logger, classifier_version=args.cazy_date)
    db_classifiers = load_genbank_data.get_classifier_db_ids(connection) # {classifier : db id}
    cazy_classifier_id = db_classifiers['CAZy']

    add_cazy_data.add_classifications(cazy_dict, cazy_classifier_id, connection, logger)





# def cache_dict(cazome_dict, time_stamp, args):
#     """Cache the dict of CAZome protein data
    
#     :param cazome_dict: {genomic_accession: {protein_accession: {tool: fams}}}
#     :param time_stamp: str, date and time program was invoked
#     :param args: cmd-line args parser
    
#     Return nothing"""
#     cache_path = Path(f"dbcan_predictions_{time_stamp}.json")

#     cahce_dict = cazome_dict

#     if args.output_dir is not None:
#         cache_path = args.output_dir / cache_path

#     json_data = convert_for_serialisation(cahce_dict, args)

#     with open(cache_path, 'w') as fh:
#         json.dump(json_data, fh)

#     return
    

# def convert_for_serialisation(protein_dict, args):
#     """Convert all data types in the dict to those suitable for JSON serialisation."""
#     for genomic_accession in protein_dict:
#         for protein_accession in protein_dict[genomic_accession]:
#             if args.cazy is not None:
#                 tools = ['hmmer', 'hotpep', 'diamond', 'dbcan', 'cazy', '#ofTools']
#             else:
#                 tools = ['hmmer', 'hotpep', 'diamond', 'dbcan', '#ofTools']

#             for tool in tools:
#                 try:
#                     if len(protein_dict[genomic_accession][protein_accession][tool]) == 0:
#                         protein_dict[genomic_accession][protein_accession][tool] = None
#                     else:
#                         protein_dict[genomic_accession][protein_accession][tool] = (
#                             str(protein_dict[genomic_accession][protein_accession][tool])
#                         )
#                 except TypeError:
#                     pass  # raised when '#ofTools' is an int

#     return protein_dict


# def add_data_to_db(cazome_dict, tax_dict, domain_dict, connection, args):
#     """Add data to the database
    
#     :param cazome_dict: {genomic_accession: {protein_accession: {tool: fams}}}
#     :param tax_dict: {genomic_accession: {'genus': str, 'species': str, 'tax_id': str}}
#     :param domain_dict: dict of domain ranges {protein: {tool: {fam: set()}}}
#     :param connection: open sqlalchemy connection to an SQLite3 db
#     :param args: cmd-line args parser
    
#     Return nothing"""
#     # add taxonomy (species) data to the db
#     add_genbank_data.add_species_data(tax_dict, connection)

#     # retrieve Taxonomies table
#     tax_table_dict = load_genbank_data.get_tax_table(connection)

#     # add genomic assemblies to the database
#     add_genbank_data.add_genomic_accessions(tax_dict, tax_table_dict, connection)

#     # retrive genomic assembly records from the db
#     assembly_dict = load_genbank_data.get_assemblies_table(connection)

#     # add proteins
#     add_genbank_data.add_proteins(cazome_dict, assembly_dict, connection)

#     # add classifiers
#     classifier_data = [
#         ('dbCAN', '2.0.11', 'March 2020'),
#         ('HMMER', 'dbCAN=2.0.11', 'March 2020'),
#         ('Hotpep', 'dbCAN=2.0.11', 'March 2020'),
#         ('DIAMOND', 'dbCAN=2.0.11', 'March 2020'),
#     ]

#     if args.cazy_date is not None:
#         classifier_data.append( ('CAZy', None, args.cazy_date) )
#     else:
#         classifier_data.append( ('CAZy', None, None) )

#     add_data.insert_data(
#         connection,
#         'Classifiers',
#         ['classifier', 'version', 'cazy_training_set_date'],
#         classifier_data,
#     )

#     # add CAZy families
#     add_genbank_data.add_families(cazome_dict, connection)

#     # load Proteins, CazyFamilies and Classifiers tables into dicts
#     protein_db_dict = load_genbank_data.get_protein_db_ids(connection)
#     classifer_db_dict = load_genbank_data.get_classifier_db_ids(connection)
#     family_db_dict = load_genbank_data.get_family_db_ids(connection)

#     # add CAZy and dbCAN domain annotations
#     add_genbank_data.add_classifications(
#         cazome_dict,
#         domain_dict,
#         protein_db_dict,
#         classifer_db_dict,
#         family_db_dict,
#         connection,
#         args,
#     )

#     return


if __name__ == "__main__":
    main()
