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
import os
import re
import sys

import pandas as pd

from datetime import datetime
from pathlib import Path
from typing import List, Optional

from Bio import SeqIO
from cazy_webscraper.sql.sql_orm import (
  CazyFamily,
  Genbank,
  Session,
)
from saintBioutils.utilities import file_io
from saintBioutils.utilities import logger
from saintBioutils.utilities.file_io import get_paths
from saintBioutils.utilities.logger import config_logger

from tqdm import tqdm

from pyrewton.sql.sql_orm import get_cazome_db_connection
from pyrewton.sql.sql_interface import add_data, load_data
from pyrewton.utilities.parsers.cmd_parser_compile_db import build_parser


def main(argv: Optional[List[str]] = None, logger: Optional[logging.Logger] = None):
    """Coordinate the retrieval of protein annotations from GenBank (.gbff) files.
    Including building parser, logger and output directory.
    Return dataframe of protein data.
    """
    time_stamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")

    if argv is None:
        parser = build_parser()
        args = parser.parse_args()
    else:
        parser = build_parser(argv)
        args = parser.parse_args()

    if logger is None:
        config_logger(args)
    logger = logging.getLogger(__package__)

    # compile path to the output CAZome database
    if args.output_db is None:
        db_path = Path(f"cazome_db_{time_stamp}.db")
    else:
        db_path = args.output_db

    if args.output_dir is not None:
        make_output_directory(args.output_dir, args.force, args.nodelete)
        db_path = args.output_dir / db_path

    if db_path.exists() and args.force is False:
            logger.error(
                f"Db at {db_path} exists\nForce is false\nNot overwriting file\n"
                "To overwrite the file use -f or --force\nTerminating program"
            )
            sys.exit(1)
    
    if args.cazy is not None:
        if args.cazy.exists() is False:
            logger.error(
                f"Path to local CAZyme db ({args.cazy})\ndoes not exist\n"
                "Check the correct path was provided\n"
                "Terminating program"
            )
            sys.exit(1)
        cazy_db_connection = get_db_connection(args.cazy, args, False)

    connection = get_cazome_db_connection(db_path, args)

    # get paths to directories containing dbCAN output
    dbcan_output_dirs = get_dir_paths(args.dbcan_dir)

    if len(dbcan_output_dirs) == 0:
        logger.error(
            f"No dirs retrieved from {args.dbcan_dir}\nCheck the correct dir was provided\n"
            "The output for each genome must be stored within a separate dir\n"
            "within a parent dir, which is provided to pyrewton.\n"
            "Terminating program"
        )
        sys.exit(1)

    # get paths to FASTA files of protein sequences parsed by dbCAN
    protein_fasta_files = get_file_paths(args.protein_dir, suffixes='fasta')

    if len(protein_fasta_files) == 0:
        logger.error(
            f"No FASTA files containing protein sequences extracted from the genomes retrieved from\n"
            f"{args.protein_dir}\nCheck the correct dir was provided\n"
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


def get_protein_data(protein_fasta_files):
    """Retrieve protein, genomic and taxonomic data from the FASTA files parsed by dbCAN.
    
    :param protein_fasta_files: list of Paths, one path per FASTA file
    
    Return two dicts
    tax_dict = {genomic_accession: {'genus': str, 'species': str, 'tax_id': str}}
    protein_dict = {genomic_accession: {protein_accession: str(sequence)}}
    """
    logger = logging.getLogger(__name__)

    protein_dict = {}
    tax_dict = {}  # {genomes: {'genus': str, 'species': str, 'tax_id': str}}

    for fasta in tqdm(protein_fasta_files, desc="Retrieving protein data from FASTA files"):
        for record in SeqIO.parse(fasta, 'fasta'):
            protein_acc = record.id

            desc = record.description.split(" ")

            assembly_acc = desc[2]
            assembly_acc = assembly_acc.replace("_", ".")
            assembly_acc = assembly_acc.replace("GCA.", "GCA_")
            assembly_acc = assembly_acc.replace("GCF.", "GCF_")

            tax_id = desc[3].replace("txid", "")
            
            genus = desc[4]
            species = desc[5]

            # add taxonomy data
            try:
                tax_dict[assembly_acc]
            
            except KeyError:
                tax_dict[assembly_acc] = {
                    'genus': genus,
                    'species': species,
                    'txid': tax_id,
                }

            # add protein data
            try:
                protein_dict[assembly_acc]

                try:
                    protein_dict[assembly_acc][protein_acc]

                    if protein_dict[assembly_acc][protein_acc] != str(record.seq):
                        logger.warning(
                            f"Multiple unique sequences retrieved for {protein_acc}\n"
                            "Retrieving the first protein sequence only"
                        )

                except KeyError:
                    protein_dict[assembly_acc][protein_acc] = str(record.seq)

            except KeyError:
                protein_dict[assembly_acc] = {
                    protein_acc: str(record.seq),
                }

    return protein_dict, tax_dict


def get_dbcan_annotations(dbcan_output_dirs, protein_dict):
    """Retrieve CAZy family annotations from dbCAN ouput.
    
    :param dbcan_output_dirs: list of Paths() to output dirs containing dbCAN output
    :param protein_dict: {genomic_accession: {protein_accession: {sequence}}
    
    Return 
    cazome_dict: {genomic_accession: {
        protein_accession: {
            'sequence': SeqIO.seq, 'diamond':set(), 'hmmer':set(), 'hotpep':set(), '#ofTools': int, 'dbcan':set()}}
    domain_dict: dict of domain ranges {protein: {tool: {fam: set()}}}
    """
    logger = logging.getLogger(__name__)

    # {genomic_accession: {
    # protein_accession: {
    # 'sequence': SeqIO.seq, 'diamond':set(), 'hmmer':set(), 'hotpep':set(), '#ofTools': int, 'dbcan':set()}}
    cazyme_dict = {}

    # {protein: {tool/classifier: {fam: set}}}
    domain_dict = {}

    for dbcan_dir in tqdm(dbcan_output_dirs, "Parsing dbCAN output"):
        # dir name format: GCA_123456789_1
        genomic_accession = f"{(dbcan_dir.name).split('_')[0]}_{(dbcan_dir.name).split('_')[1]}.{(dbcan_dir.name).split('_')[2]}"

        try:
            protein_dict[genomic_accession]
        
        except KeyError:
            logger.error(
                f"Genome {genomic_accession} parsed by dbCAN but not FASTA file of protein seqs\n"
                "found for this assembly\n"
                "Not parsing protein data for this genome"
            )
            continue

        overview_path = dbcan_dir / "overview.txt"

        try:
            with open(overview_path, 'r') as ofh:
                overview_file = ofh.read().splitlines()
        except FileNotFoundError:
            logger.error(
                f"Could not find dbCAN output file: {overview_file}\n"
                f"for genomic assembly {genomic_accession}"
                "Not retrieving CAZymes from this genome"
            )
            continue
                
        for line in overview_file[1:]:  # first line contains headings
            # separate the line content
            line = line.split("\t")

            protein_accession = line[0]

            if protein_accession.find("|") != -1:
                protein_accession = protein_accession.split("|")[1]

            try:
                protein_dict[genomic_accession][protein_accession]
            except KeyError:
                logger.warning(
                    f"{protein_accession} protein retrieved from dbCAN for {genomic_accession}\n"
                    "but not protein data retrieved from the FASTA file of protein seqs for this assembly\n"
                    "Not adding protein to the database"
                )
                continue

            # create a CazymeProteinPrediction instance each for HMMER, Hotpep and DIAMOND
            hmmer_predictions, domain_ranges = get_hmmer_prediction(line[1], protein_accession)
            hmmer_predictions = list(hmmer_predictions)
            hotpep_predictions = list(get_hotpep_prediction(line[2], protein_accession))
            diamond_predictions = list(get_diamond_prediction(line[3], protein_accession))

            # get the dbCAN consensus result
            no_of_tools = line[-1]
            
            hmmer_hotpep = list(set(hmmer_predictions) & set(hotpep_predictions))
            hmmer_diamond = list(set(hmmer_predictions) & set(diamond_predictions))
            hotpep_diamond = list(set(hotpep_predictions) & set(diamond_predictions))

            dbcan_predictions = list(set(hmmer_hotpep + hmmer_diamond + hotpep_diamond))
            
            tools = [
                ('hmmer', hmmer_predictions),
                ('hotpep', hotpep_predictions),
                ('diamond', diamond_predictions),
                ('dbcan', dbcan_predictions),
            ]

            for tool_data in tools:
                try:
                    cazyme_dict[genomic_accession]
                    try:
                        cazyme_dict[genomic_accession][protein_accession]
                        try:
                            existing = cazyme_dict[genomic_accession][protein_accession][tool_data[0]]
                            updated = existing.union(set(tool_data[1]))
                            cazyme_dict[genomic_accession][protein_accession][tool_data[0]] = updated
                        
                        except KeyError:  # tool raised KeyError
                            cazyme_dict[genomic_accession][protein_accession][tool_data[0]] = set(tool_data[1])

                    except KeyError:  # protein accession raised KeyError
                        cazyme_dict[genomic_accession][protein_accession] = {
                            tool_data[0]: set(tool_data[1]),
                            'sequence': protein_dict[genomic_accession][protein_accession],
                        }

                except KeyError:  # genomic accession raised KeyError
                    cazyme_dict[genomic_accession] = {
                        protein_accession: {
                            tool_data[0]: set(tool_data[1]),
                            'sequence': protein_dict[genomic_accession][protein_accession],
                        }
                    }

            cazyme_dict[genomic_accession][protein_accession]['#ofTools'] = str(no_of_tools)
            cazyme_dict[genomic_accession][protein_accession]['sequence'] = protein_dict[genomic_accession][protein_accession]

            for fam in domain_ranges:
                ranges = domain_ranges[fam]

                for drange in ranges:
                    try:
                        domain_dict[protein_accession]
                        try:
                            domain_dict[protein_accession][fam].add(drange)
                        except KeyError:
                            domain_dict[protein_accession][fam] = {drange}
                    except KeyError:
                        domain_dict[protein_accession] = {fam: {drange}}

    return cazyme_dict, domain_dict


def get_hmmer_prediction(hmmer_data, protein_accession):
    """Retrieve HMMER prediciton data from the dbCAN overview.txt file.

    :param hmmer_data: str, output data from HMMER written in the overview.txt file.
    :param protein_accession: str

    Return set of CAZy family predictions, dict of domain ranges {fam: set(ranges)}
    """
    logger = logging.getLogger(__name__)

    if hmmer_data == "-":
        return set(), {}
    
    cazy_fams = set()
    ranges = {}

    predicted_domains = hmmer_data.split("+")
    for domain in predicted_domains:
        # separate the predicted CAZy family from the domain range
        domain_name = domain.split("(")[0]  # CAZy (sub)family
        domain_range = domain.split("(")[1]
        domain_range = domain_range.replace(")", "")

        formated_name = None

        if domain_name.find("_") != -1:
            try: 
                re.match(r"\D{2,3}\d+?_\D", domain_name).group()  # check unusal CAZy family formating
                formated_name = domain_name.split("_")[0]

            except AttributeError:  # raised if not an usual CAZy family format
                try:
                    re.match(r"\D{2,3}\d+?_\d+", domain_name).group()  # check if a subfamily
                    formated_name = domain_name.split("_")[0]
                    
                except AttributeError:
                    logger.warning(
                        f"Unknown data type of {domain_name} for protein {protein_accession}, for HMMER.\n"
                        "Not adding as domain to CAZy family annotations for the genome"
                    )
        
        else:
            try:
                re.match(r"\D{2,3}\d+", domain_name).group()
                formated_name = domain_name
            except AttributeError:  # raised if doesn't match expected CAZy family
                logger.warning(
                    f"Unknown data type of {domain_name} for protein {protein_accession}, for HMMER.\n"
                    "Not adding as domain to CAZy family annotations for the genome"
                )

        if formated_name is not None:
            try:
                ranges[formated_name].add(domain_range)
            except KeyError:
                ranges[formated_name] = {domain_range}

        cazy_fams.add(formated_name)
    
    return cazy_fams, ranges


def get_hotpep_prediction(hotpep_data, protein_accession):
    """Retrieve Hotpep prediciton data from the dbCAN overview.txt file.

    :param hotpep_data: str, output data from HMMER written in the overview.txt file.
    :param protein_accession: str

    Return set of CAZy family predictions
    """
    logger = logging.getLogger(__name__)

    if hotpep_data == "-":
        return set()
    
    cazy_fams = set()

    predicted_domains = hotpep_data.split("+")
    for domain in predicted_domains:
        # remove k-mer group number
        domain = domain.split("(")[0]

        # check if a CAZy subfamily was predicted
        if domain.find("_") != -1:
            try:
                # check if its a CAZy subfamily
                re.match(r"\D{2,3}\d+?_\d+", domain).group()
                cazy_fams.add(domain[:domain.find("_")])

            except AttributeError:  # raised if it doesn't match the expected CAZy subfamily format
                logger.warning(
                    f"Unknown data type of {domain} for protein {protein_accession}, for Hotpep.\n"
                    "Not adding as domain to CAZy family annotations for the genome"
                )
        
        else:
            try:  # check it is a CAZy family
                re.match(r"\D{2,3}\d+", domain).group()
                cazy_fams.add(domain)

            except AttributeError:  # raised if doesn't match expected CAZy family
                logger.warning(
                    f"Unknown data type of {domain} for protein {protein_accession}, for Hotpep.\n"
                    "Not adding as domain to CAZy family annotations for the genome"
                )

    return cazy_fams


def get_diamond_prediction(diamond_data, protein_accession):
    """Retrieve DIAMOND prediciton data from the dbCAN overview.txt file.

    :param diamond_data: str, output data from HMMER written in the overview.txt file.
    :param protein_accession: str

    Return set of CAZy family predictions
    """
    logger = logging.getLogger(__name__)

    if diamond_data == "-":
        return set()
    
    cazy_fams = set()

    predicted_domains = diamond_data.split("+")
    for domain in predicted_domains:
        # check if a CAZy subfamily
        try:
            re.match(r"\D{2,3}\d+?_\d+", domain).group()
            cazy_fams.add(domain.split("_")[0])
        
        except AttributeError:
            try:
                # check if a cazy family
                re.match(r"\D{2,3}\d+", domain).group()
                cazy_fams.add(domain)
            
            except AttributeError:
                try:
                    re.match(r"\d+?\.(\d+?|\*)\.(\d+?|\*)\.(\d+?|\*)", domain). group()
                except AttributeError:
                    try:
                        re.match(r"\d+?\.(\d+?|\*)\.(\d+?|\*)\.-", domain). group()
                    except AttributeError:
                        logger.warning(
                            f"Unexpected data format of '{domain}' for protein "
                            f"{protein_accession}, for DIAMOND.\n"
                            "Not adding as domain to CAZy family annotations for the genome"
                        )

    return cazy_fams


def add_cazy_annotations(cazome_dict, connection):
    """Add CAZy family annotations from CAZy to the data dict
    
    :param cazome_dict: 
    {genomic_accession: {protein_accession: {sequence: str, tool: set, #ofTools: int}}}
    :param connection: open sqlalchemy connection to an SQLite3 db engine
    
    Return cazome_dict
    """
    genomic_accessions = list(cazome_dict.keys())
    for genomic_accession in tqdm(genomic_accessions, desc="Adding CAZy data per genome"):
        protein_accessions = list(cazome_dict[genomic_accession].keys())

        for protein_accession in tqdm(protein_accessions, desc="Retrieving CAZy annotations"):
            with Session(bind=connection) as session:
                db_query = session.query(Genbank, CazyFamily).\
                    join(CazyFamily, Genbank.families).\
                    filter(Genbank.genbank_accession==protein_accession).\
                    all()

            if len(db_query) == 0:
                cazy_data = None

            else:
                cazy_data = []
                for obj in db_query:
                    cazy_data.append(obj[1].family)
                cazy_data.sort()
            
            cazome_dict[genomic_accession][protein_accession]['cazy'] = cazy_data
    
    return cazome_dict


def cache_dict(cazome_dict, time_stamp, args):
    """Cache the dict of CAZome protein data
    
    :param cazome_dict: {genomic_accession: {protein_accession: {tool: fams}}}
    :param time_stamp: str, date and time program was invoked
    :param args: cmd-line args parser
    
    Return nothing"""
    cache_path = Path(f"dbcan_predictions_{time_stamp}.json")

    cahce_dict = cazome_dict

    if args.output_dir is not None:
        cache_path = args.output_dir / cache_path

    json_data = convert_for_serialisation(cahce_dict, args)

    with open(cache_path, 'w') as fh:
        json.dump(json_data, fh)

    return
    

def convert_for_serialisation(protein_dict, args):
    """Convert all data types in the dict to those suitable for JSON serialisation."""
    for genomic_accession in protein_dict:
        for protein_accession in protein_dict[genomic_accession]:
            if args.cazy is not None:
                tools = ['hmmer', 'hotpep', 'diamond', 'dbcan', 'cazy', '#ofTools']
            else:
                tools = ['hmmer', 'hotpep', 'diamond', 'dbcan', '#ofTools']

            for tool in tools:
                try:
                    if len(protein_dict[genomic_accession][protein_accession][tool]) == 0:
                        protein_dict[genomic_accession][protein_accession][tool] = None
                    else:
                        protein_dict[genomic_accession][protein_accession][tool] = (
                            str(protein_dict[genomic_accession][protein_accession][tool])
                        )
                except TypeError:
                    pass  # raised when '#ofTools' is an int

    return protein_dict


def add_data_to_db(cazome_dict, tax_dict, domain_dict, connection, args):
    """Add data to the database
    
    :param cazome_dict: {genomic_accession: {protein_accession: {tool: fams}}}
    :param tax_dict: {genomic_accession: {'genus': str, 'species': str, 'tax_id': str}}
    :param domain_dict: dict of domain ranges {protein: {tool: {fam: set()}}}
    :param connection: open sqlalchemy connection to an SQLite3 db
    :param args: cmd-line args parser
    
    Return nothing"""
    # add taxonomy (species) data to the db
    add_data.add_species_data(tax_dict, connection)

    # retrieve Taxonomies table
    tax_table_dict = load_data.get_tax_table(connection)

    # add genomic assemblies to the database
    add_data.add_genomic_accessions(tax_dict, tax_table_dict, connection)

    # retrive genomic assembly records from the db
    assembly_dict = load_data.get_assemblies_table(connection)

    # add proteins
    add_data.add_proteins(cazome_dict, assembly_dict, connection)

    # add classifiers
    classifier_data = [
        ('dbCAN', '2.0.11', 'March 2020'),
        ('HMMER', 'dbCAN=2.0.11', 'March 2020'),
        ('Hotpep', 'dbCAN=2.0.11', 'March 2020'),
        ('DIAMOND', 'dbCAN=2.0.11', 'March 2020'),
    ]

    if args.cazy_date is not None:
        classifier_data.append( ('CAZy', None, args.cazy_date) )
    else:
        classifier_data.append( ('CAZy', None, None) )

    add_data.insert_data(
        connection,
        'Classifiers',
        ['classifier', 'version', 'cazy_training_set_date'],
        classifier_data,
    )

    # add CAZy families
    add_data.add_families(cazome_dict, connection)

    # load Proteins, CazyFamilies and Classifiers tables into dicts
    protein_db_dict = load_data.get_protein_db_ids(connection)
    classifer_db_dict = load_data.get_classifier_db_ids(connection)
    family_db_dict = load_data.get_family_db_ids(connection)

    # add CAZy and dbCAN domain annotations
    add_data.add_classifications(
        cazome_dict,
        domain_dict,
        protein_db_dict,
        classifer_db_dict,
        family_db_dict,
        connection,
        args,
    )

    return


if __name__ == "__main__":
    main()
