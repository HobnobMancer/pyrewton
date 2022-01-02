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

import pandas as pd

from datetime import datetime
from pathlib import Path
from typing import List, Optional

from Bio import SeqIO
from saintBioutils.utilities import file_io
from saintBioutils.utilities import logger
from saintBioutils.utilities.file_io import get_paths
from saintBioutils.utilities.logger import config_logger
from tqdm import tqdm

from pyrewton.sql.sql_orm import get_db_connection
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

    if args.cazy.exists() is False:
        logger.error(
            f"Path to local CAZyme db ({args.cazy})\ndoes not exist\n"
            "Check the correct path was provided\n"
            "Terminating program"
        )

    # compile path to the output CAZome database
    if args.output_db is None:
        db_path = Path(f"cazome_db_{time_stamp}.db")
    else:
        db_path = args.output_db

    if args.output_dir is not None:
        file_io.make_output_directory(args.output_dir, args.force, args.nodelete)
        db_path = args.output_dir / db_path

    if db_path.exists() and args.force is False:
            logger.error(
                f"Db at {db_path} exists\nForce is false\nNot overwriting file\n"
                "To overwrite the file use -f or --force\nTerminating program"
            )
            sys.exit(1)

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

    # get paths to FASTA files of protein sequences parsed by dbCAN
    protein_fasta_files = get_paths.get_dir_paths(args.protein_dir)

    if len(protein_fasta_files) == 0:
        logger.error(
            f"No FASTA files containing protein sequences extracted from the genomes retrieved from\n"
            f"{args.protein_dir}\nCheck the correct dir was provided\n"
            "Terminating program"
        )
        sys.exit(1)

    # retrieve protein data from the FASTA file parsed by dbCAN
    #  {genus: {species: 'txid': str, 'genomes': {assembly_acc: {protein_acc: {'sequence': SeqIO.seq}}}}}
    all_data_dict = get_protein_data(protein_fasta_files)
    
    connection = get_db_connection(db_path)
    
    all_data_dict = get_dbcan_annotations(dbcan_output_dirs, all_data_dict)

    cache_path = Path(f"dbcan_predictions_{time_stamp}.json")

    if args.output_dir is not None:
        cache_path = args.output_dir / cache_path

    with open(cache_path, 'w') as fh:
        json.dump(all_data_dict, fh)

    prediction_df = build_prediction_df(all_data_dict)

    prediction_df_path = Path(f"dbcan_predictions_{time_stamp}.csv")

    if args.output_dir is not None:
        prediction_df_path = args.output_dir / prediction_df_path

    prediction_df.to_csv(prediction_df_path)

    # annotation_df = add_cazy_annotations(prediction_df)

    # annotation_df.to_csv('dbcan_cazy_annotations.csv')


def get_protein_data(protein_fasta_files):
    """Retrieve protein, genomic and taxonomic data from the FASTA files parsed by dbCAN.
    
    :param protein_fasta_files: list of Paths, one path per FASTA file
    
    Return dict {genus: {species: 'txid': str, 'genomes': {assembly_acc: {protein_acc: {'sequence': SeqIO.seq}}}}}
    """
    all_data_dict = {}

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

            try:
                all_data_dict[genus]

                try:
                    all_data_dict[genus][species]

                    try:
                        all_data_dict[genus][species]['genomes']

                        try:
                            all_data_dict[genus][species]['genomes'][assembly_acc]
                            
                            try:
                                all_data_dict[genus][species]['genomes'][assembly_acc][protein_acc]
                                
                                if all_data_dict[genus][species]['genomes'][protein_acc]['sequence'] != record.seq:
                                    logger.error(
                                        f"Multiple unique sequences found for {protein_acc}\n"
                                        "Additional digit added to protein accession to differentiate the proteins"
                                    )

                                    success = False
                                    while success is False:
                                        i = 0
                                        protein_acc += f"{str(i)}"
                                        try:
                                            all_data_dict[genus][species]['genomes'][protein_acc]
                                            i += 1
                                        except KeyError:
                                            all_data_dict[genus][species]['genomes'][protein_acc] = {
                                                'sequence': record.seq,
                                                'hmmer': set(),
                                                'hotpep': set(),
                                                'diamond': set(),
                                                'dbcan': set(),
                                                '#ofTools': set(),
                                            }
                                            success = True

                            except KeyError:  # protein accession raised KeyError
                                all_data_dict[genus][species]['genomes'][assembly_acc][protein_acc] = {
                                    'sequence': record.seq,
                                    'hmmer': set(),
                                    'hotpep': set(),
                                    'diamond': set(),
                                    'dbcan': set(),
                                    '#ofTools': set(),
                                }

                        except KeyError:  # assembly accession raised KeyError
                            all_data_dict[genus][species]['genomes'][assembly_acc] = {
                                protein_acc: {
                                    'sequence': record.seq,
                                    'hmmer': set(),
                                    'hotpep': set(),
                                    'diamond': set(),
                                    'dbcan': set(),
                                    '#ofTools': set(),
                                }
                            }

                    except KeyError:  # 'genomes' raised KeyError
                        all_data_dict[genus][species]['genomes'] = {assembly_acc:
                            {
                                protein_acc: {
                                    'sequence': record.seq,
                                    'hmmer': set(),
                                    'hotpep': set(),
                                    'diamond': set(),
                                    'dbcan': set(),
                                    '#ofTools': set(),
                                }
                            }
                        }

                except KeyError:  # species raised KeyError
                    all_data_dict[genus] = {species:
                        {
                            'txid': tax_id,
                            'genomes': {
                                assembly_acc: {
                                    protein_acc: {
                                        'sequence': record.seq,
                                        'hmmer': set(),
                                        'hotpep': set(),
                                        'diamond': set(),
                                        'dbcan': set(),
                                        '#ofTools': set(),
                                    }
                                }
                            },
                        }
                    }

            except KeyError:   # genus raised KeyError
                all_data_dict[genus] = {species:
                    {
                        'txid': tax_id,
                        'genomes': {
                            assembly_acc: {
                                protein_acc: {
                                    'sequence': record.seq,
                                    'hmmer': set(),
                                    'hotpep': set(),
                                    'diamond': set(),
                                    'dbcan': set(),
                                    '#ofTools': set(),
                                }
                            }
                        },
                    }
                }

    return all_data_dict


def get_dbcan_annotations(dbcan_output_dirs, all_data_dict):
    """Retrieve CAZy family annotations from dbCAN ouput.
    
    :param dbcan_output_dirs: list of Paths() to output dirs containing dbCAN output
    :param all_data_dict: 
    {genus: {species: 'txid': str, 'genomes': {assembly_acc: 
        {protein_acc: {'sequence': SeqIO.seq, 'diamond':set(), 'hmmer':set(), 'hotpep':set(), '#ofTools': int, 'dbcan':set()}}}}
    }
    
    Return dict:
    {genus: {species: 'txid': str, 'genomes': {assembly_acc: 
        {protein_acc: {'sequence': SeqIO.seq, 'diamond':set(), 'hmmer':set(), 'hotpep':set(), '#ofTools': int, 'dbcan':set()}}}}
    }
    """
    for dbcan_dir in tqdm(dbcan_output_dirs, "Parsing dbCAN output"):
        # dir name format: GCA_123456789_1
        genomic_accession = f"{(dbcan_dir.name).split('_')[0]}_{(dbcan_dir.name).split('_')[1]}.{(dbcan_dir.name).split('_')[2]}"

        # identify the genus and species
        genus, species = get_genome_tax(all_data_dict, genomic_accession)

        if genus is None and species is None:
            logger.warning(
                f"No tax data retrieved from FASTA files parsed by dbCAN retrieved from {genomic_accession}\n"
                "Not retrieving data for this genomic assembly"
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

            # create a CazymeProteinPrediction instance each for HMMER, Hotpep and DIAMOND
            hmmer_predictions = list(get_hmmer_prediction(line[1], protein_accession))
            hotpep_predictions = list(get_hotpep_prediction(line[2], protein_accession))
            diamond_predictions = list(get_diamond_prediction(line[3], protein_accession))

            # get the dbCAN consensus result
            no_of_tools = line[-1]
            
            hmmer_hotpep = list(set(hmmer_predictions) & set(hotpep_predictions))
            hmmer_diamond = list(set(hmmer_predictions) & set(diamond_predictions))
            hotpep_diamond = list(set(hotpep_predictions) & set(diamond_predictions))

            dbcan_predictions = list(set(hmmer_hotpep + hmmer_diamond + hotpep_diamond))
                        
            # add data from dbcan_dict to overview_dict which contains data for all proteins in the FASTA file
            try:
                all_data_dict[genus][species]['genomes'][genomic_accession][protein_accession]

                all_data_dict[genus][species]['genomes'][genomic_accession][
                    protein_accession]['hmmer'] += list(hmmer_predictions)
                all_data_dict[genus][species]['genomes'][genomic_accession][
                    protein_accession]['hotpep'] += list(hotpep_predictions)
                all_data_dict[genus][species]['genomes'][genomic_accession][
                    protein_accession]['diamond'] += list(diamond_predictions)
                all_data_dict[genus][species]['genomes'][genomic_accession][
                    protein_accession]['dbcan'] += list(dbcan_predictions)
                all_data_dict[genus][species]['genomes'][genomic_accession][
                    protein_accession]['#ofTools'] = no_of_tools
            except KeyError:
                logger.warning(
                    f"Retrieving dbCAN output for protein {protein_accession} from {genomic_accession}\n"
                    "but protein not retrieved from FASTA file parsed by dbCAN\n"
                    "Skipping protein"
                )
                continue

    return all_data_dict


def get_genome_tax(all_data_dict, genomic_accession):
    """Retrieve the genus and species of the source organism of a genome from the all_data_dict
    
    :param all_data_dict: 
    {genus: {species: 'txid': str, 'genomes': {assembly_acc: 
        {protein_acc: {'sequence': SeqIO.seq, 'diamond':set(), 'hmmer':set(), 'hotpep':set(), '#ofTools': int, 'dbcan':set()}}}}
    }
    :param genomic_accession: str

    Return genus (str) and species (str)
    """
    for genus in all_data_dict:
        genus_species = all_data_dict[genus]
        for species in genus_species:
            if genomic_accession in list(all_data_dict[genus][species]['genomes'].keys()):
                return genus, species

    return None, None


def get_hmmer_prediction(hmmer_data, protein_accession):
    """Retrieve HMMER prediciton data from the dbCAN overview.txt file.

    :param hmmer_data: str, output data from HMMER written in the overview.txt file.
    :param protein_accession: str

    Return set of CAZy family predictions
    """
    logger = logging.getLogger(__name__)

    if hmmer_data == "-":
        return set()
    
    cazy_fams = set()

    predicted_domains = hmmer_data.split("+")
    for domain in predicted_domains:
        # separate the predicted CAZy family from the domain range
        domain_name = domain.split("(")[0]  # CAZy (sub)family

        if domain_name.find("_") != -1:
            try: 
                re.match(r"\D{2,3}\d+?_\D", domain_name).group()  # check unusal CAZy family formating
                cazy_fams.add(domain_name.split("_")[0])

            except AttributeError:  # raised if not an usual CAZy family format
                try:
                    re.match(r"\D{2,3}\d+?_\d+", domain_name).group()  # check if a subfamily
                    cazy_fams.add(domain_name.split("_")[0])
                except AttributeError:
                    logger.warning(
                        f"Unknown data type of {domain_name} for protein {protein_accession}, for HMMER.\n"
                        "Not adding as domain to CAZy family annotations for the genome"
                    )
        
        else:
            try:
                re.match(r"\D{2,3}\d+", domain_name).group()
                cazy_fams.add(domain_name)

            except AttributeError:  # raised if doesn't match expected CAZy family
                logger.warning(
                    f"Unknown data type of {domain_name} for protein {protein_accession}, for HMMER.\n"
                    "Not adding as domain to CAZy family annotations for the genome"
                )
    
    return cazy_fams


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
                    logger.warning(
                        f"Unexpected data format of '{domain}' for protein "
                        f"{protein_accession}, for DIAMOND.\n"
                        "Not adding as domain to CAZy family annotations for the genome"
                    )

    return cazy_fams


def build_prediction_df(all_data_dict):
    """Build Pandas df of data in the dict.
    
    :param all_data_dict: 
    {genus: {species: 'txid': str, 'genomes': {assembly_acc: 
        {protein_acc: {'sequence': SeqIO.seq, 'diamond':set(), 'hmmer':set(), 'hotpep':set(), '#ofTools': int, 'dbcan':set()}}}}
    }
    
    Return Pandas df
    """
    column_names = [
        'Genus',
        'Species',
        'NCBI_txid',
        'Genomic_accession',
        'Protein_accession',
        'HMMER',
        'Hotpep',
        'DIAMOND',
        '#dbCAN_tools',
        'dbCAN',
    ]

    prediction_df = pd.DataFrame(columns=column_names)

    for genus in tqdm(all_data_dict, desc='Building prediction df per genus'):
        genus_species = all_data_dict[genus]

        for species in tqdm(genus_species, desc='Building prediction df per species'):
            tax_id = all_data_dict[genus][species]['txid']
            assemblies = all_data_dict[genus][species]['genomes']

            for genomic_accession in tqdm(assemblies, desc='Building prediction df per genome'):
                protein_data = all_data_dict[genus][species]['genomes'][genomic_accession]

                for protein_accession in tqdm(protein_data, desc='Parsing protein data'):
                    # retrieve the prediction annotations
                    hmmer = ' '.join(all_data_dict[genus][species]['genomes'][genomic_accession]['hmmer'])
                    hotpep = ' '.join(all_data_dict[genus][species]['genomes'][genomic_accession]['hotpep'])
                    diamond = ' '.join(all_data_dict[genus][species]['genomes'][genomic_accession]['diamond'])
                    dbcan = ' '.join(all_data_dict[genus][species]['genomes'][genomic_accession]['dbcan'])

                    new_row = [[
                        genus,
                        species,
                        tax_id,
                        genomic_accession,
                        protein_accession,
                        protein_accession,
                        hmmer,
                        hotpep,
                        diamond,
                        all_data_dict[genus][species]['genomes'][genomic_accession]['no#tools'],
                        dbcan,
                    ]]
                    new_row = pd.DataFrame(new_row, columns=column_names)

                    prediction_df = prediction_df.append(new_row, ignore_index=True)

    return prediction_df


def add_cazy_annotations(prediction_df):
    row_index = 0

    cazy_column = []

    for row_index in tqdm(range(len(prediction_df['HMMER'])), desc="Adding annotations from CAZy"):
        row = prediction_df.iloc[row_index]
        protein_accession = row['Protein_accession']
        
        with Session(bind=connection) as session:
            db_query = session.query(Genbank, CazyFamily).\
                join(CazyFamily, Genbank.families).\
                filter(Genbank.genbank_accession==protein_accession).\
                all()
            
        if len(db_query) == 0:
            # protein not in CAZy
            cazy_column.append([''])
        
        else:
            column_data = []
            for obj in db_query:
                column_data.append(obj[1].family)
            column_data.sort()
            cazy_column.append(' '.join(column_data))

    prediction_df['CAZy'] = cazy_column

    return prediction_df


if __name__ == "__main__":
    main()
