#!/usr/bin/env python
# -*- coding: utf-8 -*-
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
"""Retrieves all cazymes for each species in an input dataframe.

:func:

"""

import io
import re
import sys

from typing import List, Optional

import pandas as pd

from bioservices import UniProt
from tqdm import tqdm
from urllib.error import HTTPError

from pyrewton.directory_handling import input_dir_get_cazyme_annotations, output_dir_handling_main
from pyrewton.loggers.logger_pyrewton_main import build_logger
from # create parser for this script and import 

def main(argv: Optional[List[str]] = None, logger: Optional[logging.logger] = None):
    """docstring summary.

    Extra.

    Return.
    """
    # Programme preparation:
    # Parser arguments
    # Check if namespace isn't passed, if not parser command-line
    if argv is None:
        # Parse command-line
        parser =  # build parser
        args = parser.parse_args()
    else:
        args = build_parser(argv).parse_args()
    
    # Initate logger
    # Note: log file only created if specificied at cmdline
    if logger is None:
        logger = build_logger("get_cazyme_annotations", args)
    logger.info("Run initated")

    # If specified output directory, create output directory
    if args.output is not sys.stdout:
        output_dir_handling_main.make_output_directory(args.output, logger, args.force, args.nodelete)
    
    # Open input dataframe
    logger.info("Opening input dataframe")
    input_df = input_dir_get_cazyme_annotations.get_input_df(args.df_input, logger)

    # Build protein dataframe
    protein_df = build_dataframe(input_df, args, logger)

def build_dataframe(input_df, args, logger):
    """Build dataframe containing all cazymes for each species.

    :param input_df: input dataframe containing genus and species of host organisms
    :param args: parser object
    :param logger: logger object

    Return dataframe.
    """
    # Create empty dataframe to add data to
    cazyme_df = pd.DataFrame(
        columns=[
            "Genus",
            "Species",
            "NCBI Taxonomy ID",
            "UniProt entry ID",
            "UniProt entry name",
            "UniProt assigned protein names",
            "EC Numbers",
            "Length (Aa)",
            "Mass (Da)",
            "Domains",
            "Domain count",
            "UniProt linked protein families",
            "GO IDs",
            "GO molecular function",
            "G0 biological process",
            "Sequence",
        ]
    )

    # parse input_df. retrieving call cazymes for each species from UniProtKB
    df_index = 0
    for df_index in tqdm(range(len(input_df["Genus"])), desc="Retrieving Uniprot data"):
        cazyme_df = cazyme_df.append(get_uniprotkb_data(input_df.iloc[df_index], logger), ignore_index=True)
        df_index += 1
    
    # for development purposes and will be removed before release
    print("====\nUniProt protein data:\n", cazyme_df)


def get_uniprotkb_data(df_row, logger):
    """Retrieve all cazyme entries from UniProtKB.

    :param df_row: pandas sereies
    :param logger: logger object

    Return dataframe.
    """
    # Establish data to be retrieved from UniProt
    columnlist = (
        "id,entry name, protein names,length,mass,domains,domain,"
        "families,"
        "go-id,go(molecular function),go(biological process),"
        "sequence"
    )

    try:
        # open connection to UniProt(), search and convert results to panndas df
        search_result_df = pd.read_table(
            io.StringIO(
                UniProt().search(
                    f'organism:"{df_row[0]} {df_row[1]}"',
                    columns=columnlist
                )
            )
        )

    except HTTPError:
        logger.warning(
            (
                f"Network error occured when searching UniProt for locus tag:{df_row[2]}.\n"
                "Returning null value 'NA' for all UniProt data"
            )
        )
        data = {
            "Genus": df_row[0],
            "Species": df_row[1],
            "NCBI Taxonomy ID": df_row[2],
            "UniProt entry ID": "NA",
            "UniProt entry name": "NA",
            "UniProt assigned protein names": "NA",
            "EC Numbers": "NA",
            "Length (Aa)": "NA",
            "Mass (Da)": "NA",
            "Domains": "NA",
            "Domain count": "NA",
            "UniProt linked protein families": "NA",
            "GO IDs": "NA",
            "GO molecular function": "NA",
            "G0 biological process": "NA",
            "Sequence": "NA",
        }
        return pd.DataFrame(
            data,
            columns=[
                "Genus",
                "Species",
                "NCBI Taxonomy ID",
                "UniProt entry ID",
                "UniProt entry name",
                "UniProt assigned protein names",
                "Length (Aa)",
                "Mass (Da)",
                "Domains",
                "Domain count",
                "UniProt linked protein families",
                "GO IDs",
                "GO molecular function",
                "G0 biological process",
                "Sequence",
            ],
        )

    # Once fixed in get_cazyme_annotations add error exception
    # for returning empty dataframe

    # Rename columns to match dataframe design and indicate UniProt source
    search_result_df.rename(
        columns={
            "Entry": "UniProtKB Entry ID",
            "Entry name": "UniProtKB Entry Names",
            "Protein names": "UniProtKB Protein Names",
            "Length": "Length (Aa)",
            "Mass": "Mass (Da)",
            "Protein families": "UniProtKB Linked Protein Families",
        }
    )
    # Retrieve EC number from 'UniProtKB Protein Names'
    # or return 'NA' if not included
    ec_numbers = []  # tuple containing all EC numbers, each list is for a unique protein
    # create tuples so correct left to add genus and species to dataframe
    genus = []  # tuple to store genus in
    species = [] # tuple to store species

    # use for loop or .apply?
    ec_numbers = search_result_df.apply(get_ec_number, args=logger, axis=1)

    row_index = 0
    for row_index in range(len(search_result_df["UniProtKB Entry ID"])):
        ec_numbers.append(get_ec_number(search_result_df.iloc[row_index], logger))
        genus.append([df_row[0]])
        species.appned([df_row[1]])
        row_index += 1
    
    # Add EC numbers column to dataframe
    search_result_df.insert(3, "EC Number", ec_number)

    # Add Genus and Species column to dataframe
    search_result_df.insert(0, "Genus", genus)
    search_result_df.insert(1, "Species", species)

    # Remove those with no indication they are cazymes
    search_result_df = remove_non_cazyme(search_result_df, logger)

    # Add genus and species columns

    return search_result_df


def get_ec_number(df_row, logger):
    """Retrieve EC Number from UniProt search results.

    :param df_row: Pandas series
    :param logger: logger object

    Return list of all EC numbers.
    """
    # Search protein name cell for EC numbers
    ec_search = re.findall(
        r"\(EC [\d-]\d*\.[\d-]\d*\.[\d-]\d*\.[\d-]\d*\)", df_row[2]
    )
    ec_numbers = []
    if ec_search is None:
        ec_numbers.append('NA')
    else:
        # compile EC numbers together if multiple were given
        for ec in ec_search:
            ec_numbers.append(ec)

    return ec_numbers


def remove_non_cazyme(df, logger):
    """Remove proteins with no indication of being CAZymes.

    :param df: Pandas dataframe
    :param logger: logger object

    Return dataframe.
    """
    # Create empty dataframe to add potential cazyme entries too
    filtered_df = pd.DataFrame(
        columns=[
            "Genus",
            "Species",
            "NCBI Taxonomy ID",
            "UniProt entry ID",
            "UniProt entry name",
            "UniProt assigned protein names",
            "EC Numbers",
            "Length (Aa)",
            "Mass (Da)",
            "Domains",
            "Domain count",
            "UniProt linked protein families",
            "GO IDs",
            "GO molecular function",
            "G0 biological process",
            "Sequence",
        ]
    )

    # Retrieve all rows with CAZyDB link
    df.query('')

    # Retrieve all rows with EC number or GO function indicating cazyme

    return filtered_df

