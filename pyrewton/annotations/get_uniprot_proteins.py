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
"""Retrieves all proteins for each species in an input dataframe.

:func main: coordianate function of script
:func build_uniprot_df: coordinates retrieval of UniProt entries per species
:func get_uniprotkb_data: retreives entries from UniProt
"""

import io
import logging
import re
import sys

from typing import List, Optional

import pandas as pd

from bioservices import UniProt
from pandas.errors import EmptyDataError
from tqdm import tqdm
from urllib.error import HTTPError

from pyrewton import file_io
from pyrewton.loggers import build_logger
from pyrewton.parsers.parser_get_uniprot_proteins import build_parser


def main(argv: Optional[List[str]] = None, logger: Optional[logging.logger] = None):
    """Coordinate retrieval of entries from UniProtKB database.

    Extra.

    Store entries in pandase dataframe; write out dataframe to csv file.
    """
    # Programme preparation:
    # Parser arguments
    # Check if namespace isn't passed, if not parser command-line
    if argv is None:
        # Parse command-line
        parser = build_parser()
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
        file_io.make_output_directory(args, logger)

    # Open input dataframe
    # Open input dataframe
    logger.info("Opening input dataframe %s", args.df_input)
    input_df = pd.read_csv(args.df_input, header=0, index_col=0)

    # Build protein dataframe
    uniprot_protein_df = build_uniprot_df(input_df, args, logger)

    # Write out Uniprot dataframe to csv file
    file_io.write_out_dataframe(
        uniprot_protein_df, logger, args.ouput, args.force, args.nodelete,
    )


def build_uniprot_df(input_df, args, logger):
    """Retrieve all proteins in UniProt for each species in an input dataframe.

    :param input_df: input dataframe containing genus and species of host organisms
    :param args: parser object
    :param logger: logger object

    Return dataframe.
    """
    # Create empty dataframe to add data to
    uniprot_df = pd.DataFrame(
        columns=[
            "Genus",
            "Species",
            "NCBI Taxonomy ID",
            "UniProt entry ID",
            "UniProt entry name",
            "UniProt assigned protein names",
            "EC Number",
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
        uniprot_df = uniprot_df.append(
            get_uniprotkb_data(input_df.iloc[df_index], logger), ignore_index=True
        )
        df_index += 1

    # for development purposes and will be removed before release
    print("====\nUniProt protein data:\n", uniprot_df)

    return uniprot_df


def get_uniprotkb_data(df_row, logger):
    """Retrieve all cazyme entries from UniProtKB for species in dataframe series/row.

    :param df_row: pandas series from input df, where:
        df_row[0] = Genus
        df_row[1] = Species
        df_row[2] = NCBI Taxonomy ID
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

    # This dictionary will be used to populate "blank"/"empty" databases when
    # an error is thrown. Iterables are used as values to avoid problems with
    # "ValueError: If using all scalar values, you must pass an index"
    blank_data = {
        "UniProtKB Entry ID": ["NA"],
        "UniProtKB Entry Name": ["NA"],
        "UniProtKB Protein Names": ["NA"],
        "EC number": ["NA"],
        "Length (Aa)": ["NA"],
        "Mass (Da)": ["NA"],
        "Domains": ["NA"],
        "Domain count": ["NA"],
        "UniProtKB Linked Protein Families": ["NA"],
        "Gene ontology IDs": ["NA"],
        "Gene ontology (molecular function)": ["NA"],
        "Gene ontology (biological process)": ["NA"],
    }

    try:
        # open connection to UniProt(), search and convert result into pandas df
        logger.info(df_row)
        query = f'{df_row[5]} AND organism:"{df_row[0]} {df_row[1]}"'
        search_result = UniProt().search(
            query, columns=columnlist,
        )  # returns empty string for no result
        logger.info(search_result)
        search_result_df = pd.read_table(io.StringIO(search_result))

    except HTTPError:
        logger.warning(
            (
                f"Network error occured when searching UniProt for locus tag:{df_row[5]}.\n"
                "Returning null value 'NA' for all UniProt data"
            )
        )
        return pd.DataFrame(blank_data)

    except EmptyDataError:
        # No UniProt entries found for locus tag, return null data for
        logger.warning(
            (
                f"No data returned from UniProt for locus tag:{df_row[5]}.\n"
                "Returning null value 'NA' for all UniProt data"
            )
        )
        return pd.DataFrame(blank_data)

    # rename columns to match to indicate UniProtKB source of data
    search_result_df.rename(
        columns={
            "Entry": "UniProtKB Entry ID",
            "Entry name": "UniProtKB Entry Name",
            "Protein names": "UniProtKB Protein Names",
            "Length": "Length (Aa)",
            "Mass": "Mass (Da)",
            "Protein families": "UniProtKB Linked Protein Families",
        }
    )

    # Retrieve EC number from 'UniProtKB Protein Names'
    # or return 'NA' if not included
    EC_search = re.findall(
        r"\(EC [\d-]\d*\.[\d-]\d*\.[\d-]\d*\.[\d-]\d*\)", search_result_df[2]
    )
    if EC_search is None:
        EC_number = "NA"
    else:
        # compiall EC together incase multiple are given
        # and remove EC numbers from protein name
        EC_number = ""
        for EC in EC_search:
            search_result_df[2] = search_result_df[2].replace(EC, "")
            EC = EC.replace("(", "")
            EC = EC.replace(")", "")
            EC_number += EC
            # Do I want to remove the EC number becuase if multiple are given,
            # maybe helpful to see in the protein name section which EC number
            # corresponds to which function
    # Add EC number to dataframe0
    search_result_df.insert(3, "EC number", EC_number)

    # Add genus, species and NCBI taxonomy ID column, so the host organism is identifable
    # for each protein
    genus_column_data = [df_row[0]] * range(len("UniProtKB Entry ID"))
    species_column_data = [df_row[1]] * range(len("UniProtKB Entry ID"))
    tax_id_column_data = [df_row[2]] * range(len("UniProtKB Entry ID"))

    search_result_df.insert(0, "Host Genus", genus_column_data)
    search_result_df.insert(1, "Host Species", species_column_data)
    search_result_df.insert(2, "Host NCBI Taxonomy ID", tax_id_column_data)

    return search_result_df
