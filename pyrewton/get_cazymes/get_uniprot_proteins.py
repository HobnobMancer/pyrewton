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
"""Retrieves proteins from UniProt for given species,
and create fasta file of protein sequence.

:param ...

:func main: set-up script
"""

import io
import json
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


def main(argv: Optional[List[str]] = None, logger: Optional[logging.Logger] = None):
    """Coordinate retrieval of entries from UniProtKB database.

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

    # If specified output directory, create output directory to write FASTA files too
    if args.output is not sys.stdout:
        file_io.make_output_directory(args, logger)

    # Initate scripts main function


def get_uniprot_proteins(args, logger):
    """Coordinate overall function of script to retrieve protein data.

    :param args: parser arguments
    :param logger: logger object
    """
    logger.info("Run initated")

    # Open input dataframe, containing species dataframe
    logger.info("Opening input dataframe %s", args.df_input)
    input_df = pd.read_csv(args.df_input, header=0, index_col=0)

    # Open file containing query terms, saving dictionary with query field
    # as the Key and query term as the value
    with (args.query).open("r") as fh:
        query_dict = json.load(fh)

    # Iterate over species in input dataframe, to search for potential CAZymes
    df_index = 0
    for df_index in tqdm(range(len(input_df["Genus"])), desc="Retrieving Uniprot data"):
        pd_series = input_df.iloc[df_index]
        # Iterate through query fields and query terms for each field
        for key in query_dict:
            query_list = query_dict[key]
            for term in query_list:
                # Call function which coordinates call to UniProt
                build_uniprot_df(
                    f"{pd_series[0]}",
                    f"{pd_series[1]}",
                    pd_series[2],
                    key,
                    term,
                    logger,
                    args,
                )
        df_index += 1

    logger.info("Program finished")


def build_uniprot_df(genus, species, tax_id, query_field, query_term, logger, args):
    """Coordinates call to UniProt and formatting of results.

    :param genus: genus of host species
    :param species: species name
    :param tax_id: NCBI taxonomy ID of species
    :param query_field: field in UniProt to search using query term
    :param query_term: term to search for in specified
    :param logger: logger object
    :param args: parser arguments

    Returns nothing.
    """
    # Call UniProtKB and return results as dataframe
    uniprot_df = call_uniprotkb(
        genus, species, tax_id, query_field, query_term, logger, args
    )

    # Rename columns and create separate column to store EC numbers
    uniprot_df = format_search_results(uniprot_df, genus, species, tax_id, logger, args)

    return


def call_uniprotkb(genus, species, tax_id, query_field, query_term, logger, args):
    """Call to UniProt.

    If no data is retieved a default 'blank' dataframe is returned.

    :param genus: genus of host species
    :param species: species name
    :param tax_id: NCBI taxonomy ID of species
    :param query_field: field in UniProt to search using query term
    :param query_term: term to search for in specified
    :param logger: logger object
    :param args: parser arguments

    Returns dataframe of search results.
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
        "Sequences": ["NA"],
    }

    try:
        # construct query
        query = f'organism:"{genus} {species}" {query_field}:"{query_term}"'
        # open connection to UniProt(), search and convert result into pandas df
        search_result = UniProt().search(
            query, columns=columnlist,
        )  # returns empty string for no result

        return pd.read_table(io.StringIO(search_result))

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


def format_search_results(search_result_df, genus, species, tax_id, logger, args):
    """Rename columns, add EC number, genus, species and tax ID columns.

    :param search_result_df: pandas dataframe of UniProt search results
    :param genus: genus of species
    :param species: species name
    :param tax_id: NCBI taxonomy ID of species
    :param logger: logger object
    :param args: parser arguments

    Return pandas dataframe.
    """
    logger.info("Renaming column headers")
    search_result_df = search_result_df.rename(
        columns={
            "Entry": "UniProtKB Entry ID",
            "Entry name": "UniProtKB Entry Name",
            "Protein names": "UniProtKB Protein Names",
            "Length": "Length (Aa)",
            "Mass": "Mass (Da)",
            "Protein families": "UniProtKB Linked Protein Families",
        }
    )

    logger.info("Retrieving EC numbers")
    index = 0
    all_ec_numbers = []  # list, each item is a str of all EC numbers in a unique row
    for index in range(len(search_result_df["UniProtKB Entry ID"])):
        df_row = search_result_df.iloc[index]
        all_ec_numbers.append(get_ec_numbers(df_row, logger))
        # Write out protein sequence to FASTA file if sequence was retrieved from UniProtKB
        if df_row["Sequences"] is not "NA":
            write_fasta(df_row, logger, args)
        index += 1

    # Add EC numbers to dataframe0
    search_result_df.insert(3, "EC number", all_ec_numbers)

    # Add genus, species and NCBI taxonomy ID column, so the host organism is identifable
    # for each protein
    genus_column_data = genus * len(all_ec_numbers)
    species_column_data = species * len(all_ec_numbers)
    tax_id_column_data = tax_id * len(all_ec_numbers)

    search_result_df.insert(0, "Genus", genus_column_data)
    search_result_df.insert(1, "Species", species_column_data)
    search_result_df.insert(2, "NCBI Taxonomy ID", tax_id_column_data)

    return search_result_df


def get_ec_numbers(df_row, logger):
    """Retrieve EC numbers in dataframe row.

    :param df_row: pandas series
    :param logger: logger object

    Returns str of all EC numbers retrieved (human readible list).
    """
    ec_search = re.findall(r"\(EC [\d-]\d*\.[\d-]\d*\.[\d-]\d*\.[\d-]\d*\)", df_row[2])

    if ec_search is None:
        ec_numbers = "NA"

    else:
        # compiall EC numbers together incase multiple are given
        # and remove EC numbers from protein name
        ec_numbers = ""
        for ec in ec_search:
            ec = ec.replace("(", "")
            ec = ec.replace(")", "")
            ec_numbers += ec

    return ec_numbers


def write_fasta(df_row, logger, args):
    """Write out FASTA file.

    :param df_row: row from pandas df of UniProt search results
    :param logger: logger object
    :param args: parser arguments

    Returns nothing.
    """

    return


if __name__ == "__main__":
    main()

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

    # Call to UniProt
    call_uniprotkb

    return uniprot_df
