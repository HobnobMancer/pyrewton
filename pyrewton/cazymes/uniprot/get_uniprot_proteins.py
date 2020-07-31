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
import logging
import re
import sys
import yaml

from typing import List, Optional
from datetime import datetime

import pandas as pd

from bioservices import UniProt
from pandas.errors import EmptyDataError
from tqdm import tqdm
from urllib.error import HTTPError

from pyrewton.file_io import write_out_pre_named_dataframe
from pyrewton.loggers import build_logger
from pyrewton.parsers.parser_get_uniprot_proteins import build_parser


def main(argv: Optional[List[str]] = None, logger: Optional[logging.Logger] = None):
    """Set up parser, loggers, IO files and directories, then invoke scripts main function.

    Set up loggers, parsers and directories for retrieval of cazymes from UniProtKB.
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
    logger.info("Run initated")

    # Retrieve data from configuration file
    tax_ids, query_list = get_config_data(logger, args)

    # Iterate over species to search using each query from config file
    for tax_id in tax_ids:
        for query in query_list:
            build_uniprot_df(tax_id, query, logger, args)

    logger.info("Program finished")


def get_config_data(logger, args):
    """Retrieve data from configration file.

    :param logger: logger objects
    :param args: parser arguments.

    Returns list of taxonomy IDs and list of queries.
    """
    logger.info("Retrieving queries from config file")
    with open(args.input) as ifh:
        config_dict = yaml.safe_load(ifh)

    # Retrieve Taxonomy IDs from configuration data
    tax_ids = config_dict["tax_ids"]

    # Retrieve other queries from configuration data
    query_list = []
    query_list.append(config_dict["queries"])

    return tax_ids, query_list


def build_uniprot_df(tax_id, query, logger, args):
    """Build dataframe to store data retrieved from UniProtKB.

    :param tax_id: NCBI taxonomy ID of species
    :param query: term to search for in specified
    :param logger: logger object
    :param args: parser arguments

    Returns nothing.
    """
    # Call UniProtKB and return results as dataframe
    uniprot_df = call_uniprotkb(tax_id, query, logger, args)

    # Rename columns and create separate column to store EC numbers, and write
    # out sequences to FASTA files if enabled
    uniprot_df = format_search_results(uniprot_df, tax_id, logger, args)

    # write out resulting dataframe for UniProtKB query
    invalid_file_name_characters = re.compile(r'[,;./ "*#<>?|\\:]')
    query_title = re.sub(invalid_file_characters, "_", query)
    dataframe_name = f"UniProt_{tax_id}_{query_title}_{now.strftime("%Y-%m-%d_%H-%M-%S")}.csv"
    write_out_pre_named_dataframe(uniprot_df, dataframe_name, logger, args.output, args.force)

    return


def call_uniprotkb(tax_id, query, logger, args):
    """Calls to UniProt.

    If no data is retieved a default 'blank' dataframe is returned.

    :param tax_id: NCBI taxonomy ID of species
    :param query: query for UniProt
    :param logger: logger object
    :param args: parser arguments

    Returns dataframe of search results.
    """
    # Establish data to be retrieved from UniProt
    columnlist = (
        "organism,id,entry name, protein names,length,mass,domains,domain,"
        "families,"
        "go-id,go(molecular function),go(biological process),"
        "sequence"
    )

    # This dictionary will be used to populate "blank"/"empty" databases when
    # an error is thrown. Iterables are used as values to avoid problems with
    # "ValueError: If using all scalar values, you must pass an index"
    blank_data = {
        "Organism": ["NA"],
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
        "Sequence": ["NA"],
    }

    try:
        # open connection to UniProt(), search and convert result into pandas df
        search_result = UniProt().search(
            query, columns=columnlist,
        )  # returns empty string for no result

        return pd.read_table(io.StringIO(search_result))

    except HTTPError:
        logger.warning(
            (
                f"Network error occured when searching UniProt for locus tag:{tax_id}.\n"
                "Returning null value 'NA' for all UniProt data"
            )
        )
        return pd.DataFrame(blank_data)

    except EmptyDataError:
        # No UniProt entries found for locus tag, return null data for
        logger.warning(
            (
                f"No data returned from UniProt for locus tag:{tax_id}.\n"
                "Returning null value 'NA' for all UniProt data"
            )
        )
        return pd.DataFrame(blank_data)


def format_search_results(search_result_df, tax_id, logger, args):
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
    for index in tqdm(range(len(search_result_df["UniProtKB Entry ID"])), desc="Processing protein data")
        df_row = search_result_df.iloc[index]
        all_ec_numbers.append(get_ec_numbers(df_row, logger))
        # Write out protein sequence to FASTA file if sequence was retrieved from UniProtKB
        if (args.fasta) and (df_row["Sequence"] != "NA"):
            write_fasta(df_row, logger, args)
        index += 1

    # Add EC numbers to dataframe0
    search_result_df.insert(3, "EC number", all_ec_numbers)

    # Add genus, species and NCBI taxonomy ID column, so the host organism is identifable
    # for each protein
    tax_id_column_data = tax_id * len(all_ec_numbers)

    search_result_df.insert(0, "NCBI Taxonomy ID", tax_id_column_data)

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
    file_content = f">{df_row["NCBI Taxonomy ID"]} {df_row[""]}"

    return


if __name__ == "__main__":
    main()
