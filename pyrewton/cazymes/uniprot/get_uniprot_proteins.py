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
"""Retrieves proteins from UniProt and write fasta files.

Queries and Taxonomy IDs are passed to this script via a YAML
file.
Taxonomy IDs must be stored as a list under 'tax_ids'.
User defined queries must be written under 'queries' as a list,
and written in the UniProtKB query syntax (see
uniprot.ord/help/text-syntax).
Write all tax IDs and queries within quotation marks,
such as "database:(type:cazy)".

:param input: optional, path to configuration file
:param fasta: optional, enable writing out of fasta files
:param force: optional, force overwrite if output already exists
:param log: optional, enable writing out of log file
:param nodelete: optional, enable no deletion of content in output dir
:param output: optional, path to output dir
:param verbose: optional, change logger level to 'info'

:func main: set-up script, configure call to UniProtKB
:func get_config_data: retrieve data from config file
:func build_uniprot_df: build query, coordinate dataframe formating
:func call_uniprotkb: call to UniProtKB
:func format_search_result: rename columns, add EC number column
:func get_ec_numbers: retrieve EC numbers for UniProt dataframe
:func write_fasta: write out data to fasta file

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

from pyrewton.file_io import write_out_pre_named_dataframe, make_output_directory
from pyrewton.loggers import build_logger
from pyrewton.parsers.parser_get_uniprot_proteins import build_parser


def main(argv: Optional[List[str]] = None, logger: Optional[logging.Logger] = None):
    """Set up parser, loggers, and IO directories, then invoke scripts main function."""
    # Parser arguments
    # Check if namespace isn't passed, if not parser command-line
    if argv is None:
        parser = build_parser()
        args = parser.parse_args()
    else:
        args = build_parser(argv).parse_args()

    # Initate logger
    # Note: log file only created if specificied at cmdline
    if logger is None:
        logger = build_logger("get_cazyme_annotations", args)

    # If specified output directory, create output directory to write FASTA files too
    if args.output is not sys.stdout:
        make_output_directory(args, logger)

    # Initate scripts main function
    configuration(args, logger)


def configuration(args, logger):
    """Coordinate calling to UniProtKB."""
    # Retrieve data from configuration file
    tax_ids, query_list = get_config_data(logger, args)

    # Mediate iteration of tax IDs and/or queries, as retrieved from config file

    # Search by user defined query only
    if tax_ids is None:
        for query in tqdm(range(len(query_list)), desc="Querying UniProtKB"):
            build_uniprot_df(None, query, logger, args)

    # Search by taxonomy ID only
    elif query_list is None:
        for tax_id in tqdm(range(len(tax_ids)), desc="Querying UniProtKB"):
            build_uniprot_df(tax_id, None, logger, args)

    # Search every user defined query with ('AND') for each taxonomy ID
    else:
        id_index = 0
        for id_index in range(len(tax_ids)):
            query_index = 0
            for query in tqdm(range(len(query_list)), desc="Querying UniProtKB"):
                build_uniprot_df(
                    tax_ids[id_index], query_list[query_index], logger, args
                )
                query_index += 1
            id_index += 1

    logger.info("Program finished")


def get_config_data(logger, args):
    """Retrieve data from configration file.

    :param logger: logger objects
    :param args: parser arguments.

    Returns list of taxonomy IDs and list of queries.
    """
    logger.info("Retrieving queries from config file")
    with open(args.input) as ifh:
        config_dict = yaml.full_load(ifh)

    # Retrieve Taxonomy IDs from configuration data
    try:
        tax_ids = config_dict["tax_ids"]
    except (KeyError, TypeError) as e:
        logger.warning(
            (
                "No Taxonomy IDs retrieved from configuration file\n."
                "If not restricting search by organism continue.\n"
                "Else make sure Taxonomy IDs are under the heading 'tax_ids:'"
            )
        )
        tax_ids = None
    # If tax IDs key exists but no tax IDs are stored underneath return None
    if tax_ids is not None:
        if len(tax_ids) == 0:
            tax_ids = None

    # Retrieve user defined queries from configuration data
    query_list = []
    try:
        query_list.append(config_dict["queries"])
    except (KeyError, TypeError) as e:
        logger.warning(
            (
                "No queries retrieved from configuration file\n."
                "If no additional queriers wanted, continue.\n"
                "Else make sure additional queries are under the 'queries:'"
            )
        )
        query_list = None

    if query_list is not None:
        if len(query_list) == 0:
            query_list = None

    if (tax_ids is None) and (query_list is None):
        logger.error(
            (
                "Nothing retrieved from configuration file.\n"
                "Ensure correct file was passed to script, and file contains configuration data.\n"
                "Terminating program."
            ),
            sys.exit(1),
        )

    return tax_ids, query_list


def build_uniprot_df(tax_id, query, logger, args):
    """Build dataframe to store data retrieved from UniProtKB.

    :param tax_id: str, NCBI taxonomy ID of species
    :param query: str, term to search for in specified
    :param logger: logger object
    :param args: parser arguments

    Returns nothing.
    """
    # remove NCBI:txid prefix, otherwise UniProtKB will return no results
    if tax_id is not None:
        tax_id = tax_id.replace("NCBI:txid", "")

    # construct query and filestem for dataframes and FASTA files
    time_stamp = datetime.now().strftime("%Y-%m-%d--%H-%M-%S")
    if tax_id is None:
        uniprot_query = query
    elif query is None:
        uniprot_query = f'taxonomy:"{tax_id}"'
    else:
        uniprot_query = f'taxonomy:"{tax_id}" AND {query[0]}'

    filestem = f"uniprot_{uniprot_query}_{time_stamp}"

    # Call UniProtKB and return results as dataframe
    uniprot_df = call_uniprotkb(uniprot_query, logger)

    # remove characters that could make file names invalid
    invalid_file_name_characters = re.compile(r'[,;./ "*#<>?|\\:]')
    filestem = re.sub(invalid_file_name_characters, "_", filestem)

    # Rename columns and create separate column to store EC numbers, and write
    # out sequences to FASTA files if enabled
    uniprot_df = format_search_results(uniprot_df, filestem, logger, args)

    # write out resulting dataframe for UniProtKB query
    write_out_pre_named_dataframe(
        uniprot_df, f"{filestem}", logger, args.output, args.force
    )

    return


def call_uniprotkb(query, logger):
    """Calls to UniProt.

    If no data is retieved a default 'blank' dataframe is returned.

    :param query: str, query for UniProt
    :param logger: logger object

    Returns dataframe of search results.
    """
    # Establish data to be retrieved from UniProt
    columnlist = (
        "organism-id,organism,id,entry name, protein names,length,mass,domains,domain,"
        "families,"
        "go-id,go(molecular function),go(biological process),"
        "sequence"
    )

    # This dictionary will be used to populate "blank"/"empty" databases when
    # an error is thrown. Iterables are used as values to avoid problems with
    # "ValueError: If using all scalar values, you must pass an index"
    blank_data = {
        "NCBI Taxonomy ID": ["NA"],
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

    logger.info("querying uniprot, query: {query}")
    try:
        # open connection to UniProt(), search and convert result into pandas df
        search_result = UniProt().search(
            query, columns=columnlist,
        )  # returns empty string for no result

        return pd.read_table(io.StringIO(search_result))

    except HTTPError:
        logger.warning(
            (
                f"Network error occured during query: {query}\n"
                "Returning null value 'NA' for all UniProt data"
            )
        )
        return pd.DataFrame(blank_data)

    except EmptyDataError:
        # No UniProt entries found for locus tag, return null data for
        logger.warning(
            (
                f"No data returned from UniProt during query: {query}\n"
                "Returning null value 'NA' for all UniProt data"
            )
        )
        return pd.DataFrame(blank_data)


def format_search_results(search_result_df, filestem, logger, args):
    """Rename columns, add EC number, genus, species and tax ID columns.

    :param search_result_df: pandas dataframe of UniProt search results
    :param tax_id: str, NCBI taxonomy ID of species
    :param filestem: str, FASTA file name
    :param logger: logger object
    :param args: parser arguments

    Return pandas dataframe.
    """
    logger.info("Renaming column headers")
    search_result_df = search_result_df.rename(
        columns={
            "Organism ID": "NCBI Taxonomy ID",
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
    all_ec_numbers = []  # list, each item as a str of all EC numbers in a unique row
    for index in tqdm(
        range(len(search_result_df["UniProtKB Entry ID"])),
        desc="Processing protein data",
    ):
        df_row = search_result_df.iloc[index]
        all_ec_numbers.append(get_ec_numbers(df_row, logger))
        # Write out protein sequence to FASTA file if sequence was retrieved from UniProtKB
        if (args.fasta) and (df_row["Sequence"] != "NA"):
            write_fasta(df_row, filestem, logger, args)
        index += 1

    # Add EC numbers to dataframe0
    try:
        search_result_df.insert(3, "EC number", all_ec_numbers)
    except ValueError:
        logger.warning(
            ("Failed to insert EC numbers into dataframe,\n" "column already present")
        )

    return search_result_df


def get_ec_numbers(df_row, logger):
    """Retrieve EC numbers in dataframe row.

    :param df_row: pandas series
    :param logger: logger object

    Returns str of all EC numbers retrieved (human readible list).
    """
    ec_search = re.findall(r"\(EC [\d-]\d*\.[\d-]\d*\.[\d-]\d*\.[\d-]\d*\)", df_row[4])

    if ec_search is None:
        ec_numbers = "NA"

    else:
        # compiall EC numbers together incase multiple are given
        # and remove EC numbers from protein name
        ec_numbers = ""
        for ec in ec_search:
            ec = ec.replace("(", "")
            ec = ec.replace(")", "")
            ec_numbers += ec + ", "

    # remove terminal ", "
    return ec_numbers[:-2]


def write_fasta(df_row, filestem, logger, args):
    """Write out FASTA file.

    :param df_row: row from pandas df of UniProt search results
    :param filestem: str, FASTA file name
    :param logger: logger object
    :param args: parser arguments

    Returns nothing.
    """
    # FASTA sequences have 60 characters per line, add line breakers into protein sequence
    # to match FASTA format
    sequence = df_row["Sequence"]
    sequence = "\n".join([sequence[i : i + 60] for i in range(0, len(sequence), 60)])

    # Retrieve Taxonomy ID and ensure NCBI prefix is present
    uniprot_tax_id = str(df_row["NCBI Taxonomy ID"])
    if uniprot_tax_id.startswith("NCBI:txid") is False:
        uniprot_tax_id = "NCBI:txid" + uniprot_tax_id
    uniprot_tax_id.replace(" ", "")

    # Retrieve organism name and protein id
    organism = df_row["Organism"]
    protein_id = df_row["UniProtKB Entry ID"]

    file_content = f">{protein_id} {uniprot_tax_id} {organism} \n{sequence}\n"

    # Create output path
    if args.output is not sys.stdout:
        output_path = args.output / f"{filestem}.fasta"
    else:
        output_path = args.output

    # Write out data to Fasta file
    with open(output_path, "a") as fh:
        fh.write(file_content)

    return


if __name__ == "__main__":
    main()
