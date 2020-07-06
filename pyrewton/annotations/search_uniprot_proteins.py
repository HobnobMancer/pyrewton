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
"""Search dataframe containing UniProt entries for cazymes.

Cazymes identified by UniProt entry being linekd to a CAZyDB entry.
Potential cazymes identifed by EC number indicating cazyme functionality,
and by GO annotated function inferring cazyme functionality.

:func main: coordianate searching for cazymes in local UniProt dataframe.

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

    # Retrieve UniProt entries with link to CAZy database
    cazy_linked_df = get_cazy_proteins(input_df, logger)
    # write out df to csv file

    # Retrieve UniProt entries whose EC number infers cazyme functionality
    ec_inferred_df =
    # write out df to csv file

    # Retrieve UniProt entries whose GO annotated function infers cazyme funcitonality
    go_inferred_df =
    # write out df to csv file


    # Build protein dataframe
    uniprot_protein_df = build_uniprot_df(input_df, args, logger)

    # Write out Uniprot dataframe to csv file
    file_io.write_out_dataframe(
        uniprot_protein_df, logger, args.ouput, args.force, args.nodelete,
    )


def get_cazy_cazymes(input_df, logger):
    """Retrieve subset of entries with link to CAZy database.

    :func input_df: pandas dataframe
    :func logger: logger, object

    Return pandas dataframe.
    """
    # retrieve indexes of entries with link to CAZy database
    logger.info("Retrieving entries with link to CAZy database")
    cazy_link_indexes = input_df["UniProt linked protein families"].str.contains(
        r"cazy", flags=re.IGNORECASE, regex=True, na=False
    )
    return input_df[cazy_link_indexes]


def get_ec_cazymes(input_df, logger):
    """Retrieve subset of entries with indicated cazyme functionality
    from EC number(s).

    The EC numbers are searched for in number order.

    Return pandas dataframe.
    """
    ec_series = [] # store indexed pandas series results of EC search results
    ec_series.append(input_df["EC number"].str.contains(
        r"3.2.1.4 ", flags=re.IGNORECASE, regex=True, na=False
    ))
    ec_series.append(input_df["EC number"].str.contains(
        r"3.2.1.4 ", flags=re.IGNORECASE, regex=True, na=False
    ))
    ec_series.append(input_df["EC number"].str.contains(
        r"3.2.1.91", flags=re.IGNORECASE, regex=True, na=False
    ))
    ec_series.append(input_df["EC number"].str.contains(
        r"3.2.1.21", flags=re.IGNORECASE, regex=True, na=False
    ))
    ec_series.append(input_df["EC number"].str.contains(
        r"3.2.1.6", flags=re.IGNORECASE, regex=True, na=False
    ))
    ec_series.append(input_df["EC number"].str.contains(
        r"3.2.1.7", flags=re.IGNORECASE, regex=True, na=False
    ))
    ec_series.append(input_df["EC number"].str.contains(
        r"3.2.1.25", flags=re.IGNORECASE, regex=True, na=False
    ))
    ec_series.append(input_df["EC number"].str.contains(
        r"3.2.1.4", flags=re.IGNORECASE, regex=True, na=False
    ))
    ec_series.append(input_df["EC number"].str.contains(
        r"3.2.1.8", flags=re.IGNORECASE, regex=True, na=False
    ))
    ec_series.append(input_df["EC number"].str.contains(
        r"3.2.1.37", flags=re.IGNORECASE, regex=True, na=False
    ))
    ec_series.append(input_df["EC number"].str.contains(
        r"3.2.1.156", flags=re.IGNORECASE, regex=True, na=False
    ))
    ec_series.append(input_df["EC number"].str.contains(
        r"3.2.1.15", flags=re.IGNORECASE, regex=True, na=False
    ))
    ec_series.append(input_df["EC number"].str.contains(
        r"3.2.1.67", flags=re.IGNORECASE, regex=True, na=False
    ))
    
    # Create empty dataframe to append search results to
    ec_result_df = pd.DataFrame(
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
    
    # Retrieve rows from input dataframe whose EC numbers infer cazyme functionality
    # and combine into single dataframe, and remove duplicates
    for pandas_series in ec_series:
        ec_result_df = ec_result_df.append(input_df[pandas_series])
    ec_result_df = ec_result_df[ec_result_df.duplicated(subset="UniProt entry ID", keep='first')
    return ec_result_df


def get_go_cazymes(input_df, logger):
    return go_cazyme_df

# "Gene ontology (molecular function)": ["NA"],
# "Gene ontology (biological process)": ["NA"],