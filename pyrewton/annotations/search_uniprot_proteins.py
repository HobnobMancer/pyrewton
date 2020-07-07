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
and by GO (Gene Ontology) annotated function inferring cazyme functionality.

:func main: coordianate searching for cazymes in local UniProt dataframe.
:func get_cazy_proteins: retrieve entries from input dataframe with link to CAZy DB
:func get_ec_cazymes: retrieve entries from input datafame with EC
    number(s) indicated cazyme functionality
:func get_go_cazymes: retrieve entries from input datafame with GO annotated function(s)
    numbers indicated cazyme functionality
:func retrieve_df_subset: coordinate retrieval of rows from input dataframe

Writes out 4 dataframes to csv files:
[1] cazy_linked_cazymes_df
[2] ec_num_only_cazymes
[3] go_fun_only_cazymes
[4] ec_go_cazymes
"""

import io
import logging
import re
import sys

from typing import List, Optional

import pandas as pd

from tqdm import tqdm

from pyrewton import file_io
from pyrewton.loggers import build_logger
from pyrewton.parsers.parser_get_uniprot_proteins import build_parser


def main(argv: Optional[List[str]] = None, logger: Optional[logging.Logger] = None):
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
    logger.info("Opening input dataframe %s", args.df_input)
    input_df = pd.read_csv(args.df_input, header=0, index_col=0)

    # Separate entries with and without CAZy links
    cazyme_dfs = get_cazyme_subset_df(input_df, logger)
    # Write out cazy_linked df
    file_io.write_out_pre_named_dataframe(
        cazyme_dfs[0],
        "cazy_linked_cazymes_df",
        logger,
        args.ouput,
        args.force,
        args.nodelete,
    )

    # Compare dataframes (CAZy, EC number and GO function) to retrieve dataframes of
    # (1) only EC number, (2) only GO function and (3) GO+EC inferred cazyme functionality
    ec_go_cazyme_dfs = compare_cazyme_dfs(cazyme_dfs, logger)

    # Write out dfs to .csv files
    file_names = ["go_fun_only_cazymes", "ec_num_only_cazymes", "ec_go_cazymes"]
    fname_index = 0
    for df in ec_go_cazyme_dfs:
        file_io.write_out_pre_named_dataframe(
            df, file_names[fname_index], logger, args.output, args.force, args.nodelete
        )
        fname_index += 1

    logger.info("Program finished.")


def get_cazyme_subset_df(input_df, logger):
    """Coordinate retrieval of entries with indicated cazyme functionality.

    :param input_df: pandas dataframe
    :param logger: logger object

    Return 3 pandas dataframes.
    """
    # Retrieve UniProt entries with link to CAZy database
    cazy_linked_df = get_cazy_cazymes(input_df, logger)

    # Search non-CAZy linked entries to retrieve UniProt entries whose EC
    # number inferred cazyme functionality
    ec_inferred_df = get_ec_cazymes(cazy_linked_df[1], logger)

    # Search non-CAZy linked entries to retrieve UniProt entries whose GO
    # function annotated inferred cazyme functionality
    go_inferred_df = get_go_cazymes(cazy_linked_df[1], logger)

    # Return entries with CAZy link, and those with EC number or GO function
    # inferred cazyme functionality
    return cazy_linked_df[0], ec_inferred_df, go_inferred_df


def get_cazy_cazymes(input_df, logger):
    """Separate entries into those with and without a link to the CAZy database.

    :func input_df: pandas dataframe
    :func logger: logger, object

    Return 2 pandas dataframes: with link and without links to CAZy database.
    """
    # retrieve indexes of entries with link to CAZy database
    logger.info("Retrieving entries with link to CAZy database")
    cazy_link_indexes = input_df["UniProt linked protein families"].str.contains(
        r"cazy", flags=re.IGNORECASE, regex=True, na=False
    )
    return input_df[cazy_link_indexes], input_df[~cazy_link_indexes]


def get_ec_cazymes(non_cazy_input_df, logger):
    """Retrieve subset of entries with indicated cazyme functionality
    from EC number(s), and no CAZy database link.

    :param non_cazy_input_df: pandas dataframe of entries from input df with
        no link to CAZy database
    :param logger: logger object

    Return pandas dataframe.
    """
    logger.info("Retrieving rows whose EC number indicates cazyme functionality")
    ec_series = []  # store indexed pandas series results of EC search results
    search_terms = [
        "3.1.1.11",
        "3.1.1.72",
        "3.1.1.73",
        "3.2.1.4",
        "3.2.1.6",
        "3.2.1.7",
        "3.2.1.8",
        "3.2.1.15",
        "3.2.1.21",
        "3.2.1.25",
        "3.2.1.37",
        "3.2.1.40",
        "3.2.1.55",
        "3.2.1.67",
        "3.2.1.91",
        "3.2.1.131",
        "3.2.1.139",
        "3.2.1.156",
    ]

    for term in tqdm((search_terms), desc="Searching EC numbers"):
        ec_series.append(
            non_cazy_input_df["EC number"].str.contains(
                rf"{term}", flags=re.IGNORECASE, regex=True, na=False
            )
        )

    # Retrieve search corresponding rows from input dataframe using search results
    return retrieve_df_subset(non_cazy_input_df, ec_series, logger)


def get_go_cazymes(non_cazy_input_df, logger):
    """Retrieve subset of entries with indicated cazyme functionality
    from GO annotated function(s).

    :param non_cazy_input_df: pandas dataframe of entries from input df with
        no link to CAZy database
    :param logger: logger object

    Return pandas dataframe.
    """
    logger.info("Retrieving rows whose EC number indicates cazyme functionality")
    go_series = []  # store indexed pandas series results of EC search results
    search_terms = [
        "arabino",
        "arabinofuranosidase",
        "cellulase",
        "feruloyl",
        "galacturonan",
        "glucanase",
        "glucosidase",
        "glucuronisdase",
        "glucuronoyl",
        "lichenase",
        "mannan",
        "pectin",
        "xylan",
        "xylo",
        "xylosidase",
    ]

    for term in tqdm((search_terms), desc="Searching GO (Gene Ontology) functions"):
        go_series.append(
            non_cazy_input_df["EC number"].str.contains(
                rf"{term}", flags=re.IGNORECASE, regex=True, na=False
            )
        )

    # Retrieve search corresponding rows from input dataframe using search results
    return retrieve_df_subset(non_cazy_input_df, go_series, logger)


def retrieve_df_subset(non_cazy_input_df, df_index_list, logger):
    """Retrieve subset of input dataframe, as a single dataframe.
    
    :param input_df: pandas df, original input dataframe
    :param df_index_list: list, results from searching input df using .str.contains.
    :param logger: logger object
    
    Return pandas dataframe.
    """
    # Create empty dataframe to add subset of dataframes to
    cazyme_subset_df = pd.DataFrame(
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
    logger.info("Using search results to retrieve rows from input dataframe.")
    for pandas_series in df_index_list:
        cazyme_subset_df = cazyme_subset_df.append(non_cazy_input_df[pandas_series])
    cazyme_subset_df = cazyme_subset_df[
        cazyme_subset_df.duplicated(subset="UniProt entry ID", keep="first")
    ]
    return cazyme_subset_df


def compare_cazyme_dfs(cazyme_dfs, logger):
    """"Identify common entries across EC number and GO function inferred cazymes dfs.
    
    :param cazyme_dfs: 3 pandas dataframes (CAZy, EC number and GO inferred cazymes)
    :param logger: logger object
    
    Return 3 pandas dataframes of entries with GO only, EC only and GO+EC inferred cazyme functionality.
    """
    # create single dataframe of all EC number and GO inferred cazyme functionality.
    non_cazy_cazyme_df = cazyme_dfs[1].append(cazyme_dfs[2])
    non_cazy_cazyme_df = non_cazy_cazyme_df[
        non_cazy_cazyme_df.duplicated(subset="UniProt entry ID", keep="first")
    ]

    # Retrieve entries with GO only inferred cazyme function
    ec_indexes = non_cazy_cazyme_df["EC Number"].isna()
    go_only_cazymes = non_cazy_cazyme_df[
        ~ec_indexes
    ]  # removes entries with no EC number
    # Retrieve entries with GO function inferred cazyme function
    go_indexes = non_cazy_cazyme_df["GO molecular function"].isna()
    ec_only_cazymes = non_cazy_cazyme_df[
        ~go_indexes
    ]  # removes entries with no GO number

    # Retreve entries with GO function AND EC number inferred cazyme function and no CAZy database link
    ec_go_cazymes = non_cazy_cazyme_df.append(ec_only_cazymes)
    ec_go_cazymes = ec_go_cazymes.append(go_only_cazymes)
    ec_go_cazymes = ec_go_cazymes[
        ec_go_cazymes.duplicated(subset="UniProt entry ID", keep="first")
    ]

    return go_only_cazymes, ec_only_cazymes, ec_go_cazymes


if __name__ == "__main__":
    main()
