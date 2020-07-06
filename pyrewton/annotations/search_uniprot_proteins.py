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

