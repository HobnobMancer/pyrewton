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
"""Retrieves datasets for evaluating the CAZyme prediction tools: dbCAN, CUPP and eCAMI.

:cmd_args   :

:func    :

Writes out a FASTA per candidate species, containing all the protein sequences to analysed by the
CAZymes prediction tools. The datasets contain an equal number of CAZymes to non-CAZymes.
"""

import gzip
import logging
import re
import sys

from pathlib import Path
from typing import List, Optional

import pandas as pd

from pyrewton.utilities import build_logger
from pyrewton.utilities.cmd_get_evaluation_dataset import build_parser
from pyrewton.utilities.file_io import make_output_directory


def main(argv: Optional[List[str]] = None, logger: Optional[logging.Logger] = None):
    """Coordinate preparation for script, and terminating the programme when finished."""
    # programme preparation

    # build cmd-line arguments parser
    # Check if namepsace isn't passed, if not parse command-line
    if argv is None:
        # Parse command-line
        parser = build_parser()
        args = parser.parse_args()
    else:
        parser = build_parser(argv)
        args = parser.parse_args()
    
    # build logger
    # Note: log file is only created if specified at cmd-line
    if logger is None:
        logger = build_logger("get_ncbi_genomes", args)
    logger.info("Run initated")

    # retrieve paths to FASTA files
    fasta_files = get_fasta_file_paths()

    # create dataset per FASTA file
    for fasta in fasta_files:


    # terminate


def get_fasta_file_paths(args, logger):
    """Retrieve paths to call FASTA files in input dir.
    :param args: parser object
    :param logger: logger object
    Returns list of paths to fasta files.
    """
    # create empty list to store the file entries, to allow checking if no files returned
    fasta_file_paths = []

    # retrieve all files from input directory
    files_in_entries = (
        entry for entry in Path(args.input).iterdir() if entry.is_file()
    )
    # retrieve paths to fasta files
    for item in files_in_entries:
        # search for fasta file extension
        if item.name.endswith(".fasta") or item.name.endswith(".fa"):
            fasta_file_paths.append(item)

    # check files were retrieved from input directory
    if len(fasta_file_paths) == 0:
        logger.warning(
            (
                "No FASTA files retrieved from input directory.\n"
                "Check the path to the input directory is correct.\n"
                "Terminanting program."
            )
        )
        sys.exit(1)

    return fasta_file_paths


def build_protein_dataframe():
    """Build a dataframe containing the protein data for the current working input FASTA file."""
    return


def get_cazy_classification():
    """For each protein check if classified as a CAZyme by CAZy, and its annotated CAZy families."""
    return


def get_dataset():
    """Retrieve dataset of equal number of CAZymes and non-CAZymes."""
    return


def write_out_dataset():
    """Write out the dataset of CAZymes and non-CAZymes to a single FASTA file."""
    return
