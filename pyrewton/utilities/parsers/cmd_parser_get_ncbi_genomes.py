#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author:
# Emma E. M. Hobbs

# Contact
# eemh1@st-andrews.ac.uk

# Emma E. M. Hobbs,
# Biomolecular Sciences Building,
# University of St Andrews,
# North Haugh Campus,
# St Andrews,
# KY16 9ST
# Scotland,
# UK

# The MIT License
"""Build CLI for invoke_dbcan.py"""


import multiprocessing
import sys

from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser, _SubParsersAction
from pathlib import Path
from typing import List, Optional

from pyrewton.genbank.get_ncbi_genomes import get_ncbi_genomes


def build_parser(
    subps: _SubParsersAction, parents: Optional[List[ArgumentParser]] = None
) -> None:
    """Return ArgumentParser parser for script."""
    parser = subps.add_parser(
        "download_genomes", formatter_class=ArgumentDefaultsHelpFormatter
    )

    # Add positional arguments to parser
    # Add user email for Entrez
    parser.add_argument(
        "user_email",
        type=str,
        metavar="user email address",
        default=None,
        help="Email address of user, this must be provided for Entrez",
    )

    # Add optional arguments to parser

    # Add dataframe write out options
    # if not given dataframe written to STDOUT
    parser.add_argument(
        "-d",
        "--dataframe",
        type=Path,
        default=sys.stdout,
        help="Location of file for species table to be written to, include .csv extention",
    )
    # Add option to force file over writting
    parser.add_argument(
        "-f",
        "--force",
        dest="force",
        action="store_true",
        default=False,
        help="Force file over writting",
    )
    # Add option to disable pull down of GenBank files
    parser.add_argument(
        "-g",
        "--genbank",
        dest="genbank",
        action="store_true",
        default=True,
        help="Disable pulldown of GenBank (.gbff) files",
    )
    # Add input file name option
    # If not given input will be taken from STDIN
    parser.add_argument(
        "-i",
        "--input_file",
        type=Path,
        metavar="input file name",
        default=sys.stdin,
        help="Input filename",
    )
    # Add log file name option
    # If not given, no log file will be written out
    parser.add_argument(
        "-l",
        "--log",
        type=Path,
        metavar="log file name",
        default=None,
        help="Defines log file name and/or path",
    )
    # Add option to prevent over writing of existing files
    # and cause addition of files to output directory
    parser.add_argument(
        "-n",
        "--nodelete",
        dest="nodelete",
        action="store_true",
        default=False,
        help="enable/disable deletion of exisiting files",
    )
    # Add output file name option
    # Must include file extension
    # If not given, output will be written to STDOUT
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        metavar="output file name",
        default=sys.stdout,
        help="Output filename",
    )
    # Add custom maximum number of retries if network error is encountered
    parser.add_argument(
        "-r",
        "--retries",
        type=int,
        metavar="maximum number of retries",
        default=10,
        help="Defines the maximum number of retries if network errors are encountered",
    )
    # Add option to alter timeout allowance before cancelling downloading of files
    parser.add_argument(
        "-t",
        "--timeout",
        dest="timeout",
        action="store",
        default=10,
        help="Timeout for URL connections, in seconds",
    )
    # Add option for more detail (verbose) logging
    parser.add_argument(
        "-v",
        "--verbose",
        dest="verbose",
        action="store_true",
        default=False,
        help="Set logger level to 'INFO'",
    )

    parser.set_defaults(func=get_ncbi_genomes.main)
