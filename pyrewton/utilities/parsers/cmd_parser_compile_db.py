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
"""Build parser for 'calculate_stats.py'"""

import argparse
import os

from pathlib import Path
from typing import List, Optional


def build_parser(argv: Optional[List] = None):
    """Return ArgumentParser parser for script."""
    # Create parser object
    parser = argparse.ArgumentParser(
        prog="calculate_stats.py",
        description="Programme to evaluate third-party CAZyme prediction tools",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    # Add positional arguments to parser

    # Add path to input fasta files
    parser.add_argument(
        "protein_dir",
        type=Path,
        help="Path to directory containing FASTA files of protein seqs extract from genomic assemblies",
    )
    parser.add_argument(
        "dbcan_dir",
        type=Path,
        help="Path to directory containing output from dbCAN, with one child per genome",
    )

    # Add optional arguments to parser

    parser.add_argument(
        "-c",
        "--cazy",
        type=Path,
        default=None,
        help=(
            "Path to local CAZyme db of CAZy annotations of proteins created using cazy_webscraper\n"
            "Will add CAZy annotations to the output."
        ),
    )

    parser.add_argument(
        "--cazy_date",
        type=str,
        default=None,
        help="Date CAZy data was pulled down",
    )

    parser.add_argument(
        "-d",
        "--output_db",
        type=Path,
        default=None,
        help="Name of resulting database",
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

    # Add option to specify output directory
    # This will enable creation of a new output directory
    parser.add_argument(
        "-o",
        "--output_dir",
        type=Path,
        metavar="output directory",
        default=Path(os.getcwd()),
        help="Directory to which all outputs are written. Defauly, write to cwd",
    )

    parser.add_argument(
        "--sql_echo",
        dest="sql_echo",
        action="store_true",
        default=False,
        help="Set SQLite echo to True, adds verbose SQL messaging",
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

    if argv is None:
        # parse command-line
        return parser
    else:
        # return namespace
        return parser.parse_args(argv)
