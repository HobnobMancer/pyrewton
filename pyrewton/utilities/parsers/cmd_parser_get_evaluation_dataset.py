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
"""Build parser for 'get_evaluation_dataset.py'"""

import argparse
import sys

from pathlib import Path
from typing import List, Optional


def build_parser(argv: Optional[List] = None):
    """Return ArgumentParser parser for script."""
    # Create parser object
    parser = argparse.ArgumentParser(
        prog="get_evaluation_dataset.py",
        description="Build test sets for evaluating CAZyme prediction tools",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # Add positional arguments to parser

    parser.add_argument(
        "input",
        type=Path,
        help="Path to directory containing FASTA files of proteomes",
    )

    parser.add_argument(
        "database",
        type=Path,
        help="Path to local CAZy database",
    )

    # Add optional arguments to parser

    # Add option to change test set sample size
    parser.add_argument(
        "-s",
        "--sample",
        type=int,
        default=400,
        help=(
            "Total sample size of each test set.\n"
            "50-percent of the test set is positive controls, 50-percent negative controls"
        ),

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
    # Add option to specific directory for log to be written out to
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
    # Add option to specify directory for output to be written to
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        metavar="output directory for fasta files",
        default=sys.stdout,
        help="Path to directory yo whicg FASTA files are written",
    )
    # Add option to specify verbose logging
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
