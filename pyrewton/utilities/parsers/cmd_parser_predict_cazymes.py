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
"""Build parser for 'predict_cazymes.py'"""

import argparse
import sys

from pathlib import Path
from typing import List, Optional


def build_parser(argv: Optional[List] = None):
    """Return ArgumentParser parser for script."""
    # Create parser object
    parser = argparse.ArgumentParser(
        prog="predict_cazymes.py",
        description=(
            "Programme to invoke and evaluate third-party CAZyme prediction tools\n"
            "To be able to access dbCAN, CUPP and eCAMI please invokve script when cwd is\n"
            "pyrewton/cazymes/evaluation within the pyrewton programm."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # Add positional arguments to parser

    # Add path to input fasta files
    parser.add_argument(
        "input",
        type=Path,
        metavar="input directory",
        help="Path to directory containing FASTA files for prediction tools",
    )
    parser.add_argument(
        "output",
        type=Path,
        metavar="output directory",
        help="Directory to which all outputs are written",
    )

    # Add optional arguments to parser
    # Add option to define F-beta weighting
    parser.add_argument(
        "-b",
        "--beta",
        dest="beta",
        type=int,
        default=1,
        help="Weighting of Fbeta-stat when calculating Fbeta score of CAZy family prediction",
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
