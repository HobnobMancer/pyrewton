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
"""Build parser for 'process_cazyme_predictions.py'"""

import argparse
import sys

from pathlib import Path
from typing import List, Optional


def build_parser(argv: Optional[List] = None):
    """Return ArgumentParser parser for script."""
    # Create parser object
    parser = argparse.ArgumentParser(
        prog="process_cazyme_predictions.py",
        description="Programme to formate dbCAN, CUPP and eCAMI output into dataframes",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # Add positional arguments to parser
    parser.add_argument(
        "input_df",
        type=Path,
        metavar="input dataframe",
        help="Path to the input dataframe",
    )

    # Add path to the directory containing the output file of the CAZyme prediction tool
    parser.add_argument(
        "input", type=Path, metavar="input file name", help="Path to the CAZyme"
    )

    # Add string to define the CAZyme prediction tool whose output is to be processed
    parser.add_argument(
        "tool",
        choices=["dbcan", "cupp", "ecami"],
        help="The CAZyme prediction tool whose output is to be processed",
    )

    # Add additional optional arguments
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
        "--output",
        type=Path,
        metavar="output directory name",
        default=sys.stdout,
        help="Directory output files are written to",
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
