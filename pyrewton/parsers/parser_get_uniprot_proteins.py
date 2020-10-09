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
"""Build parser for 'get_uniprot_proteins.py'"""

import argparse
import sys

from pathlib import Path
from typing import List, Optional


def build_parser(argv: Optional[List] = None):
    """Return ArgumentParser parser for script."""
    # Create parser object
    parser = argparse.ArgumentParser(
        prog="get_uniprot_proteins.py",
        description="Retrieve protein data from UniProtKB",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # Add positional arguments to parser
    # Specify input config file
    parser.add_argument(
        "input",
        type=Path,
        metavar="configuration file",
        help="Path to configuration file",
    )

    # Add optional arguments to parser
    # Add option to enable/disable fasta file writing
    parser.add_argument(
        "-fa",
        "--fasta",
        dest="fasta",
        action="store_true",
        default=False,
        help="Enable fasta file writing",
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
        metavar="output directory path",
        default=sys.stdout,
        help="output directory path",
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
