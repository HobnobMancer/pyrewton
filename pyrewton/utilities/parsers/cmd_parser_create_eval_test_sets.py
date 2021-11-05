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
"""Build parser for 'create_evaluation_tests_set_from_*.py'"""


import argparse

from pathlib import Path
from typing import List, Optional


def build_parser_dict(argv: Optional[List] = None):
    """Return ArgumentParser parser for script."""
    # Create parser object
    parser = argparse.ArgumentParser(
        prog="create_test_set_from_dict.py",
        description="Build test sets for evaluating CAZyme prediction tools",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # Add positional arguments to parser

    parser.add_argument(
        "email",
        type=str,
        help="User email address, required by NCBI to use Entrez",
    )

    parser.add_argument(
        "yaml",
        type=Path,
        help="Path to yaml file containing genomic assembly accessions",
    )

    parser.add_argument(
        "cazy",
        type=Path,
        help="Path to JSON file containing dict of GenBank accesions and CAZy familes from CAZy",
    )

    parser.add_argument(
        "output",
        type=Path,
        help="Path to output directory",
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
    parser.add_argument(
        "-g",
        "--genomes",
        type=Path,
        default=None,
        help="Path to dir containing already downloaded genomes",
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
    # Add option to change the sample size (this is the number of CAZymes so the total sample size will be 2x sample_size)
    parser.add_argument(
        "-s",
        "--sample_size",
        type=int,
        default=100,
        help=(
            "Number of CAZymes to be included in the sample size."
            "Total sample size is twice this, becuase the test set "
            "includes equal numbers of CAZymes and non-CAZymes"
        ),
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


def build_parser_db(argv: Optional[List] = None):
    """Return ArgumentParser parser for script."""
    # Create parser object
    parser = argparse.ArgumentParser(
        prog="create_test_set_from_db.py",
        description="Build test sets for evaluating CAZyme prediction tools",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # Add positional arguments to parser

    parser.add_argument(
        "email",
        type=str,
        help="User email address, required by NCBI to use Entrez",
    )

    parser.add_argument(
        "yaml",
        type=Path,
        help="Path to yaml file containing genomic assembly accessions",
    )

    parser.add_argument(
        "cazy",
        type=Path,
        help="Path to local CAZyme database",
    )

    parser.add_argument(
        "output",
        type=Path,
        help="Path to output directory",
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
    parser.add_argument(
        "-g",
        "--genomes",
        type=Path,
        default=None,
        help="Path to dir containing already downloaded genomes",
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
    # Add option to change the sample size (this is the number of CAZymes so the total sample size will be 2x sample_size)
    parser.add_argument(
        "-s",
        "--sample_size",
        type=int,
        default=100,
        help=(
            "Number of CAZymes to be included in the sample size."
            "Total sample size is twice this, becuase the test set "
            "includes equal numbers of CAZymes and non-CAZymes"
        ),
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
