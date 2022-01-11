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
        "prediction_dir",
        type=Path,
        help="Path to directory containing outputs from prediction tools and test sets fasta files",
    )
    parser.add_argument(
        "testset_dir",
        type=Path,
        help=(
            "Path to directory containing test set dir from create_test_sets_*.py,\n"
            "containing the dirs 'alignment_scores' and 'test_sets'"
        ),
    )
    parser.add_argument(
        "cazy",
        type=Path,
        help="Path to JSON file or local CAZyme db of CAZy annotations of proteins (ground truths)",
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


    # Add option to change the number rounds of bootstrapping
    parser.add_argument(
        "-b",
        "--bs_resampling",
        type=int,
        default=100,
        help="Number of rounds of bootstrap resampling",
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
        metavar="output directory",
        default=Path(os.getcwd()),
        help="Directory to which all outputs are written",
    )

    parser.add_argument(
        "-r",
        "--recombine_tools",
        type=Path,
        default=None,
        help="Path to txt file containing list of recombined tool groups. Three tools per row, separted by a single comma.",
    )
    # Add option to change the number of test sets used for bootstrapping
    parser.add_argument(
        "-s",
        "--bs_sample_size",
        type=int,
        default=6,
        help="Number of test sets used for bootstrapping CAZyme/non-CAZyme predictions",
    )
    # add option for evaluating performance per taxonomy group
    parser.add_argument(
        "-t",
        "--tax_groups",
        type=Path,
        default=None,
        help="Path to YAML file containing taxonomy groups of the source genomes",
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
