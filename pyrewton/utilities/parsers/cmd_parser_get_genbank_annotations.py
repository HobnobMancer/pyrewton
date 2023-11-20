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
"""Build parser for 'get_genbank_annotations.py'"""


import multiprocessing
import sys

from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser, _SubParsersAction
from pathlib import Path
from typing import List, Optional

from pyrewton.genbank.get_genbank_annotations import get_genbank_annotations


def build_parser(
    subps: _SubParsersAction, parents: Optional[List[ArgumentParser]] = None
) -> None:
    """Return ArgumentParser parser for script."""
    parser = subps.add_parser(
        "extract_protein_seqs", formatter_class=ArgumentDefaultsHelpFormatter
    )

    # Add positional arguments to parser
    # Specify path to input dataframe
    parser.add_argument(
        "input_df",
        type=Path,
        metavar="input dataframe name",
        help="Path to input dataframe",
    )
    # Specify path to directory containing GenBank files
    parser.add_argument(
        "genome_directory",
        type=Path,
        help="Path to directory containing compressed genomic assemblies in GenBank Flat File Format (.gbff.gz)",
    )

    # Add optional arguments to parser

    # Add option to specify path for output datafame to be written to
    parser.add_argument(
        "-d",
        "--output_df",
        type=Path,
        default=sys.stdout,
        help="path to output directory to write FASTA files to",

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
        default=sys.stdout,
        help="Output directory. Path to directory to which FASTA files are written",
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

    parser.set_defaults(func=get_genbank_annotations.main)
