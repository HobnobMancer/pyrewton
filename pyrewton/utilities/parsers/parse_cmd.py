#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) University of St Andrews 2022
# (c) University of Strathclyde 2022
# (c) James Hutton Institute 2022
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
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
"""Build CLI for pyrewton"""


from asyncio import subprocess
import sys

from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter, Namespace
from typing import List, Optional

from pathlib import Path
from typing import List, Optional

from pyrewton import __version__, __citation__
from pyrewton.utilities.parsers import (
    cmd_parser_get_ncbi_genomes,
    cmd_parser_get_genbank_annotations,
    cmd_parser_predict_cazymes,
    cmd_parser_compile_db,
    cmd_parser_uniprot,
    cmd_parser_extract_db_seqs,
    cmd_parser_gather_seqs,
    cmd_parser_cluster_seqs,
    cmd_parser_get_cluster_seqs,
)


def build_parser(argv: Optional[List] = None) -> Namespace:
    """Parse command-line arguments for script.

    :param argv: Namesapce, command-line arguments

    The script offers a single main, parser with subcommands for the actions.
    """
    # Create parser object
    parser_main = ArgumentParser(
        prog="pyrewton",
        description="Compile and explore CAZymes and CAZomes",
        formatter_class=ArgumentDefaultsHelpFormatter,
    )
    subparsers = parser_main.add_subparsers(
        title="subcommands", description="Valid subcommands",
    )

    # add args to the main parser
    parser_main.add_argument(
        "--version",
        action="store_true",
        default=False,
        help="Print version number"
    )
    parser_main.add_argument(
        "--citation",
        action="store_true",
        default=False,
        help="Print citation information"
    )

    # Download genomes and get protein seqs
    cmd_parser_get_ncbi_genomes.build_parser(subparsers)
    
    # Extract protein seqs from .gbff.gz files
    cmd_parser_get_genbank_annotations.build_parser(subparsers)

    # Automate running CAZyme classifiers over multi-sequence FASTA files in an input dir
    cmd_parser_predict_cazymes.build_parser(subparsers)

    # build and add data to a local CAZome database
    cmd_parser_compile_db.build_parser(subparsers)

    cmd_parser_uniprot.build_parser(subparsers)

    cmd_parser_extract_db_seqs.build_parser(subparsers)

    cmd_parser_gather_seqs.build_parser(subparsers)

    cmd_parser_cluster_seqs.build_parser(subparsers)

    cmd_parser_get_cluster_seqs.build_parser(subparsers)

    # Parse arguments
    # The list comprehension is to allow PosixPaths to be defined and passed in testing
    if argv is None:
        if len(sys.argv) == 1:
            argv = ["-h"]
        else:
            argv = sys.argv[1:]
    return parser_main.parse_args([str(_) for _ in argv])
