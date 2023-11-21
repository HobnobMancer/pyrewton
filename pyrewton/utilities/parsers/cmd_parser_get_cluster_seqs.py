#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) University of St Andrews 2020-2021
# (c) University of Strathclyde 2020-2021
# Author:
# Emma E. M. Hobbs
#
# Contact
# eemh1@st-andrews.ac.uk
#
# Emma E. M. Hobbs,
# Biomolecular Sciences Building,
# University of St Andrews,
# North Haugh Campus,
# St Andrews,
# KY16 9ST
# Scotland,
# UK
#
# The MIT License
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#
"""Build parser for 'add_uniprot_data.py'"""


import sys

from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser, _SubParsersAction, Action
from pathlib import Path
from typing import List, Optional

from pyrewton.uniprot import get_uniprot_proteins

def build_parser(
    subps: _SubParsersAction, parents: Optional[List[ArgumentParser]] = None
) -> None:
    """Return ArgumentParser parser for script."""
    parser = subps.add_parser(
        "cluster_cazymes", formatter_class=ArgumentDefaultsHelpFormatter
    )

    # Add positional arguments to parser

    # Add path to input fasta files
    parser.add_argument(
        "fasta",
        type=Path,
        help="Path to FASTA file of protein sequence to cluster",
    )

    # Add optional arguments to parser
    parser.add_argument(
        "--mmseqs_db",
        type=Path,
        default="mmseqs/mmseqs_db",
        help="Path to write out mmseqs db. Default ./mmseqs/mmseqs_db",
    )
    parser.add_argument(
        "--mmseqs_out",
        type=Path,
        default="mmseqs/mmseqs_out",
        help="Path to write out mmseqs output file. Default ./mmseqs/mmseqs_out",
    )
    parser.add_argument(
        "--out_tsv",
        type=Path,
        default="./mmseqs_out.tsv",
        help="Path to write out mmseqs TSV file with cluster inforamtion. Default ./mmseqs_out.tsv",
    )
    parser.add_argument(
        "--pident",
        type=float,
        default=0.7,
        help="Percentage identity cutoff (as DECIMAL)",
    )
    parser.add_argument(
        "--cov",
        type=float,
        default=0.7,
        help="Converage cutoff (as DECIMAL)",
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

    parser.set_defaults(func=get_uniprot_proteins.main)
