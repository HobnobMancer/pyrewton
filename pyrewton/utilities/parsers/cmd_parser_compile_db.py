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
"""Build parser for 'calculate_stats.py'"""


import sys

from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser, _SubParsersAction, Action
from pathlib import Path
from typing import List, Optional

from pyrewton.cazymes.annotate_cazome import build_cazome_db


class ValidateProteinDir(Action):
    """Check the protein directory exists."""
    def __call__(self, parser, args, value, option_string=None):
        invalid = False
        if (value.exists() is False) or (value.is_dir() is False):
            raise ValueError(f'Could not find directory at "{value}". Please check the path is correct.')
        if invalid:
            sys.exit(1)
        setattr(args, self.dest, value)


class ValidateConfig(Action):
    """Check the config file exists."""
    def __call__(self, parser, args, value, option_string=None):
        invalid = False
        if (value.exists() is False):
            invalid = True
            raise ValueError(f'Could not find config file at "{value}". Please check the path is correct.')
        if invalid:
            sys.exit(1)
        setattr(args, self.dest, value)


class ValidateCazy(Action):
    """Check the cazy database file exists."""
    def __call__(self, parser, args, value, option_string=None):
        invalid = False
        if (value.exists() is False):
            invalid = True
            raise ValueError(f'Could not find local CAZyme (CAZy) database file at "{value}". Please check the path is correct.')
        if invalid:
            sys.exit(1)
        setattr(args, self.dest, value)


def build_parser(
    subps: _SubParsersAction, parents: Optional[List[ArgumentParser]] = None
) -> None:
    """Return ArgumentParser parser for script."""
    parser = subps.add_parser(
        "compile_cazome_db", formatter_class=ArgumentDefaultsHelpFormatter
    )

    parser.add_argument(
        "config_file",
        type=Path,
        action=ValidateConfig,
        help=(
            "Path to directory containing FASTA files of protein seqs extract from genomic assemblies"
        ),
    )

    # args to build new db or add data to an existing one 
    parser.add_argument(
        "--new_db",
        type=Path,
        default=None,
        help=(
            "Path to build a new CAZome database. Path to create database file.\n"
            "CANNOT be used at same time as --db"
        ),
    )
    parser.add_argument(
        "--db",
        type=Path,
        default=None,
        help=(
            "Path to existing CAZome database to add data to.\n"
            "CANNOT be used at same time as --db"
        ),
    )

    parser.add_argument(
        "--genome_csv",
        type=Path,
        action=ValidateConfig,
        help=(
            "Path to CSV file listing taxonomies and genomic accessions. Created using pyrewton download_genomes subcommand"
        ),
    )
    

    # Add optional arguments to parser
    parser.add_argument(
        "--protein_dir",
        type=Path,
        action=ValidateProteinDir,
        default=None,
        help=(
            "Path to directory containing FASTA files of protein seqs extract from genomic assemblies"
        ),
    )

    parser.add_argument(
        "--cazy",
        type=Path,
        default=None,
        action=ValidateCazy,
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

    parser.set_defaults(func=build_cazome_db.main)
