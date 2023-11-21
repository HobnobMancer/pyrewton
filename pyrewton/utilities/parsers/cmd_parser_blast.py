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
"""Build parser for running all versus all diamond or blast"""


import sys

from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser, _SubParsersAction, Action
from pathlib import Path
from typing import List, Optional

from pyrewton.select_candidates import run_all_vs_all_blast


class ValidateMethod(Action):
    """Check the method is ok."""
    def __call__(self, parser, args, value, option_string=None):
        invalid = False
        valid_methods = ['BLAST', 'DIAMOND']
        if value not in valid_methods:
            invalid = True
            raise ValueError(f'{value} is not a valid method. Valid methods: {valid_methods}')
        if invalid:
            sys.exit(1)
        setattr(args, self.dest, value.upper())


def build_parser(
    subps: _SubParsersAction, parents: Optional[List[ArgumentParser]] = None
) -> None:
    """Return ArgumentParser parser for script."""
    parser = subps.add_parser(
        "run_all_v_all", formatter_class=ArgumentDefaultsHelpFormatter
    )

    # Add positional arguments to parser

    parser.add_argument(
        "method",
        type=Path,
        action=ValidateMethod,
        help=(
            "Method for all versus all sequence alignemnt. BLAST or DIAMOND. Not case sensitivity"
        ),
    )

    # Add optional arguments to parser
    parser.add_argument(
        "--db_path",
        type=Path,
        default="./diamond.db",
        help="Path to build DIAMOND database",
    )

    # Add optional arguments to parser
    parser.add_argument(
        "--evalue",
        type=float,
        default=10.0,
        help="E-value threshold",
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

    # Add option for more detail (verbose) logging
    parser.add_argument(
        "-v",
        "--verbose",
        dest="verbose",
        action="store_true",
        default=False,
        help="Set logger level to 'INFO'",
    )

    parser.set_defaults(func=run_all_vs_all_blast.main)
