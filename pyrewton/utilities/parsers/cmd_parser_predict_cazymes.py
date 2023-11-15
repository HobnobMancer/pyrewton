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


import multiprocessing
import sys

from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser, _SubParsersAction, Action
from pathlib import Path
from typing import List, Optional

from pyrewton.cazymes.annotate_cazome import predict_cazymes


class ValidateTools(Action):
    """Check the user has provided valid cazyme classifiertool names."""
    def __call__(self, parser, args, values, option_string=None):
        valid_formats = ("cupp", "ecami", "dbcan")
        invalid = False
        for value in values:
            if value.lower() not in valid_formats:
                invalid = True
                raise ValueError(f'Invalid tool name "{value}" provided. Accepted tools: {valid_formats} (case insensitive)')
        if invalid:
            sys.exit(1)
        setattr(args, self.dest, values)


def build_parser(
    subps: _SubParsersAction, parents: Optional[List[ArgumentParser]] = None
) -> None:
    """Return ArgumentParser parser for script."""
    parser = subps.add_parser(
        "predict_cazymes", formatter_class=ArgumentDefaultsHelpFormatter
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
        "tools",
        nargs='+',
        action=ValidateTools,
        choices=["cupp","ecami","dbcan"],
        type=str,
        help="CAZyme classifiers to run. Space separted list. Pick as many as wanted. Case insensitive",
    )
    parser.add_argument(
        "tool_dir",
        type=Path,
        help=(
            "Path to parent directory where CAZyme classifiers are installed.\n"
            "E.g. point to the parent directory containing the directory called 'dbcan/'"
        ),
    )

    # Add optional arguments to parser
    
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
        metavar="output directory",
        default=sys.stdout,
        help="Directory to which all outputs are written",
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

    parser.set_defaults(func=predict_cazymes.main)
