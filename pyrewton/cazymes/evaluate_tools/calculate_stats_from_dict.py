#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) University of St Andrews 2020-2021
# (c) University of Strathclyde 2020-2021
# (c) James Hutton Institute 2020-2021
#
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
"""Predict CAZymes from protein sequences and evaluate prediction accuracy.

:cmd_args input: directory containing all prediction output to be parse
:cmd --beta: change the value of beta used when calculating Fbeta-score (default = 1)
:cmd --force: force writing out to existing output directory
:cmd --log: path to file to write out log to
:cmd --nodelete: do not delete data in an existing output directory
:cmd --output: path to output directory
:cmd --bs_sample_size: number of test sets included in bootstrap evaluation
:cmd --bs_resampling: number of times to perform bootstrap resampling
:cmd --verbose: enable verbose logging

Creates dataframes of CAZyme predictions and report
summarising statistical analsis of prediction accuracy.
"""

import logging
import os

from datetime import datetime
from pathlib import Path
from saintBioutils.utilities.file_io import make_output_directory
from saintBioutils.utilities.logger import config_logger
from typing import List, Optional

from pyrewton.cazymes.evaluate_tools import stats
from pyrewton.utilities.file_io import io_create_eval_testsets
from pyrewton.utilities.parsers.cmd_parser_calc_stats import build_parser


def main(argv: Optional[List[str]] = None, logger: Optional[logging.Logger] = None):
    """Set up parser and logger, and coordinate prediction of CAZymes."""
    if argv is None:
        # Parse command-line
        parser = build_parser()
        args = parser.parse_args()
    else:
        parser = build_parser(argv)
        args = parser.parse_args()

    # Initiate logger
    # Note: log file only created if specified at cmdline
    if logger is None:
        logger = logging.getLogger(__name__)
        config_logger(args)

    # class_classification_evaluation_<date>.csv from class_predictions_classifications_<date>.csv

    # If specified output directory for output files, create output directory
    if args.output != Path(os.getcwd()):
        make_output_directory(args.output, args.force, args.nodelete)
    
    # make dir to store binary CAZyme classification csv files
    binary_classification_dir = args.output / "binary_classifications"
    make_output_directory(binary_classification_dir, args.force, args.nodelete)

    # open the CAZy dict
    cazy_dict = io_create_eval_testsets.get_cazy_dict(args.cazy)

    # retrieve paths to all dirs
    predictions = stats.get_predictions(args)

    # perform stats evaluations
    ground_truths_fam_df = stats.evaluate_performance(predictions, cazy_dict, 'dict', args)

    time_stamp = datetime.now().strftime("%Y_%m_%d")
    stats.get_fam_freq(args, cazy_dict, time_stamp, 'dict', ground_truths_fam_df)  # USED IN R EVALUATION


if __name__ == "__main__":
    main()