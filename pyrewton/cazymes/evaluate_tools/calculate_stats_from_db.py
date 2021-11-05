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

:class TestSet: represents a single test parsed by the prediction tools

:func main: co-ordinate entire module function
:func get_cazy_dict: retrieve dictionary of CAZy classifications of proteins
:func get_predictions: retrieve the paths to the directories contain the prediction tool predictions

Creates dataframes of CAZyme predictions and report
summarising statistical analsis of prediction accuracy.
"""

import logging
import re
import sys

from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import List, Optional

from scraper.sql.sql_orm import get_db_connection
from saintBioutils.utilities.logger import config_logger
from saintBioutils.utilities.file_io import make_output_directory
from tqdm import tqdm

from pyrewton.cazymes.evaluate_tools import stats
from pyrewton.utilities.parsers.cmd_parser_calc_stats import build_parser


@dataclass
class TestSet:
    """TestSet from a genomic assembly.

    :param fasta: path, path to query fasta file
    :param tax_id: str, NCBI taxonomy ID of host species
    :param source: str, source of proteins within fasta file
    :param prediction_paths: path, path to which prediction tool output is written
    """

    fasta: Path  # path to FASTA file containing proteins for prediction
    tax_id: str  # NCBI taxonomy id, prefix NCBI:txid
    source: str  # source of protein sequences, genomic assembly or database
    prediction_paths: dict()  # contains path to outputs

    def __str__(self):
        """Representation of object"""
        return(
            (
                f"<tax-id={self.tax_id} source={self.source} "
                f"fasta={self.fasta} prediction_paths={self.prediction_paths}"
            )
        )


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

    # If specified output directory for genomic files, create output directory
    if args.output is not sys.stdout:
        make_output_directory(args.output, args.force, args.nodelete)

    # open the CAZy dict
    db_connection = get_db_connection(args.cazy, args, False)

    # retrieve paths to all dirs
    predictions = get_predictions(args.input)

    # perform stats evaluations
    stats.evaluate_performance(predictions, db_connection, 'db', args)

    time_stamp = datetime.now().strftime("%Y_%m_%d")
    stats.get_fam_freq(args, db_connection, time_stamp, 'db')  # USED IN R EVALUATION


def get_predictions(prediction_dir):
    """Retrieve create class instances to represent output for each test set.

    :param prediction_dir: path to dir containing all prediction tool outputs.

    Return list containing class instances, one instance per test set.
    """
    logger = logging.getLogger(__name__)

    output_dirs = [f for f in prediction_dir.iterdir() if f.is_dir()]  # store the paths to the output dirs

    predictions = []  # store TestSet instances

    for output_dir in tqdm(output_dirs, desc="Building TestSet instances"):
        # get the genomic assembly accession
        search_result = re.search(r"GCA_\d+?(_|.)\d+?\.\d+?.", str(output_dir), re.IGNORECASE)
        try:
            genomic_accession = search_result.group()[:-1]
        except AttributeError as e:
            logger.error(
                f"Could not get genomic asseccion from {output_dir}.\n"
                "Not including test set in the evaluation. The following error was raised:\n"
                f"{e}"
            )
            continue

        # get the taxonomy id
        search_result = re.search(r"txid\d+?_", str(output_dir), re.IGNORECASE)
        try:
            tax_id = search_result.group()[:-1]
        except AttributeError as e:
            logger.error(
                f"Could not get taxonomy ID from {output_dir}.\n"
                "Not including test set in the evaluation. The following error was raised:\n"
                f"{e}"
            )
            continue

        # build path to the FASTA file
        fasta_path = f"genbank_proteins_{tax_id}_{genomic_accession}._test_set.fasta"
        fasta_path = prediction_dir / fasta_path

        # build dict of output paths
        output_paths = {}
        output_paths["dir"] = output_dir

        # build TestSet instance
        test_set = TestSet(fasta_path, tax_id, genomic_accession, output_paths)
        predictions.append(test_set)

    return predictions


if __name__ == "__main__":
    main()