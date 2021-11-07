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
"""Invokes CAZyme prediction tools to predict CAZymes from protein sequences in the test sets."""

import logging
import os
import re
import shutil
import sys

from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from saintBioutils.utilities.file_io import get_paths, make_output_directory
from saintBioutils.utilities.logger import config_logger
from tqdm import tqdm
from typing import List, Optional

from pyrewton.cazymes.evaluate_tools.tools import invoke_prediction_tools
from pyrewton.utilities.parsers.cmd_parser_predict_cazymes import build_parser


@dataclass
class Proteome:
    """Proteome from a genomic assembly.

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
    # build parser
    # Parse arguments
    # Check if namepsace isn't passed, if not parse command-line
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

    # check current working directory, to make sure can access the CAZyme prediction tools
    check_cwd()

    # If specified output directory for genomic files, create output directory
    make_output_directory(args.output, args.force, args.nodelete)

    # invoke prediction tools and build prediciton Proteome instances
    get_predictions(args)

    logger.info("Program finished, and no terminating.")


def check_cwd():
    """Check user invoked script in correct directory so can access the prediciton tools.

    If user is not within any of pyrewton/cazymes/evaluate_tools module directories then
    raise error and terminate program. Excepted directories for starting at:
    pyrewton/cazymes/evaluate_tools
    pyrewton/cazymes/evaluate_tools/tools

    Return nothing.
    """
    logger = logging.getLogger(__name__)
    current_path = os.getcwd()

    if current_path.endswith("pyrewton/cazymes/evaluate_tools"):
        return

    elif current_path.endswith("pyrewton/cazymes/evaluate_tools/tools"):
        logger.warning(
            (
                "Script invoked in 'tools' directory.\n"
                "Changed current working directory to pyrewton/cazymes/evaluate_tools."
            )
        )
        os.chdir('..')
        return

    else:
        logger.error(
            (
                "Script invokved in wrong directory.\n"
                f"CWD: {current_path}\n"
                "To be able to access dbCAN, CUPP and eCAMI please invokve script when cwd is\n"
                "pyrewton/cazymes/evaluate_tools within the pyrewton programm. Terminating program."
            )
        )
        sys.exit(1)


def get_predictions(args):
    """Build prediction queries and invoke prediction tools for each query.

    :param args: parser object

    Return list of Proteome class objects, queries to prediction tools.
    """
    # create list of paths to all fasta files in input directory
    all_fasta_paths = get_paths.get_file_paths(args.input, suffixes=[".fasta", ".fa"])

    # for each FASTA file invoke dbCAN, CUPP and eCAMI
    for file_path in tqdm(all_fasta_paths, desc="Invoking tools for FASTA file"):  # make tqdm
        # retrieve data on source of protein sequences and species taxonomy ID
        protein_source = get_protein_source(file_path)
        tax_id = get_tax_id(file_path)

        time_stamp = datetime.now().strftime("%Y-%m-%d--%H-%M-%S")

        # name output dir to store prediction output and statistical evaluation
        if tax_id is None:
            outdir_name = f"cazyme_predictions_{protein_source}_{time_stamp}"
        else:
            outdir_name = f"cazyme_predictions_{tax_id}_{protein_source}_{time_stamp}"

        # create output_dir for given input FASTA file within user specified parent directory
        output_path = args.output / outdir_name
        make_output_directory(output_path, args.force, args.nodelete)

        # create Proteome class object to store data on the query made to the prediction tools
        prediction_tool_query = Proteome(file_path, tax_id, protein_source, {"dir": output_path})

        # invoke prediction tools and retrieve paths to the prediction tools outputs
        full_outdir_path = invoke_prediction_tools(prediction_tool_query)

        # copy parsed fasta file into the prediction output dir
        shutil.copy(file_path, full_outdir_path)

    return


def get_protein_source(file_path):
    """Retrieve source of protein sequences from FASTA file path.

    :param file_path: path, path to fasta file

    Return string, source of protein sequences.
    """
    # Attempt to retrieve genomic assembly from FASTA file path
    search_result = re.search(r"GC(A|F)_\d+(_|.)\d+?", str(file_path), re.IGNORECASE)

    try:
        protein_source = search_result.group()
        return protein_source
    except AttributeError:
        # search for uniprot query used to retrieve proteins
        search_result = re.search(r"uniprot(.+?\.)", str(file_path), re.IGNORECASE)
        try:
            protein_source = search_result.group()[:-1]  # remove terminal '.'
            return protein_source
        except AttributeError:
            protein_source = "unknown_source"
            return protein_source


def get_tax_id(file_path):
    """Retrieve taxonomy ID from fasta file path.

    :param file_path: path, path to fasta file path

    Return string, taxonomy ID of proteins' host species.
    """
    # search for first taxonomy ID format
    search_result = re.search(r"txid\d+?\D", str(file_path), re.IGNORECASE)

    try:
        tax_id = search_result.group()[:-1]
        return tax_id
    except AttributeError:      
        # search for other taxonomy ID format
        search_result = re.search(r"taxonomy__\d+?__", str(file_path), re.IGNORECASE)
        try:
            tax_id = search_result.group()[:-2]
            return tax_id
        except AttributeError:
            return


if __name__ == "__main__":
    main()
