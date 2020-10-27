#!/usr/bin/env python
# -*- coding: utf-8 -*-
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
"""Predict CAZymes from protein sequences and evaluate prediction accuracy.

:cmd_args   :

:func    :

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

from pyrewton.cazymes.prediction.tools import invoke_prediction_tools
from pyrewton.file_io import make_output_directory
from pyrewton.loggers import build_logger
from pyrewton.parsers.parser_predict_cazymes import build_parser


@dataclass
class Query:
    """Data of prediction tool query.

    :param fasta: path, path to query fasta file
    :param tax_id: str, NCBI taxonomy ID of host species
    :param source: str, source of proteins within fasta file
    :param prediction_dir: path, path to which prediction tool output is written
    """

    fasta: Path  # path to FASTA file containing proteins for prediction
    tax_id: str  # NCBI taxonomy id, prefix NCBI:txid
    source: str  # source of protein sequences, genomic assembly or database
    prediction_dir: Path  # path to dir containing outputs from prediciton tools

    def __str__(self):
        """Representation of object"""
        return f"{self.tax_id} {self.source} fasta: {self.fasta}, pred_dir: {self.prediction_dir}"


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
        args = build_parser(argv).parse_args()

    # Initiate logger
    # Note: log file only created if specified at cmdline
    if logger is None:
        logger = build_logger("predict_cazymes", args)

    # create list of paths to all fasta files in input directory
    all_fasta_paths = get_fasta_paths(args, logger)

    # create empty list to store all instances of Prediction class objects
    predictions = []

    # for each FASTA file invoke dbCAN, CUPP and eCAMI
    for file_path in all_fasta_paths:  # make tqdm
        # retrieve data on source of protein sequences and species taxonomy ID
        protein_source = get_protein_source(file_path, args, logger)
        tax_id = get_tax_id(file_path, args, logger)

        time_stamp = datetime.now().strftime("%Y-%m-%d--%H-%M-%S")

        # name output dir to store prediction output and statistical evaluation
        if tax_id is None:
            outdir_name = f"cazyme_predictions_{protein_source}_{time_stamp}"
        else:
            outdir_name = f"cazyme_predictions_{tax_id}_{protein_source}_{time_stamp}"

        # create output_dir
        make_output_directory(outdir_name, logger, args.force, args.nodelete)
        outdir_path = Path(outdir_name)

        # create FastaFile class object to store data for fasta file
        prediction_tool_query = Query(file_path, tax_id, protein_source, outdir_path)

        # pass FASTA file path and outdir_path to invoke prediction tools
        invoke_prediction_tools(prediction_tool_query, args, logger)

    # standardist output from prediction tools for the prediction output per
    # FASTA file
    # for prediction in predictions:  # make tqdm
    #

    # statistical evaluate CAZyme predictions per FASTA file input
    # for prediction in predictions:  # make tqdm

    # calculate overall statistcal evaluation of prediction

    # write reports of statistical evaluation
    # for prediction in predictions:  # make tqdm


def get_fasta_paths(args, logger):
    """Retrieve paths to call FASTA files in input dir.

    :param args: parser object
    :param logger: logger object

    Returns list of paths to fasta files.
    """
    # create empty list to store the file entries, to allow checking if no files returned
    fasta_file_paths = []

    # retrieve all files from input directory
    files_in_entries = (
        entry for entry in Path(args.input).iterdir() if entry.is_file()
    )
    # retrieve paths to fasta files
    for item in files_in_entries:
        # search for fasta file extension
        if item.name.endswith(".fasta") or item.name.endswith(".fa"):
            fasta_file_paths.append(item)

    # check files were retrieved from input directory
    if len(fasta_file_paths) == 0:
        logger.warning(
            (
                "No FASTA files retrieved from input directory.\n"
                "Check the path to the input directory is correct.\n"
                "Terminanting program."
            )
        )
        sys.exit(1)

    return fasta_file_paths


def get_protein_source(file_path, args, logger):
    """Retrieve source of protein sequences from FASTA file path.

    :param file_path: path, path to fasta file
    :param args: parser objects
    :param logger: logger object

    Return string, source of protein sequences.
    """
    # Attempt to retrieve genomic assembly from FASTA file path
    accession_pattern = re.compile(r"GC(A|F)\d+?_\d+?")
    search_result = re.search(accession_pattern, str(file_path), re.IGNORECASE)

    try:
        protein_source = search_result.group()
        return protein_source
    except AttributeError:
        # search for uniprot query used to retrieve proteins
        query_pattern = re.compile(r"uniprot(.+?\.)")
        search_result = re.search(query_pattern, str(file_path), re.IGNORECASE)
        try:
            protein_source = search_result.group()[:-1]  # remove terminal '.'
            return protein_source
        except AttributeError:
            protein_source = "unknown_source"
            return protein_source


def get_tax_id(file_path, args, logger):
    """Retrieve taxonomy ID from fasta file path.

    :param file_path: path, path to fasta file path
    :param args: parser object
    :param logger: logger object

    Return string, taxonomy ID of proteins' host species.
    """
    # search for first taxonomy ID format
    tax_pattern = re.compile(r"ncbi(-|_)txid\d+?")
    search_result = re.search(tax_pattern, str(file_path), re.IGNORECASE)

    try:
        tax_id = search_result.group()
        return tax_id
    except AttributeError:
        # search for other taxonomy ID format
        tax_pattern = re.compile(r"taxonomy__\d+?__")
        search_result = re.search(tax_pattern, str(file_path), re.IGNORECASE)
        try:
            tax_id = search_result.group()
            return tax_id
        except AttributeError:
            return


if __name__ == "__main__":
    main()
