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
import subprocess
import sys

from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import List, Optional

from pyrewton.file_io import make_output_directory


@dataclass
class SpeciesData:
    """Store species and associated FASTA file data.

    Contains all data for a given species, including path to
    FASTA file containing species proteins, NCBI taxonomy ID and
    the source of the protein sequences, i.e the accession number
    of the genomic assembly or the remote database, e.g. UniProt.
    """

    tax_id: str  # NCBI taxonomy id, prefix NCBI:txid
    fasta: Path  # path to FASTA file containing species protein
    source: str  # source of protein sequences, genomic assembly or database
    prediction_dir: Path  # path to output dir of prediciton tools

    def __str__(self):
        """Representation of object"""
        return f"{self.tax_id} {self.source} fasta: {self.fasta}, pred_dir: {self.prediction_dir}"


def main(argv: Optional[List[str]] = None, logger: Optional[logging.Logger] = None):
    """Set up parser and logger, and coordinate prediction of CAZymes."""
    # build parser

    # build logger

    # create list of paths to all fasta files in input directory
    all_fasta_paths = get_fasta_paths(args, logger)

    # for each FASTA file invoke dbCAN, CUPP and eCAMI
    for file_path in all_fasta_paths:
        # retrieve taxonomy ID from file path
        tax_pattern = re.compile(r"ncbi-txid\d+?", re.IGNORECASE)
        match_result = re.match(tax_pattern, str(file_path), re.IGNORECASE)
        tax_id = match_result.group()

        # retrieve protein source name, e.g. genomic assembly or UniProt
        # try to find accession number of source genomic assembly
        accession_pattern = re.compile(r"GC(A|F)\d+?_\d+?", re.IGNORECASE)
        match_result = re.match(accession_pattern, str(file_path), re.IGNORECASE)
        if match_result.group() is None:
            protein_source = "UniProt"
        else:
            protein_source = match_result.group()

        # name output dir to store prediction output and statistical evaluation
        time_stamp = datetime.now().strftime("%Y-%m-%d--%H-%M-%S")
        outdir_name = f"cazyme_predictions_{tax_id}_{protein_source}_{time_stamp}"

        # create output_dir
        make_output_directory(outdir_name, logger, args.force, args.nodelete)
        outdir_path = Path(outdir_name)

        # create class object to store species, FASTA file and output data
        class_data = SpeciesData(tax_id, file_path, protein_source, outdir_path)

        # pass path to FASTA file to prediction tools and invoke CAZyme prediction
        invoke_prediction_tools()  # what inputs do each need

        # possibly "<fasta file name (minus the file extension)>_<tool name>_<time/date stamp>"
        out_dir_name = f""
        outdirs.append(out_dir_name)
        make_output_directory()  # ADD PARAM IN!!!!

        # run dbCAN
        subprocess.run()
        # run CUPP
        subprocess.run()
        # run eCAMI
        subprocess.run()

    # parse output from prediction tools to create standardised output format
    # call functions from pyrewton.cazymes.prediction.parse submodule
    # for directory in outdirs:
    # parse dbCAN output

    # parse CUPP output

    # parse eCAMI output


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


if __name__ == "__main__":
    main()
