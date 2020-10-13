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

import subprocess
import sys

from dataclasses import dataclass
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
        return f"{self.tax_id} {self.source} fasta: {self.fasta}, pred_dir: {prediction_dir}"


def main():
    """Set up utilities and coordinate prediction of CAZymes."""
    # build parser

    # build logger

    # retrieve paths to fasta files
    fasta_file_paths = get_fasta_paths(args, logger)

    # empty list to store the names of the output directories
    outdirs = []

    # for each FASTA file invoke dbCAN, CUPP and eCAMI
    for file_path in fasta_file_paths:
        # create name of output dir that will house all raw output
        # from prediction tools and parsed outputs
        # possibly "<fasta file name (minus the file extension)>_<tool name>_<time/date stamp>"
        out_dir_name = f""
        outdirs.append(out_dir_name)
        make_output_directory() # ADD PARAM IN!!!!

        # run dbCAN
        subprocess.run()
        # run CUPP
        subprocess.run()
        # run eCAMI
        subprocess.run()

    # parse output from prediction tools to create standardised output format
    # call functions from pyrewton.cazymes.prediction.parse submodule
    for directory in outdirs:
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
    files_in_entries = (entry for entry in Path(args.input).iterdir if entry.is_file())

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
