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
"""Handle input directory checking, opening and parsing for get_cazyme_annotations.py.

:func

"""

import sys

import pandas as pd

from pathlib import Path


def get_genbank_file(accession, args, logger):
    """Retrieve GenBank file for accession number in local dir.

    :param accession: str, accession number of GenBank file
    :param args: parser arguments
    :param logger: logger object

    Return list of length 1, containing path to GenBank file.
    """
    # replace '.' with '_' to match format in GenBank file name
    file_stem = accession.replace(".", "_")

    # create empty list to store file entries, to allow checking if multiple files were retrieved
    gb_file = []

    # retrieve all files from directory
    files_in_entries = (
        entry for entry in Path(args.genbank).iterdir() if entry.is_file()
    )
    for item in files_in_entries:
        # search for accession number's GenBank file
        if item.name.startswith(f"{file_stem}") and item.name.endswith(".gbff.gz"):
            gb_file.append(item)

    # check file was retrieved, not multiple or none
    if len(gb_file) == 0:
        logger.warning(
            (
                f"Retrieved 0 files for {accession}.\n"
                "Returning null ('NA') value for all protein data"
            )
        )
        return None

    elif len(gb_file) > 1:
        logger.warning(
            (
                f"Retrieved multiple files for {accession}.\n"
                "Returning null ('NA') value for all protein data"
            )
        )
        return None

    # check if files is empty
    if gb_file[0].stat().st_size == 0:
        logger.warning(
            (
                f"GenBank file retrieved for {accession} is empty.\n"
                "Returning null ('NA' value for all protein data"
            )
        )
        return None

    return gb_file
