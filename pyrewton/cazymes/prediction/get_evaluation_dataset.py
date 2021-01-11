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
"""Retrieves datasets for evaluating the CAZyme prediction tools: dbCAN, CUPP and eCAMI.

:cmd_args   :

:func    :

Writes out a FASTA per candidate species, containing all the protein sequences to analysed by the
CAZymes prediction tools. The datasets contain an equal number of CAZymes to non-CAZymes.
"""

import gzip
import logging
import re
import sys

from pathlib import Path
from typing import List, Optional

import pandas as pd

from pyrewton.utilities import build_logger
from pyrewton.utilities.cmd_get_evaluation_dataset import build_parser
from pyrewton.utilities.file_io import make_output_directory


def main():
    """Coordinate preparation for script, and terminating the programme when finished."""
    return


def get_the_fasta_file_paths():
    """Retrieve the paths to all the FASTA files in the input directory."""
    return


def build_protein_dataframe():
    """Build a dataframe containing the protein data for the current working input FASTA file."""
    return


def get_cazy_classification():
    """For each protein check if classified as a CAZyme by CAZy, and its annotated CAZy families."""
    return


def get_dataset():
    """Retrieve dataset of equal number of CAZymes and non-CAZymes."""
    return


def write_out_dataset():
    """Write out the dataset of CAZymes and non-CAZymes to a single FASTA file."""
    return
