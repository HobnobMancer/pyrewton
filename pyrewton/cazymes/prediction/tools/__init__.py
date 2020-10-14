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
"""Module invoking CAZyme prediction tools: dbCAN, CUPP and eCAMI"""


def invoke_prediction_tools(query, args, logger):
    """Pass paramaters to CAZyme prediction tools and invoke tool run."""
    fasta_path = query.fasta
    outdir = query.prediction_dir

    # args required by run_dbCAN
    # <fasta_path> protein --output_dir <outdir_path>

    # args required by ecami for prediction.py
    # -input <fasta_path> -kmer_db CAZyme -output <path to output file, including .text ext>

    # args required by cupp prediction.py
    # -query <fasta_path> -output_path <path to output fasta file including .fasta ext>

