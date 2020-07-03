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

"""Tests for pyrewton genbank:get_ncbi_genomes submodule, retrieval
of accession number and GenBank file from NCBI.

These tests are inteded to be run from the root repository using:
pytest -v
"""
from argparse import Namespace

import pytest

from Bio import Entrez

from pyrewton.genbank.get_ncbi_genomes import get_ncbi_genomes

Entrez.email = "dumpy.email@my.domain"


@pytest.fixture
def input_df(gt_ncbi_gnms_test_inputs):
    row_data = []
    row_data.append(gt_ncbi_gnms_test_inputs[4])
    row_data.append(gt_ncbi_gnms_test_inputs[5])
    row_data.append(gt_ncbi_gnms_test_inputs[1])
    return row_data


@pytest.fixture
def args():
    argsdict = {"args": Namespace(genbank=False, retries=10, timeout=10)}
    return argsdict


@pytest.mark.run(order=13)
def test_accession_number_retrieval(input_df, null_logger, args):
    """Tests multiplpe Entrez calls to NCBI to retrieve accession numbers."""
    get_ncbi_genomes.get_accession_numbers(input_df, null_logger, args["args"])


@pytest.mark.run(order=14)
def test_compiling_url(null_logger):
    """Test generation of URL for downloading GenBank files."""
    get_ncbi_genomes.compile_url("test_accession", "test_name", null_logger, "suffix")


@pytest.mark.run(order=15)
def test_genbank_download(args, gt_ncbi_gnms_targets, null_logger):
    """Test downloading of GenBank file."""
    get_ncbi_genomes.download_file(
        "http://httpbin.org/get",
        args["args"],
        gt_ncbi_gnms_targets[2],
        null_logger,
        "test_accession",
        "test_file",
    )
