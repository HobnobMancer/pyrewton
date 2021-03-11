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

"""Tests for parsing output from eCAMI.

These tests are inteded to be run from the root repository using:
pytest -v
"""

import pytest

import numpy as np

from pathlib import Path

from pyrewton.cazymes.prediction.parse import parse_cupp_output


@pytest.fixture
def input_dir(test_dir):
    path_ = test_dir / "test_inputs" / "parsing_predictions_test_inputs" / "cupp"
    return path_


@pytest.fixture
def fasta_path(input_dir):
    path_ = input_dir / "parse_cupp.fasta"
    return path_


@pytest.fixture
def prediction_text():
    line = 'ATY67037.1 A9K55_000868	GT41	GT41:16.1	56.5	59.9	1355	1234..1496	213'
    prediction_text = line.split("\t")
    return prediction_text


def test_parse_cupp_output(input_dir, fasta_path):
    """Test the function that handles the overall parsing of the CUPP output."""
    cupp_log_file = input_dir / "parse_cupp_output.log"

    parse_cupp_output.parse_cupp_output(cupp_log_file, fasta_path)


def test_cupp_domain_range_retrieval(prediction_text):
    """Test retrieving domain ranges from CUPP output."""
    output = parse_cupp_output.get_cupp_domain_range(prediction_text)

    assert ['1234..1496'] == output


def test_cupp_ec_number_retrieval(prediction_text):
    """Test retrieval of EC numbers from CIPP output."""
    output = parse_cupp_output.get_cupp_ec_number(prediction_text)

    assert np.isnan(output)


def test_adding_noncazyme_to_cupp(fasta_path):
    """Test adding non-CAZymes to CUPP predictions."""
    cupp_predictions = {"ATY67037.1 A9K55_000868": "CAZymes"}

    parse_cupp_output.add_non_cazymes(fasta_path, cupp_predictions)
