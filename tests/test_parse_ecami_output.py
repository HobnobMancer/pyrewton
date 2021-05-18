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

from pathlib import Path

from pyrewton.cazymes.evaluate_tools.parse import parse_ecami_output


@pytest.fixture
def input_dir(test_dir):
    path_ = test_dir / "test_inputs" / "parsing_predictions_test_inputs" / "ecami"
    return path_


@pytest.fixture
def fasta_path(input_dir):
    path_ = input_dir / "parse_ecami.fasta"
    return path_

# test parse_ecami_output


def test_parse_ecami_output(input_dir, fasta_path):
    """Test parse_ecami_output()."""
    ecami_output_file = input_dir / "parse_ecami_output.txt"

    parse_ecami_output.parse_ecami_output(ecami_output_file, fasta_path)


def test_get_subfamilies_ec_numbers():
    subfam_group = "GH5_1:2|GH18:65|CBM18:48|CBM50:43|CBM66:1|3.2.1.14:1|1.2.3.4:5|ding"

    subfam, ec = parse_ecami_output.get_subfamily_ec_numbers(subfam_group, "GH5")

    assert subfam == "GH5_1"
    assert ec == ["3.2.1.14", "1.2.3.4"]


def test_adding_noncazymes(fasta_path):
    """Test adding non-cazymes from the FASTA file to the prediciton dict."""
    ecami_predictions = {"ATY59188.1 A9K55_002493": "cazyme_prediction"}

    output = parse_ecami_output.add_non_cazymes(ecami_predictions, fasta_path)

    assert list(output.keys()) == ['ATY59188.1 A9K55_002493', 'ATY61648.1 A9K55_009369']
