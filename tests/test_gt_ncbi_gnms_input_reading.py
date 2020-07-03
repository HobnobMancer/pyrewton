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

"""Tests for pyrewton genbank:get_ncbi_genomes submodule's
read of input files.

These tests are inteded to be run from the root repository using:
pytest -v
"""

import pytest

from pyrewton.genbank.get_ncbi_genomes import get_ncbi_genomes


@pytest.fixture
def test_input_file(gt_ncbi_gnms_input_dir):
    input_reading_path = gt_ncbi_gnms_input_dir / "gt_ncbi_gnms_reading_test_input.txt"
    return input_reading_path


@pytest.mark.run(order=4)
def test_reading_input_file(test_input_file, null_logger, gt_ncbi_gnms_test_inputs):
    """Tests script can open and read supplied input file."""
    get_ncbi_genomes.parse_input_file(
        test_input_file, null_logger, gt_ncbi_gnms_test_inputs[1]
    )
