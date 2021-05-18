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

from pyrewton.cazymes.evaluate_tools.parse import(
    parse_dbcan_output,
    CazymeProteinPrediction,
    CazymeDomain,
)


@pytest.fixture
def input_dir(test_dir):
    path_ = test_dir / "test_inputs" / "parsing_predictions_test_inputs" / "dbcan"
    return path_


@pytest.fixture
def fasta_path(input_dir):
    path_ = input_dir / "parse_dbcan.fasta"
    return path_


@pytest.fixture
def hotpep_file(input_dir):
    path_ = input_dir / "Hotpep.out"
    return path_


@pytest.fixture
def hotpep_file_short(input_dir):
    path_ = input_dir / "Hotpep_short.out"
    return path_


@pytest.fixture
def overview_file(input_dir):
    path_ = input_dir / "overview.txt"
    return path_


@pytest.fixture
def mock_cazymes():
    cazymes = []

    cazyme = CazymeProteinPrediction("test", "ATY58396.1", 1)
    domain = CazymeDomain("test", "ATY58396.1", "AA6", np.nan, ["1.6.5.6"])
    cazyme.cazyme_domains.append(domain)
    cazymes.append(cazyme)

    for acc in ["ATY58475.1", "ATY58478.1", "ATY58486.1"]:
        cazyme = CazymeProteinPrediction("test", acc, 1)
        domain = CazymeDomain("test", acc, "GH3", "GH3_1")
        cazyme.cazyme_domains.append(domain)
        cazymes.append(cazyme)

    return cazymes


def test_parse_all_output(fasta_path, hotpep_file, input_dir, overview_file):
    """Test all functions for parsing dbCAN output."""
    hmmer, hotpep, diamond, dbcan = parse_dbcan_output.parse_dbcan_output(
        overview_file,
        hotpep_file,
        fasta_path,
    )

    assert len(hmmer) == 9287
    assert len(hotpep) == 9287
    assert len(diamond) == 9287
    assert len(dbcan) == 9287


def test_parsing_hmmer():
    """Test parsing HMMER data from the overview.txt file."""
    data = "GT2_Glyco_tranf_2_3(265-498)+CBM66(1301-1461)"

    parse_dbcan_output.get_hmmer_prediction(data, "test_accession")


def test_getting_hmmer_domain_predictions():
    """Test getting HMMER domain predictions from the overview.txt data."""
    data = ["CBM66(1301-1461)", "GT2_Glyco_tranf_2_3(265-498)"]

    parse_dbcan_output.get_hmmer_cazyme_domains(data, "test_accession")


def test_getting_cazyfam_from_hmmer():
    """Test retrieving the CAZy family from the HMMER data."""
    data = "CBM66"

    o1, o2 = parse_dbcan_output.get_hmmer_cazy_family(data, "test_accession")

    assert o1 == "CBM66"
    assert np.isnan(o2)


def test_parsing_all_hotpep_data():
    """Test coordinating all Hotpep data parsing functions."""
    data = "GH18(92)+CBM1(13)+CBM19(1)"

    parse_dbcan_output.get_hotpep_prediction(data, "test_accession")


def test_getting_hotpep_domains():
    """Test building CazymeDomain instances from Hotpep data."""
    data = "GH18(92)+CBM1(13)+CBM19(1)".split("+")

    parse_dbcan_output.get_hotpep_cazyme_domains(data, "test_accession")


def test_adding_ecs_to_hotpep(hotpep_file_short, mock_cazymes):
    """Test adding Hotpep EC number predictions to CazymeDomain instances."""
    overview_dict = {
        "ATY58396.1": {"Hotpep": mock_cazymes[0]},
        "ATY58475.1": {"Hotpep": mock_cazymes[1]},
        "ATY58478.1": {"Hotpep": mock_cazymes[2]},
        "ATY58486.1": {"Hotpep": mock_cazymes[3]},
    }
    parse_dbcan_output.get_hotpep_ec_numbers(overview_dict, hotpep_file_short)


def test_get_diamond_prediction():
    """Test get_diamond_prediction(), which coordinates parsing all DIAMOND data."""
    data = "CBM1+GH18"

    parse_dbcan_output.get_diamond_prediction(data, "test_accession")


def test_get_diamond_predicted_domains():
    """Test creating CazymeDomain instances from DIAMOND data."""
    data = "CBM1+GH18".split("+")

    parse_dbcan_output.get_diamond_predicted_domains(data, "test_accession")
