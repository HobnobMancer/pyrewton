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

"""Tests for pyrewton genbank:get_ncbi_genomes submodule units.

These tests are inteded to be run from the root repository using:
pytest -v
"""

import pytest

from Bio import Entrez

from pyrewton.genbank.get_ncbi_genomes import get_ncbi_genomes

# Define dummy email for Entrez
Entrez.email = "my.email@my.domain"


@pytest.mark.run(order=8)
def test_reading_input_file(
    test_ncbi_species_file, null_logger, gt_ncbi_gnms_test_inputs
):
    """Tests script can open and read supplied input file."""
    get_ncbi_genomes.parse_input_file(
        test_ncbi_species_file, null_logger, gt_ncbi_gnms_test_inputs[1]
    )


# order = 9
@pytest.mark.skip(reason="mocking database call still under development")
def test_scientific_name_retrieval(
    gt_ncbi_gnms_test_inputs, gt_ncbi_gnms_targets, null_logger, monkeypatch
):
    """Tests Entrez call to NCBI to retrieve scientific name from taxonomy ID.
    Tests that correct output is returned from from get_genus_species_name()
    function in Extract_genomes_NCBI.py.
    """

    def mock_entrez_sci_call(*args, **kwargs):
        """Mocks call to Entrez."""
        return [{"ScientificName": "Aspergillus niger"}]

    monkeypatch.setattr(get_ncbi_genomes, "entrez_retry", mock_entrez_sci_call)

    assert gt_ncbi_gnms_targets[0] == get_ncbi_genomes.get_genus_species_name(
        gt_ncbi_gnms_test_inputs[1],
        null_logger,
        gt_ncbi_gnms_test_inputs[3],
        gt_ncbi_gnms_test_inputs[0],
    )


# order 10
@pytest.mark.skip(reason="mocking database call still under development")
def test_taxonomy_id_retrieval(
    gt_ncbi_gnms_test_inputs, gt_ncbi_gnms_targets, null_logger, monkeypatch
):
    """Tests Entrez call to NCBI to retrieve taxonomy ID from scientific name."""

    def mock_entrez_txid_call(*args, **kwargs):
        """Mocks call to Entrez."""
        return {"IdList": ["162425"]}

    monkeypatch.setattr(get_ncbi_genomes, "entrez_retry", mock_entrez_txid_call)

    assert gt_ncbi_gnms_targets[1] == get_ncbi_genomes.get_genus_species_name(
        gt_ncbi_gnms_test_inputs[2],
        null_logger,
        gt_ncbi_gnms_test_inputs[3],
        gt_ncbi_gnms_test_inputs[0],
    )


# order = 11
@pytest.mark.skip(reason="mocking database call still under development")
def test_accession_number_retrieval(input_ncbi_df, null_logger, ncbi_args):
    """Tests multiplpe Entrez calls to NCBI to retrieve accession numbers."""
    get_ncbi_genomes.get_accession_numbers(
        input_ncbi_df, null_logger, ncbi_args["args"]
    )


@pytest.mark.run(order=12)
def test_compiling_url(null_logger):
    """Test generation of URL for downloading GenBank files."""
    get_ncbi_genomes.compile_url("test_accession", "test_name", null_logger, "suffix")


# order = 13
@pytest.mark.skip(reason="mocking database call still under development")
def test_genbank_download(ncbi_args, gt_ncbi_gnms_targets, null_logger):
    """Test downloading of GenBank file."""
    get_ncbi_genomes.download_file(
        "http://httpbin.org/get",
        ncbi_args["args"],
        gt_ncbi_gnms_targets[2],
        null_logger,
        "test_accession",
        "test_file",
    )
