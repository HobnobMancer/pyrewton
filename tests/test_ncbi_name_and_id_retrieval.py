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
retrieval of scientific name and taxonomy ID from NCBI.

These tests are inteded to be run from the root repository using:
pytest -v
"""

from Bio import Entrez
import pytest

from pyrewton.genbank.get_ncbi_genomes import get_ncbi_genomes

# Define dummy email for Entrez
Entrez.email = "my.email@my.domain"


@pytest.mark.run(order=11)
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


@pytest.mark.run(order=12)
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
