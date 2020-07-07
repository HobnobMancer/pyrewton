  
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
