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

import json
import pytest

from argparse import Namespace

from Bio import Entrez

from pyrewton.genbank.get_ncbi_genomes import get_ncbi_genomes

# Define dummy email for Entrez
Entrez.email = "my.email@my.domain"


@pytest.fixture
def gt_ncbi_gnms_input_dir(test_dir):
    input_dir = test_dir / "test_inputs" / "gt_ncbi_gnms_test_inputs"
    return input_dir


@pytest.fixture
def gt_ncbi_gnms_test_inputs(gt_ncbi_gnms_input_dir):
    with (gt_ncbi_gnms_input_dir / "gt_ncbi_gnms_test_inputs.json").open("r") as ifh:
        test_inputs = json.load(ifh)
        inputs = []
        inputs.append(test_inputs["retries"])
        inputs.append(test_inputs["taxonomy_id"])
        inputs.append(test_inputs["genus_species_name"])
        inputs.append(test_inputs["line_number"])
        inputs.append(test_inputs["df_row_genus"])
        inputs.append(test_inputs["df_row_species"])
    return inputs


@pytest.fixture
def test_ncbi_species_file(gt_ncbi_gnms_input_dir):
    input_reading_path = gt_ncbi_gnms_input_dir / "gt_ncbi_gnms_reading_test_input.txt"
    return input_reading_path


@pytest.fixture
def tax_id_line():
    """Example line from input file."""
    return "NCBI:txid5061"


@pytest.fixture
def species_line():
    """Example line from onput file containing scientific name."""
    return "Aspergillus niger"


# define entrez mocking here!!!!


@pytest.fixture
def input_ncbi_df(gt_ncbi_gnms_test_inputs):
    row_data = []
    row_data.append(gt_ncbi_gnms_test_inputs[4])
    row_data.append(gt_ncbi_gnms_test_inputs[5])
    row_data.append(gt_ncbi_gnms_test_inputs[1])
    return row_data


@pytest.fixture
def na_df_row():
    """Create list to present pandas series containing only 'NA'."""
    mock_df_row = ["NA", "NA", "NA"]
    return mock_df_row


@pytest.fixture
def ncbi_args():
    argsdict = {"args": Namespace(genbank=False, retries=10, timeout=10)}
    return argsdict


@pytest.fixture
def gt_ncbi_gnms_targets(test_dir):
    target_dir = test_dir / "test_targets" / "gt_ncbi_gnms_test_targets"
    with (target_dir / "gt_ncbi_gnms_test_targets.json").open("r") as ifht:
        test_targets = json.load(ifht)
        target_sci_name = test_targets["target_genus_species_name"]
        target_tax_id = test_targets["target_taxonomy_id"]
    return target_sci_name, target_tax_id, target_dir


@pytest.mark.run(order=8)
def test_reading_input_file(
    test_ncbi_species_file, null_logger, gt_ncbi_gnms_test_inputs
):
    """Tests script can open and read supplied input file."""
    get_ncbi_genomes.parse_input_file(
        test_ncbi_species_file, null_logger, gt_ncbi_gnms_test_inputs[1]
    )


@pytest.mark.run(order=9)
def test_input_file_check(test_dir, null_logger, gt_ncbi_gnms_test_inputs):
    """Test get_ncbi_genomes ability to detect when no input file available."""
    get_ncbi_genomes.parse_input_file(test_dir, null_logger, gt_ncbi_gnms_test_inputs)


@pytest.mark.run(order=10)
def test_line_parsing(tax_id_line, species_line, null_logger, monkeypatch):
    """Test function parse_line ability to parse lines of input file."""

    def mock_ncbi_data_retrieval(*args, **kwargs):
        """Mock results from functions with retrieve data from ncbi."""
        return "mock line"

    monkeypatch.setattr(
        get_ncbi_genomes, "get_genus_species_name", mock_ncbi_data_retrieval
    )
    monkeypatch.setattr(get_ncbi_genomes, "get_tax_id", mock_ncbi_data_retrieval)

    get_ncbi_genomes.parse_line(tax_id_line, null_logger, 1, 1)
    get_ncbi_genomes.parse_line(species_line, null_logger, 1, 1)


@pytest.mark.run(order=11)
def test_scientific_name_retrieval(
    gt_ncbi_gnms_test_inputs, gt_ncbi_gnms_targets, null_logger, monkeypatch
):
    """Tests Entrez call to NCBI to retrieve scientific name from taxonomy ID.
    Tests that correct output is returned from from get_genus_species_name()
    function in Extract_genomes_NCBI.py.
    """

    def mock_entrez_sci_call(*args, **kwargs):
        """Mocks call to Entrez to retrieve scientific name."""
        return [{"ScientificName": "Aspergillus niger"}]

    monkeypatch.setattr(get_ncbi_genomes, "Entrez.efetch", mock_entrez_sci_call)

    assert gt_ncbi_gnms_targets[0] == get_ncbi_genomes.get_genus_species_name(
        gt_ncbi_gnms_test_inputs[1],
        null_logger,
        gt_ncbi_gnms_test_inputs[3],
        gt_ncbi_gnms_test_inputs[0],
    )


@pytest.mark.run(order=12)
def test_scientific_name_retrieval_indexerror_catch(
    gt_ncbi_gnms_test_inputs, null_logger, monkeypatch
):
    """Tests get_scientific name retrieval handling indexError catching."""
        def mock_entrez_sci_call(*args, **kwargs):
        """Mocks call to Entrez to retrieve scientific name."""
        empty_list = []
        return empty_list

    monkeypatch.setattr(get_ncbi_genomes, "Entrez.efetch", mock_entrez_sci_call)

    get_ncbi_genomes.get_genus_species_name(
        gt_ncbi_gnms_test_inputs[1],
        null_logger,
        gt_ncbi_gnms_test_inputs[3],
        gt_ncbi_gnms_test_inputs[0],
    )


@pytest.mark.run(order=13)
def test_taxonomy_id_retrieval(
    gt_ncbi_gnms_test_inputs, gt_ncbi_gnms_targets, null_logger, monkeypatch
):
    """Tests Entrez call to NCBI to retrieve taxonomy ID from scientific name."""

    def mock_entrez_txid_call(*args, **kwargs):
        """Mocks call to Entrez to retrieve taxonomy ID."""
        return {"IdList": ["162425"]}

    monkeypatch.setattr(get_ncbi_genomes, "Entrez.esearch", mock_entrez_txid_call)

    assert gt_ncbi_gnms_targets[1] == get_ncbi_genomes.get_genus_species_name(
        gt_ncbi_gnms_test_inputs[2],
        null_logger,
        gt_ncbi_gnms_test_inputs[3],
        gt_ncbi_gnms_test_inputs[0],
    )


@pytest.mark.run(order=14)
def test_tax_id_check(null_logger):
    """Tests searching of input scientific name for digits"""
    get_ncbi_genomes.get_tax_id("5061", null_logger, 1, 1)


@pytest.mark.run(order=15)
def test_tax_id_retrieval_indexerror_catch(
    gt_ncbi_gnms_test_inputs, gt_ncbi_gnms_targets, null_logger, monkeypatch
):
    """Tests handling index Error when retrieving tax ID"""

    def mock_entrez_txid_call(*args, **kwargs):
        """Mocks call to Entrez to retrieve taxonomy ID."""
        return {"IdList": []}

    monkeypatch.setattr(get_ncbi_genomes, "Entrez.esearch", mock_entrez_txid_call)

    get_ncbi_genomes.get_genus_species_name(
        gt_ncbi_gnms_test_inputs[2],
        null_logger,
        gt_ncbi_gnms_test_inputs[3],
        gt_ncbi_gnms_test_inputs[0],
    )


@pytest.mark.run(order=16)
def test_df_cell_content_check(na_df_row, null_logger):
    """Test get_accession_numbers checking of dataframe cell content."""
    get_ncbi_genomes.get_accession_numbers(na_df_row, null_logger, "args")









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
