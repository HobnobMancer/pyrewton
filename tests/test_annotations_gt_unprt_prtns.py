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

"""Tests for pyrewton.annotations get_uniprot_proteins.py script.

These tests are inteded to be run from the root repository using:
pytest -v
"""

import pytest

from argparse import Namespace

from bioservices import UniProt

from pyrewton.annotations import get_uniprot_proteins


@pytest.fixture
def genbank_df_path(test_dir):
    genbank_df = (
        test_dir / "test_inputs" / "gt_gnbnk_anntns_test_inputs" / "test_input_df.csv"
    )
    return genbank_df


@pytest.fixture
def genbank_files_dir(test_dir):
    genbank_dir = (
        test_dir / "test_inputs" / "gt_gnbnk_anntns_test_inputs" / "test_genbank_files"
    )
    return genbank_dir


@pytest.fixture
def gnbnk_anno_args(genbank_files_dir):
    argsdict = {"args": Namespace(genbank=genbank_files_dir)}
    return argsdict


@pytest.fixture
def genbank_df(genbank_df_path):
    input_df = pd.read_csv(genbank_df_path, header=0, index_col=0)
    return input_df


@pytest.fixture
def formated_uniprot_results(test_dir):
    df_dir = (
        test_dir
        / "test_inputs"
        / "gt_unprt_prtns_test_inputs"
        / "mocked_formated_uniprot_results.csv"
    )
    df = pd.read_csv(df_dir, header=0, index_col=0)
    return df


@pytest.fixture
def unformated_uniprot_results(test_dir):
    uniprot_results = (
        test_dir
        / "test_inputs"
        / "gt_unprt_prtns_test_inputs"
        / "mocked_unformated_uniprot_results.txt"
    )
    with open(uniprot_results, "r") as results:
        data = results
    return data


@pytest.fixture
def uniprot_query_source(formated_uniprot_results):
    """Return single row, a pandas series."""
    return formated_uniprot_results.head(1)


@pytest.mark.run(order=15)
def test_uniprot_df_creation(
    genbank_df, null_logger, unformated_uniprot_results, monkeypatch
):
    """Test creation and processing of dataframe of UniProt results."""

    def mock_bioservices_uniprot(*args, **kwargs):
        """Mocks call to UniProtKB database."""
        return unformated_uniprot_results

    monkeypatch.setattr(UniProt, "search", mock_bioservices_uniprot)

    get_uniprot_proteins.build_uniprot_df(genbank_df, null_logger)


@pytest.mark.run(order=16)
def test_format_uniprot_results(
    unformated_uniprot_results, uniprot_query_source, null_logger
):
    """Test formating of UniProt search results."""
    get_uniprot_proteins.format_search_results(
        unformated_uniprot_results, uniprot_query_source, null_logger
    )
