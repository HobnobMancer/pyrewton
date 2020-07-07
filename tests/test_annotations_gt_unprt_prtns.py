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

from pyrewton.annotations import get_uniprot_proteins


@pytest.mark.run(order=15)
def uniprot_df_creation(genbank_df, null_logger, formated_uniprot_results, monkeypatch):
    """Test creation and processing of dataframe of UniProt results."""

    def mock_uniprot_call(*args, **kwargs):
        """Mocks call to UniProtKB database."""
        result = formated_uniprot_results
        return result

    monkeypatch.setattr(get_uniprot_proteins, "get_uniprotkb_data", mock_uniprot_call)

    get_uniprot_proteins.build_uniprot_df(genbank_df, null_logger)


@pytest.mark.run(order=16)
def format_uniprot_results(
    unformated_uniprot_results, uniprot_query_source, null_logger
):
    """Test formating of UniProt search results."""
    get_uniprot_proteins.format_search_results(
        unformated_uniprot_results, uniprot_query_source, null_logger
    )
