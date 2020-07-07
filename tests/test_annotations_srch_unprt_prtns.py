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

"""Tests for pyrewton.annotations search_uniprot_proteins.py script.

These tests are intended to be run from the root repository using:
pytest -v
"""

import pytest

from pyrewton.annotations import search_uniprot_proteins


@pytest.mark.run(order=17)
def test_uniprot_cazyme_filter(formated_uniprot_results, null_logger):
    """Test cazyme filstering of UniProt search results."""
    search_uniprot_proteins.get_cazyme_subset_df(formated_uniprot_results, null_logger)


@pytest.mark.run(order=18)
def test_ec_go_cazyme_comparison(formated_uniprot_results, null_logger):
    """Test cazyme filstering of UniProt search results."""
    search_uniprot_proteins.compare_cazyme_dfs(formated_uniprot_results, null_logger)
