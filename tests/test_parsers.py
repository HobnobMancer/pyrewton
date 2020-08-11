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

"""Tests for building of parsers for pyrewton.

These tests are inteded to be run from the root repository using:
pytest -v
"""

import pytest

from pyrewton.parsers import (
    parser_get_ncbi_genomes,
    parser_get_genbank_annotations,
    parser_get_uniprot_proteins,
)


@pytest.mark.run(order=1)
def test_parser_gt_ncb_gnms():
    """Tests building of parser"""
    parser_get_ncbi_genomes.build_parser()


@pytest.mark.run(order=2)
def test_parser_gt_czym_anno():
    """Tests building of parser"""
    parser_get_genbank_annotations.build_parser()


@pytest.mark.run(order=3)
def test_parser_gt_unprt_prtns():
    """Tests building of parser"""
    parser_get_uniprot_proteins.build_parser()
