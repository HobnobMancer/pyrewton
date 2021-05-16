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

"""Tests for build for the prediction.parse.__init__.py file pyrewton.

These tests are inteded to be run from the root repository using:
pytest -v
"""

import pytest

from pyrewton.cazymes.evaluate_tools import parse


# test building instances of CazymeDomain


def test_cazyme_domain_no_optionals():
    """Test building CazymeDomain when no of the options are passed."""

    domain = parse.CazymeDomain(
        prediction_tool="test_tool",
        protein_accession="test_accession",
        cazy_family="test_family",
    )

    assert "-CazymeDomain in test_accession, fam=test_family, subfam=nan-" == str(domain)
    domain


def test_cazyme_domain_with_optionals():
    """Test building CazymeDomain when options are passed."""

    domain = parse.CazymeDomain(
        prediction_tool="test_tool",
        protein_accession="test_accession",
        cazy_family="test_family",
        cazy_subfamily="test_subfamily",
        ec_numbers=["EC1", "EC2"],
        domain_range=["DR1"],
    )

    assert "-CazymeDomain in test_accession, fam=test_family, subfam=test_subfamily-" == str(domain)
    domain
