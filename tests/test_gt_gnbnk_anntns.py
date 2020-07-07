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

"""Tests for pyrewton genbank:get_genbank_annotations submodule units.

These tests are inteded to be run from the root repository using:
pytest -v
"""

import pytest

from pyrewton.genbank.get_genbank_annotations import get_genbank_annotations


@pytest.mark.run(order=14)
def test_genbank_annotation_df_creation(genbank_df, gnbnk_anno_args, null_logger):
    """Test the function to build a dataframe of all protein annotations
    from GenBank files."""
    get_genbank_annotations.create_dataframe(
        genbank_df, gnbnk_anno_args["args"], null_logger
    )
