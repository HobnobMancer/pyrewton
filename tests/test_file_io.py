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

"""Tests for pyrewton file_io module.

These tests are inteded to be run from the root repository using:
pytest -v
"""

import pytest

from pyrewton import file_io


@pytest.mark.run(order=5)
def test_output_dir_creation(file_io_args, null_logger):
    """Tests function for creating output dir"""
    file_io.make_output_directory(file_io_args["args"], null_logger)


@pytest.mark.run(order=6)
def test_writing_df(testing_df, null_logger, testing_df_output_file):
    """Tests function for writing out created dataframe"""
    file_io.write_out_dataframe(
        testing_df, null_logger, testing_df_output_file, True, False
    )


@pytest.mark.run(order=7)
def test_writing_named_df(testing_df, null_logger, testing_df_output_dir):
    """Tests function for writing out a prenamed dataframe"""
    file_io.write_out_pre_named_dataframe(
        testing_df, "unitesting", null_logger, testing_df_output_dir, True, False
    )
