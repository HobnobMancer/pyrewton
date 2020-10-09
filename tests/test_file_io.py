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

from argparse import Namespace

import pytest

import pandas as pd

from pyrewton import file_io


# Create general fixtures for tests in this file


@pytest.fixture
def testing_output_dir(test_dir):
    test_outputs = test_dir / "test_targets" / "file_io_test_targets"
    return test_outputs


# Create fixtures for testing make_output_directory


@pytest.fixture
def making_output_dir(testing_output_dir):
    output_dir = testing_output_dir / "test_output_dir"
    return output_dir


@pytest.fixture
def args_make_dir_nd_true(making_output_dir):
    return [making_output_dir, True, True]


@pytest.fixture
def args_make_dir_nd_false(making_output_dir):
    return [making_output_dir, False, True]


# Create fixtures for testing writing out a dataframe


@pytest.fixture
def testing_df():
    df_data = [["A", "B", "C"]]
    df = pd.DataFrame(df_data, columns=["C1", "C2", "C3"])
    return df


@pytest.fixture
def df_output_file(test_dir):
    df_output = (
        test_dir / "test_targets" / "file_io_test_targets" / "test_writing_df.csv"
    )
    return df_output


# Test make_output_directory


@pytest.mark.run(order=6)
def test_output_dir_creation_nd_true(args_make_dir_nd_true, null_logger):
    """Test creation of output dir when args.nodelete is false"""
    file_io.make_output_directory(
        args_make_dir_nd_true[0],
        null_logger,
        args_make_dir_nd_true[1],
        args_make_dir_nd_true[2],
    )


@pytest.mark.run(order=7)
def test_output_dir_creation_n_false(args_make_dir_nd_false, null_logger):
    """Test creation of output dir when args.nodelete is true"""

    file_io.make_output_directory(
        args_make_dir_nd_false[0],
        null_logger,
        args_make_dir_nd_false[1],
        args_make_dir_nd_false[2],
    )


# Test write_out_dataframe


@pytest.mark.run(order=8)
def test_writing_df_f_true(testing_df, null_logger, df_output_file):
    """Tests function for writing out created dataframe when force is true"""
    file_io.write_out_dataframe(testing_df, null_logger, df_output_file, True)


@pytest.mark.run(order=9)
def test_writing_df_f_false(testing_df, null_logger, df_output_file):
    """Tests function for writing out created dataframe when force is false"""
    file_io.write_out_dataframe(testing_df, null_logger, df_output_file, False)


# Test write_out_pre_named_dataframe


@pytest.mark.run(order=10)
def test_writing_named_df_f_true(testing_df, null_logger, making_output_dir):
    """Tests function for writing out a prenamed dataframe"""
    file_io.write_out_pre_named_dataframe(
        testing_df, "test_writing_df.csv", null_logger, making_output_dir, True
    )


@pytest.mark.run(order=11)
def test_writing_named_df_f_false(testing_df, null_logger, making_output_dir):
    """Tests function for writing out a prenamed dataframe"""
    file_io.write_out_pre_named_dataframe(
        testing_df, "test_writing_df.csv", null_logger, making_output_dir, False
    )
