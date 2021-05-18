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

from pyrewton.utilities import file_io


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


def test_output_dir_creation_nd_true(making_output_dir):
    """Test creation of output dir when nodelete is false"""
    file_io.make_output_directory(making_output_dir, True, False)


def test_output_dir_creation_nd_false(making_output_dir):
    """Test creation of output dir when nodelete is true"""
    file_io.make_output_directory(making_output_dir, True, True)


def test_make_existing_dir(testing_output_dir):
    """Test creation of output dir when it already exists."""
    path_ = testing_output_dir / "test_dir"
    # run twice to ensure directory exists
    file_io.make_output_directory(path_, True, True)
    file_io.make_output_directory(path_, True, True)


# Test write_out_dataframe


def test_writing_df_f_true(testing_df, df_output_file):
    """Tests function for writing out created dataframe when force is true"""
    file_io.write_out_dataframe(testing_df, df_output_file, True)


def test_writing_df_f_false(testing_df,  df_output_file):
    """Tests function for writing out created dataframe when force is false"""
    file_io.write_out_dataframe(testing_df, df_output_file, False)


def test_writing_df_no_df(df_output_file):
    """Tests function for writing out created dataframe when no df is given"""
    file_io.write_out_dataframe(None, df_output_file, False)


# Test write_out_pre_named_dataframe


def test_writing_named_df_f_true(testing_df, making_output_dir):
    """Tests function for writing out a prenamed dataframe"""
    file_io.write_out_pre_named_dataframe(
        testing_df, "test_writing_df.csv", making_output_dir, True
    )


def test_writing_named_df_f_false(testing_df, making_output_dir):
    """Tests function for writing out a prenamed dataframe"""
    file_io.write_out_pre_named_dataframe(
        testing_df, "test_writing_df.csv", making_output_dir, False
    )

def test_writing_named_df_no_df(df_output_file):
    """Tests function for writing out created dataframe when no df is given"""
    file_io.write_out_pre_named_dataframe(None, "df_name.csv", df_output_file, False)
