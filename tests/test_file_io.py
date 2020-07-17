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


@pytest.fixture
def testing_df():
    df_data = [["A", "B", "C"]]
    df = pd.DataFrame(df_data, columns=["C1", "C2", "C3"])
    return df


@pytest.fixture
def testing_df_output_file(test_dir):
    df_output = (
        test_dir / "test_targets" / "file_io_test_targets" / "test_writing_df.csv"
    )
    return df_output


@pytest.fixture
def testing_df_output_dir(test_dir):
    df_output = test_dir / "test_targets" / "file_io_test_targets"
    return df_output


@pytest.fixture
def file_io_args_n_false(testing_df_output_dir):
    argsdict = {
        "args": Namespace(output=testing_df_output_dir, nodelete=False, force=True)
    }
    return argsdict


@pytest.fixture
def file_io_args_n_true(testing_df_output_dir):
    argsdict = {
        "args": Namespace(output=testing_df_output_dir, nodelete=True, force=True)
    }
    return argsdict


@pytest.mark.run(order=6)
def test_output_dir_creation_n_false(file_io_args_n_false, null_logger):
    """Test creation of output dir when args.nodelete is false"""
    file_io.make_output_directory(file_io_args_n_false["args"], null_logger)


@pytest.mark.run(order=7)
def test_output_dir_creation_n_true(file_io_args_n_true, null_logger):
    """Test creation of output dir when args.nodelete is true"""
    file_io.make_output_directory(file_io_args_n_true["args"], null_logger)


@pytest.mark.run(order=8)
def test_writing_df_f_false(testing_df, null_logger, testing_df_output_file):
    """Tests function for writing out created dataframe when force is false"""
    file_io.write_out_dataframe(
        testing_df, null_logger, testing_df_output_file, False, False
    )


@pytest.mark.run(order=9)
def test_writing_df_f_true(testing_df, null_logger, testing_df_output_file):
    """Tests function for writing out created dataframe when force is true"""
    file_io.write_out_dataframe(
        testing_df, null_logger, testing_df_output_file, True, False
    )


@pytest.mark.run(order=10)
def test_writing_named_df_f_false(testing_df, null_logger, testing_df_output_dir):
    """Tests function for writing out a prenamed dataframe"""
    file_io.write_out_pre_named_dataframe(
        testing_df, "unitesting", null_logger, testing_df_output_dir, False, False
    )


@pytest.mark.run(order=11)
def test_writing_named_df_f_true(testing_df, null_logger, testing_df_output_dir):
    """Tests function for writing out a prenamed dataframe"""
    file_io.write_out_pre_named_dataframe(
        testing_df, "unitesting", null_logger, testing_df_output_dir, True, False
    )
