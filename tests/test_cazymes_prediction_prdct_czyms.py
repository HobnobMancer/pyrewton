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

"""Tests for pyrewton cazymes:prediction:predict_cazymes.

These tests are inteded to be run from the root repository using:
pytest -v
"""

import pytest

import pandas as pd

from argparse import Namespace, ArgumentParser
from collections import namedtuple

from pyrewton.cazymes.prediction import predict_cazymes


# create pytest fixtures for tests

@pytest.fixture
def input_dir(test_dir):
    path = test_dir / "test_inputs" / "prdct_czyms_test_inputs"
    return path

@pytest.fixture
def output_dir(test_dir):
    path = test_dir / "test_targets" / "prdct_czyms_test_targets"
    return path


# test function 'main'


def test_main(output_dir, null_logger, input_dir, monkeypatch):
    """Test function 'main'."""

    def mock_built_parser(*args, **kwargs):
        parser_args = ArgumentParser(
            prog="predict_cazymes.py",
            usage=None,
            description="Predict cazymes from protein sequences",
            conflict_handler="error",
            add_help=True,
        )
        return parser_args

    def mock_parser(*args, **kwargs):
        parser = Namespace(
            force=True, output=output_dir, nodelete=True, input=input_dir
        )
        return parser

    def mock_build_logger(*args, **kwargs):
        return null_logger

    def mock_cwd_check(*args, **kwargs):
        return

    def mock_making_dir(*args, **kwargs):
        return

    def mock_fasta_retrieval(*args, **kwargs):
        return ["fasta_1", "fasta_2", "fasta_3"]

    def mock_protein_source_retrieval(*args, **kwargs):
        return "protein_souce"

    def mock_tax_id_retrieval(*args, **kwargs):
        return "tax_id"

    def mock_prediction_tool_invoking(*args, **kwargs):
        return

    monkeypatch.setattr(predict_cazymes, "build_parser", mock_built_parser)
    monkeypatch.setattr(ArgumentParser, "parse_args", mock_parser)
    monkeypatch.setattr(predict_cazymes, "build_logger", mock_build_logger)
    monkeypatch.setattr(predict_cazymes, "check_cwd", mock_cwd_check)
    monkeypatch.setattr(predict_cazymes, "make_output_directory", mock_making_dir)
    monkeypatch.setattr(predict_cazymes, "get_fasta_paths", mock_fasta_retrieval)
    monkeypatch.setattr(predict_cazymes, "get_protein_source", mock_protein_source_retrieval)
    monkeypatch.setattr(predict_cazymes, "get_tax_id", mock_tax_id_retrieval)
    monkeypatch.setattr(predict_cazymes, "invoke_prediction_tools", mock_prediction_tool_invoking)

    predict_cazymes.main()


def test_main_argv_tax_id_none(output_dir, null_logger, input_dir, monkeypatch):
    """Test function 'main'."""

    def mock_built_parser(*args, **kwargs):
        parser_args = ArgumentParser(
            prog="predict_cazymes.py",
            usage=None,
            description="Predict cazymes from protein sequences",
            conflict_handler="error",
            add_help=True,
        )
        return parser_args

    def mock_parser(*args, **kwargs):
        parser = Namespace(
            force=True, output=output_dir, nodelete=True, input=input_dir
        )
        return parser

    def mock_build_logger(*args, **kwargs):
        return null_logger

    def mock_cwd_check(*args, **kwargs):
        return

    def mock_making_dir(*args, **kwargs):
        return

    def mock_fasta_retrieval(*args, **kwargs):
        return ["fasta_1", "fasta_2", "fasta_3"]

    def mock_protein_source_retrieval(*args, **kwargs):
        return "protein_souce"

    def mock_tax_id_retrieval(*args, **kwargs):
        return

    def mock_prediction_tool_invoking(*args, **kwargs):
        return

    monkeypatch.setattr(predict_cazymes, "build_parser", mock_built_parser)
    monkeypatch.setattr(ArgumentParser, "parse_args", mock_parser)
    monkeypatch.setattr(predict_cazymes, "build_logger", mock_build_logger)
    monkeypatch.setattr(predict_cazymes, "check_cwd", mock_cwd_check)
    monkeypatch.setattr(predict_cazymes, "make_output_directory", mock_making_dir)
    monkeypatch.setattr(predict_cazymes, "get_fasta_paths", mock_fasta_retrieval)
    monkeypatch.setattr(predict_cazymes, "get_protein_source", mock_protein_source_retrieval)
    monkeypatch.setattr(predict_cazymes, "get_tax_id", mock_tax_id_retrieval)
    monkeypatch.setattr(predict_cazymes, "invoke_prediction_tools", mock_prediction_tool_invoking)

    predict_cazymes.main(["argv"])

# argsv is not None
# tax id is None