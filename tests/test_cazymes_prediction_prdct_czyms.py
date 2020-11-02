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

import os
import pytest

import pandas as pd

from argparse import Namespace, ArgumentParser

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


@pytest.fixture
def fasta_args(input_dir):
    path = input_dir / "input_fasta"
    parser = Namespace(input=path)
    return parser


@pytest.fixture
def fasta_args_fail(input_dir):
    parser = Namespace(input=input_dir)
    return parser


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
    """Test function 'main', when argv is not None and Tax ID is None."""

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


# test check_cwd()


def test_check_cwd_pass(null_logger, monkeypatch):
    """Test check_cwd() when it should fully pass."""

    def mock_getcwd(*args, **kwargs):
        return "pyrewton/cazymes/prediction"

    monkeypatch.setattr(os, "getcwd", mock_getcwd)

    predict_cazymes.check_cwd(null_logger)


def test_check_cwd_move_dir(null_logger, monkeypatch):
    """Test check_cwd() when invokved in tools/ dir."""

    def mock_getcwd(*args, **kwargs):
        return "pyrewton/cazymes/prediction/tools"

    def mock_chdir(*args, **kwargs):
        return

    monkeypatch.setattr(os, "getcwd", mock_getcwd)
    monkeypatch.setattr(os, "chdir", mock_chdir)

    predict_cazymes.check_cwd(null_logger)


def test_check_cwd_fail(null_logger, monkeypatch):
    """Test check_cwd() when script should be terminated because in incorrect dir."""

    def mock_getcwd(*args, **kwargs):
        return "pyrewton/cazymes/"

    monkeypatch.setattr(os, "getcwd", mock_getcwd)

    with pytest.raises(SystemExit) as pytest_wrapped_e:
        predict_cazymes.check_cwd(null_logger)
    assert pytest_wrapped_e.type == SystemExit


# test get_fasta_paths()


def test_get_fasta_paths(fasta_args, null_logger):
    """Test retrieval of fasta files by get_fasta_paths."""

    assert len(predict_cazymes.get_fasta_paths(fasta_args, null_logger)) == 3


def test_get_fasta_paths_sysexit(fasta_args_fail, null_logger):
    """Test get_fasta_paths when no fasta files returned and sys.exit should be raised."""

    with pytest.raises(SystemExit) as pytest_wrapped_e:
        predict_cazymes.get_fasta_paths(fasta_args_fail, null_logger)
    assert pytest_wrapped_e.type == SystemExit


# test get_protein_source()


def test_protein_source_genomic_accession(null_logger):
    """Test get_protein_source() when the protein source is a genomic accession number."""

    file_path = "GCA_123_1.fasta"
    assert predict_cazymes.get_protein_source(file_path, null_logger) == "GCA_123_1"


def test_get_protein_source_uniprot(null_logger):
    """Test get_protein_source() when the protein source is UniProtKB."""

    file_path = "UniProt_fasta.fasta"
    assert predict_cazymes.get_protein_source(file_path, null_logger) == "UniProt_fasta"


def test_get_protein_source_unknown(null_logger):
    """Test get_protein_source() when the protein source is unknown."""

    file_path = "test_test_test.fasta"
    assert predict_cazymes.get_protein_source(file_path, null_logger) == "unknown_source"


# test get_tax_id()


def test_get_tax_id(null_logger):
    """Test get_tax_id when retrieved with ncbi-txid prefix."""

    file_path = "ncbi_txid123.fasta"
    assert predict_cazymes.get_tax_id(file_path, null_logger) == "ncbi_txid123"


def test_get_tax_id_taxonomy(null_logger):
    """Test get_tax_id when retrieved with taxonomy prefix."""

    file_path = "taxonomy__123__.fasta"
    assert predict_cazymes.get_tax_id(file_path, null_logger) == "taxonomy__123"


def test_get_tax_id_fail(null_logger):
    """Test get_tax_id when no tax ID is returned."""

    file_path = "test_test_test.fasta"
    assert predict_cazymes.get_tax_id(file_path, null_logger) == None
