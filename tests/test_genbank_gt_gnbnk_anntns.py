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

import pandas as pd

from argparse import Namespace

from pyrewton.genbank.get_genbank_annotations import get_genbank_annotations


@pytest.fixture
def test_input_dir(test_dir):
    input_dir = test_dir / "test_inputs" / "gt_gnbnk_anntns_test_inputs"
    return input_dir


@pytest.fixture
def test_input_df_path(test_input_dir):
    input_df = test_input_dir / "test_input_df.csv"
    return input_df


@pytest.fixture
def test_input_df(test_input_df_path):
    df = pd.read_csv(test_input_df_path)
    return df


@pytest.fixture
def annotation_df(test_input_dir):
    input_path = test_input_dir / "test_annotation_df.csv"
    input_df = pd.read_csv(input_path)
    return input_df


@pytest.fixture
def pandas_series():
    mock_pandas_series = ["genes", "species", "txid", "A, B, C"]
    return mock_pandas_series


@pytest.fixture
def test_accession():
    accession = "test_accession"
    return accession


@pytest.fixture
def gb_file_dir(test_input_dir):
    gb_dir = test_input_dir / "test_genbank_files"
    return gb_dir


@pytest.fixture
def test_gb_file(gb_file_dir):
    gb_file = gb_file_dir / "GCA_test####_genomic.gbff.gz"
    return gb_file


@pytest.fixture
def test_gb_file_no_location(gb_file_dir):
    gb_file = gb_file_dir / "GCA_test####_no_location_genomic.gbff.gz"
    return gb_file


@pytest.fixture
def test_gb_file_no_translation(gb_file_dir):
    gb_file = gb_file_dir / "GCA_test####_no_translation_genomic.gbff.gz"
    return gb_file


@pytest.fixture
def coordination_args(test_input_df_path, test_dir, gb_file_dir):
    argsdict = {
        "args": Namespace(
            output=test_dir,
            force=False,
            genbank=gb_file_dir,
            df_input=test_input_df_path,
        )
    }
    return argsdict


@pytest.fixture
def no_gb_args(test_dir, gb_file_dir):
    argsdict = {"args": Namespace(output=test_dir, force=False, genbank=test_dir)}
    return argsdict


# Test coordination of script


@pytest.mark.run(order=36)
def test_genbank_anno_coordination(
    null_logger, coordination_args, test_input_df, monkeypatch
):
    """Test coordination of GenBank protein annotation retrieval."""

    def mock_create_dataframe(*args, **kwargs):
        df = test_input_df
        return df

    def mock_write_out_dataframe(*args, **kwargs):
        return

    monkeypatch.setattr(
        get_genbank_annotations, "create_dataframe", mock_create_dataframe
    )
    monkeypatch.setattr(
        get_genbank_annotations, "write_out_dataframe", mock_write_out_dataframe
    )

    get_genbank_annotations.retrieve_genbank_annotations(
        null_logger, coordination_args["args"]
    )


# Test the creation of the dataframe


@pytest.mark.run(order=37)
def test_dataframe_creation(
    test_input_df, null_logger, coordination_args, annotation_df, monkeypatch
):
    """Test creation of dataframe."""

    def mock_annotation_retrieval(*args, **kwargs):
        df = annotation_df
        return df

    monkeypatch.setattr(
        get_genbank_annotations, "get_genbank_annotations", mock_annotation_retrieval
    )

    get_genbank_annotations.create_dataframe(
        test_input_df, coordination_args["args"], null_logger
    )


# Test retrieval of genbank annotations


@pytest.mark.run(order=38)
def test_get_annotations_no_data(
    pandas_series, null_logger, coordination_args, monkeypatch
):
    """Test coordination of annotation retrieval when no data retrieved for a given protein."""

    def mock_get_anno(*args, **kwargs):
        protein_data = []
        return protein_data

    monkeypatch.setattr(get_genbank_annotations, "get_annotations", mock_get_anno)

    get_genbank_annotations.get_genbank_annotations(
        pandas_series, coordination_args["args"], null_logger
    )


@pytest.mark.run(order=39)
def test_get_annotations_data_returned(
    pandas_series, null_logger, coordination_args, monkeypatch
):
    """Test coordination of annotation retrieval for a given protein."""

    def mock_get_anno(*args, **kwargs):
        annotations = [["NA", "NA", "NA", "NA", "NA"]]
        return annotations

    monkeypatch.setattr(get_genbank_annotations, "get_annotations", mock_get_anno)

    get_genbank_annotations.get_genbank_annotations(
        pandas_series, coordination_args["args"], null_logger
    )


# Test retrieval of accessions


@pytest.mark.run(order=40)
def test_get_annotations_na(null_logger, coordination_args):
    """Test get_annotations when accession number is 'NA'."""
    accession = "NA"
    get_genbank_annotations.get_annotations(
        accession, coordination_args["args"], null_logger
    )


@pytest.mark.run(order=41)
def test_get_annotations_file_none(
    test_accession, null_logger, coordination_args, monkeypatch
):
    """Test get_annotations when gb_file is None."""

    def mock_get_gb_file(*args, **kwargs):
        return

    monkeypatch.setattr(get_genbank_annotations, "get_genbank_file", mock_get_gb_file)

    get_genbank_annotations.get_annotations(
        test_accession, coordination_args["args"], null_logger
    )


@pytest.mark.run(order=42)
def test_get_annotations_all_data_na(
    test_gb_file, test_accession, null_logger, coordination_args, monkeypatch
):
    """Test get_annotations when all returned protein data is 'NA'."""

    def mock_get_gb_file(*args, **kwargs):
        gb_file = test_gb_file
        return gb_file

    def mock_get_record(*args, **kwargs):
        returned_data = "NA"
        return returned_data

    monkeypatch.setattr(get_genbank_annotations, "get_genbank_file", mock_get_gb_file)
    monkeypatch.setattr(get_genbank_annotations, "get_record_feature", mock_get_record)

    get_genbank_annotations.get_annotations(
        test_accession, coordination_args["args"], null_logger
    )


@pytest.mark.run(order=43)
def test_get_annotations_successful(
    test_gb_file, test_accession, null_logger, coordination_args, monkeypatch
):
    """Test get_annotations when all returned protein data is 'NA'."""

    def mock_get_gb_file(*args, **kwargs):
        gb_file = test_gb_file
        return gb_file

    monkeypatch.setattr(get_genbank_annotations, "get_genbank_file", mock_get_gb_file)

    get_genbank_annotations.get_annotations(
        test_accession, coordination_args["args"], null_logger
    )


@pytest.mark.run(order=44)
def test_get_annotations_not_5(
    test_gb_file, test_accession, null_logger, coordination_args, monkeypatch
):
    """Test get_annotations when length of protein data is not 5."""

    def mock_get_gb_file(*args, **kwargs):
        gb_file = test_gb_file
        return gb_file

    def mock_get_record(*args, **kwargs):
        return

    monkeypatch.setattr(get_genbank_annotations, "get_genbank_file", mock_get_gb_file)
    monkeypatch.setattr(get_genbank_annotations, "get_record_feature", mock_get_record)

    get_genbank_annotations.get_annotations(
        test_accession, coordination_args["args"], null_logger
    )


# Test when retrieval of location data fails


@pytest.mark.run(order=45)
def test_anno_retrieval_no_location(
    test_gb_file_no_location,
    test_accession,
    null_logger,
    coordination_args,
    monkeypatch,
):
    """Test get_record_feature when there is no location."""

    def mock_get_gb_file(*args, **kwargs):
        gb_file = test_gb_file_no_location
        return gb_file

    monkeypatch.setattr(get_genbank_annotations, "get_genbank_file", mock_get_gb_file)

    get_genbank_annotations.get_annotations(
        test_accession, coordination_args["args"], null_logger
    )


@pytest.mark.run(order=46)
def test_anno_retrieval_no_qualifier(
    test_gb_file_no_translation,
    test_accession,
    null_logger,
    coordination_args,
    monkeypatch,
):
    """Test get_record_feature when retrieval of qualifier fails."""

    def mock_get_gb_file(*args, **kwargs):
        gb_file = test_gb_file_no_translation
        return gb_file

    monkeypatch.setattr(get_genbank_annotations, "get_genbank_file", mock_get_gb_file)

    get_genbank_annotations.get_annotations(
        test_accession, coordination_args["args"], null_logger
    )


# Test gb_file retrieval


@pytest.mark.run(order=47)
def test_get_file_success(coordination_args, null_logger):
    """Test successful retrieval of single gb_file using get_gb_file."""
    accession = "GCA_test####"
    get_genbank_annotations.get_genbank_file(
        accession, coordination_args["args"], null_logger
    )


@pytest.mark.run(order=48)
def test_get_file_no_file(no_gb_args, null_logger):
    """Test get_gb_file when no files were retrieved."""
    accession = "GCA_test####"
    get_genbank_annotations.get_genbank_file(accession, no_gb_args["args"], null_logger)


@pytest.mark.run(order=49)
def test_get_file_multiple(coordination_args, null_logger):
    """Test get_gb_file when multiple files are retrieved."""
    accession = "GCA_testmultiple"
    get_genbank_annotations.get_genbank_file(
        accession, coordination_args["args"], null_logger
    )


@pytest.mark.run(order=50)
def test_get_file_empty(coordination_args, null_logger):
    """Test get_gb_file when the returned file is empty."""
    accession = "GCA_testempty"
    get_genbank_annotations.get_genbank_file(
        accession, coordination_args["args"], null_logger
    )
