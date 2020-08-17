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

"""Tests for pyrewton cazyme:uniprot:get_uniprot_proteins.

These tests are inteded to be run from the root repository using:
pytest -v
"""

import pytest

import pandas as pd

from argparse import Namespace, ArgumentParser
from bioservices import UniProt
from collections import namedtuple

from pyrewton.cazymes.uniprot import get_uniprot_proteins


# create pytest fixtures for tests


@pytest.fixture
def tax_ids():
    id_list = ["5061", "51453"]
    return id_list


@pytest.fixture
def queries():
    queries_list = ["query_1"]
    return queries_list


@pytest.fixture
def mock_uniprot_result():
    row_data = [[1, 2, 3], [4, 5, 6]]
    df = pd.DataFrame(row_data, columns=["Genus", "Species", "NCBI Taxonomy ID"])
    return df


@pytest.fixture
def input_path(test_dir):
    input_dir = test_dir / "test_inputs" / "gt_unprt_prtns_test_inputs"
    return input_dir


@pytest.fixture
def args_config_no_tax_ids(input_path):
    """Args containing path to config file with only user defined queries."""
    path = input_path / "config_no_tax_ids.yaml"
    argsdict = {"args": Namespace(input=path)}
    return argsdict


@pytest.fixture
def args_config_no_queries(input_path):
    """Args containing path to config file with only tax IDs."""
    path = input_path / "config_no_queries.yaml"
    argsdict = {"args": Namespace(input=path)}
    return argsdict


@pytest.fixture
def args_config_no_tax_ids_queries(input_path):
    """Args containing path toconfig file with both tax IDs and user defined queries."""
    path = input_path / "config_tax_ids_queries.yaml"
    argsdict = {"args": Namespace(input=path)}
    return argsdict


@pytest.fixture
def args_build_uniprot_df():
    argsdict = {"args": Namespace(outdir=None, force=False)}
    return argsdict


@pytest.fixture
def raw_uniprot_result():
    result = (
        "Organism ID     Organism        Entry   Entry name      "
        "Protein names   Length  Mass    Domains Domain count    Protein families    "
        "Gene ontology IDs       Gene ontology (molecular function)      "
        "Gene ontology (biological process) Sequence\n"
        "425011  Aspergillus niger (strain CBS 513.88 / FGSC A1513)      "
        "A2R5N0  EGLD_ASPNC      "
        "Probable endo-beta-1,4-glucanase D (Endoglucanase D) (EC 3.2.1.4) "
        "(Carboxymethylcellulase D) (Cellulase D)  "
        "412     41,981          0       Glycosyl hydrolase 61 family        "
        "GO:0005576; GO:0008810; GO:0030245; GO:0030248  cellulase activity "
        "[GO:0008810]; cellulose binding [GO:0030248]     "
        "cellulose catabolic process [GO:0030245]        "
        "MKTTTYSLLALAAASKLASAHTTVQAVWINGEDQGLGNSADGYIRSPPSNSPVTDVTSTDMTCNVNGDQAASKTLSVKA"
        "GDVVTFEWHHSDRSDSDDIIASSHKGPVQVYMAPTAKGSNGNNWVKIAEDGYHKSSDEWATDILIANKGKHNITVPDVP"
        "AGNYLFRPEIIALHEGNREGGAQFYMECVQFKVTSDGSSELPSGVSIPGVYTATDPGILFDIYNSFDSYPIPGPDVWDG"
        "SSSGSSSGSSSAAAAATTSAAVAATTPATQAAVEVSSSAAAVVESTSSAAAATTEAAAPVVSSAAPVQQATSAVTSQAQ"
        "APTTFATSSKSSKTACKNKTKSKSKVAASSTEAVVAPAPTSSVVPAVSASASASAGGVAKMYERCGGINHTGPTTCESG"
        "SVCKKWNPYYYQCVASQ"
    )
    return result


@pytest.fixture
def unformated_uniprot_results(input_path):
    path = input_path / "mocked_uniprot_results.csv"
    df = pd.read_csv(path, header=0, index_col=0)
    return df


@pytest.fixture
def formated_uniprot_results(input_path):
    path = input_path / "mocked_formated_uniprot_results.csv"
    df = pd.read_csv(path, header=0, index_col=0)
    return df


@pytest.fixture
def output_dir(test_dir):
    path = test_dir / "test_targets" / "gt_unprt_prtns_test_targets"
    return path


@pytest.fixture
def args_fasta(output_dir):
    argsdict = {"args": Namespace(fasta=True, outdir=output_dir)}
    return argsdict


@pytest.fixture
def df_series(unformated_uniprot_results):
    df_row = unformated_uniprot_results.iloc[0]
    return df_row


@pytest.fixture
def df_series_for_writing(input_path):
    path = input_path / "mocked_full_formated_uniprot_results.csv"
    df = pd.read_csv(path, header=0, index_col=0)
    df_row = df.iloc[0]
    return df_row


# Test function 'configuration'


def test_config_no_tx_ids(
    queries, null_logger, mock_uniprot_result, args_config_no_tax_ids, monkeypatch
):
    """Test config function when no tax ids given."""

    def mock_get_config_data(*args, **kwargs):
        return None, queries

    def mock_build_uniprot_df(*args, **kwargs):
        return mock_uniprot_result

    monkeypatch.setattr(get_uniprot_proteins, "get_config_data", mock_get_config_data)
    monkeypatch.setattr(get_uniprot_proteins, "build_uniprot_df", mock_build_uniprot_df)

    get_uniprot_proteins.read_configuration(args_config_no_tax_ids["args"], null_logger)


def test_config_no_queries(
    tax_ids, null_logger, args_config_no_tax_ids, mock_uniprot_result, monkeypatch
):
    """Test config function when no queries given."""

    def mock_get_config_data(*args, **kwargs):
        return tax_ids, None

    def mock_build_uniprot_df(*args, **kwargs):
        return mock_uniprot_result

    monkeypatch.setattr(get_uniprot_proteins, "get_config_data", mock_get_config_data)
    monkeypatch.setattr(get_uniprot_proteins, "build_uniprot_df", mock_build_uniprot_df)

    get_uniprot_proteins.read_configuration(args_config_no_tax_ids["args"], null_logger)


def test_config_ids_and_queries(
    tax_ids,
    queries,
    null_logger,
    mock_uniprot_result,
    args_config_no_tax_ids,
    monkeypatch,
):
    """Test config function when tax ids and queries given."""

    def mock_get_config_data(*args, **kwargs):
        return tax_ids, queries

    def mock_build_uniprot_df(*args, **kwargs):
        return mock_uniprot_result

    monkeypatch.setattr(get_uniprot_proteins, "get_config_data", mock_get_config_data)
    monkeypatch.setattr(get_uniprot_proteins, "build_uniprot_df", mock_build_uniprot_df)

    get_uniprot_proteins.read_configuration(args_config_no_tax_ids["args"], null_logger)


# test function "get_config_data"


def test_get_config_no_tax_ids(args_config_no_tax_ids, null_logger):
    """Test 'get_config_data' when only user defined queries are provided."""
    get_uniprot_proteins.get_config_data(null_logger, args_config_no_tax_ids["args"])


def test_get_config_no_queries(args_config_no_queries, null_logger):
    """Test 'get_config_data' when only tax IDs are provided."""
    get_uniprot_proteins.get_config_data(null_logger, args_config_no_queries["args"])


def test_get_config_no_tax_or_queries(args_config_no_tax_ids_queries, null_logger):
    """Test 'get_config_data' when no tax IDs or user defined queries are provided."""
    with pytest.raises(SystemExit) as pytest_wrapped_e:
        get_uniprot_proteins.get_config_data(
            null_logger, args_config_no_tax_ids_queries["args"]
        )
    assert pytest_wrapped_e.type == SystemExit


# test function 'build_uniprot_df'


def test_build_unprt_df_no_tax_id(
    queries, mock_uniprot_result, null_logger, args_build_uniprot_df, monkeypatch
):
    """Test build_uniprot_df when no tax id provided."""

    def mock_call_uniprotkb(*args, **kwargs):
        return mock_uniprot_result

    def mock_writing_out_df(*args, **kwargs):
        return

    UniProtQuery = namedtuple("UniProtQuery", "query taxid")

    monkeypatch.setattr(get_uniprot_proteins, "call_uniprotkb", mock_call_uniprotkb)
    monkeypatch.setattr(
        get_uniprot_proteins, "format_search_results", mock_call_uniprotkb
    )
    monkeypatch.setattr(
        get_uniprot_proteins, "write_out_pre_named_dataframe", mock_writing_out_df
    )

    get_uniprot_proteins.build_uniprot_df(
        UniProtQuery(queries[0], None), null_logger, args_build_uniprot_df["args"]
    )


def test_build_unprt_df_no_query(
    tax_ids, mock_uniprot_result, null_logger, args_build_uniprot_df, monkeypatch
):
    """Test build_uniprot_df when no query provided."""

    def mock_call_uniprotkb(*args, **kwargs):
        return mock_uniprot_result

    def mock_writing_out_df(*args, **kwargs):
        return

    UniProtQuery = namedtuple("UniProtQuery", "query taxid")

    monkeypatch.setattr(get_uniprot_proteins, "call_uniprotkb", mock_call_uniprotkb)
    monkeypatch.setattr(
        get_uniprot_proteins, "format_search_results", mock_call_uniprotkb
    )
    monkeypatch.setattr(
        get_uniprot_proteins, "write_out_pre_named_dataframe", mock_writing_out_df
    )

    get_uniprot_proteins.build_uniprot_df(
        UniProtQuery(None, tax_ids[0]), null_logger, args_build_uniprot_df["args"]
    )


def test_build_unprt_df(
    tax_ids,
    queries,
    mock_uniprot_result,
    null_logger,
    args_build_uniprot_df,
    monkeypatch,
):
    """Test build_uniprot_df when tax id and query provided."""

    def mock_call_uniprotkb(*args, **kwargs):
        return mock_uniprot_result

    def mock_writing_out_df(*args, **kwargs):
        return

    UniProtQuery = namedtuple("UniProtQuery", "query taxid")

    monkeypatch.setattr(get_uniprot_proteins, "call_uniprotkb", mock_call_uniprotkb)
    monkeypatch.setattr(
        get_uniprot_proteins, "format_search_results", mock_call_uniprotkb
    )
    monkeypatch.setattr(
        get_uniprot_proteins, "write_out_pre_named_dataframe", mock_writing_out_df
    )

    get_uniprot_proteins.build_uniprot_df(
        UniProtQuery(queries[0], tax_ids[0]), null_logger, args_build_uniprot_df["args"]
    )


# test function 'call_uniprotkb'


def test_call_uniprot_success(null_logger, raw_uniprot_result, monkeypatch):
    """Test calling UniProtKB when protien data is returned."""

    def mock_uniprot_result(*args, **kwargs):
        return raw_uniprot_result

    monkeypatch.setattr(UniProt, "search", mock_uniprot_result)

    get_uniprot_proteins.call_uniprotkb("test_query", null_logger)


def test_call_uniprot_no_data(null_logger, monkeypatch):
    """Test call UniProtKB when no data is returned."""

    empty_result = ""

    def mock_uniprot_result(*args, **kwargs):
        return empty_result

    monkeypatch.setattr(UniProt, "search", mock_uniprot_result)

    get_uniprot_proteins.call_uniprotkb("test_query", null_logger)


# test function 'formate_search_results'


def test_formating_results(
    unformated_uniprot_results, null_logger, args_fasta, monkeypatch
):
    """Test 'formate_search_results'."""

    def mock_ec_numbers(*args, **kwargs):
        ec_numbers = "EC 1.2.3.4, EC 2.3.4.1, "
        return ec_numbers

    def mock_writing_fasta(*args, **kwargs):
        return

    monkeypatch.setattr(get_uniprot_proteins, "get_ec_numbers", mock_ec_numbers)
    monkeypatch.setattr(get_uniprot_proteins, "write_fasta", mock_writing_fasta)

    get_uniprot_proteins.format_search_results(
        unformated_uniprot_results, "filestem", null_logger, args_fasta["args"]
    )


def test_formating_results_error_catch(
    formated_uniprot_results, null_logger, args_fasta, monkeypatch
):
    """Test 'formate_search_results' when fails to insert new column."""

    def mock_ec_numbers(*args, **kwargs):
        ec_numbers = "EC 1.2.3.4, EC 2.3.4.1, "
        return ec_numbers

    def mock_writing_fasta(*args, **kwargs):
        return

    monkeypatch.setattr(get_uniprot_proteins, "get_ec_numbers", mock_ec_numbers)
    monkeypatch.setattr(get_uniprot_proteins, "write_fasta", mock_writing_fasta)

    get_uniprot_proteins.format_search_results(
        formated_uniprot_results, "filestem", null_logger, args_fasta["args"]
    )


# test function 'get_ec_numbers'


def test_ec_number_retrieval(df_series, null_logger):
    """Test 'get_ec_numbers'."""
    get_uniprot_proteins.get_ec_numbers(df_series, null_logger)


# test function 'write_fast'


def test_writing_out_fasta_file(df_series_for_writing, null_logger, args_fasta):
    """Test writing out fasta file."""
    get_uniprot_proteins.write_fasta(
        df_series_for_writing, "filestem", null_logger, args_fasta["args"]
    )


# test function 'main'


def test_main(output_dir, null_logger, monkeypatch):
    """Test function 'main'."""

    def mock_built_parser(*args, **kwargs):
        parser_args = ArgumentParser(
            prog="get_uniprot_proteins.py",
            usage=None,
            description="Retrieve protein data from UniProtKB",
            conflict_handler="error",
            add_help=True,
        )
        return parser_args

    def mock_parser(*args, **kwargs):
        parser = Namespace(force=True, outdir=output_dir, nodelete=True)
        return parser

    def mock_build_logger(*args, **kwargs):
        return null_logger

    def mock_making_dir(*args, **kwargs):
        return

    def mock_configuration(*args, **kwargs):
        return

    monkeypatch.setattr(get_uniprot_proteins, "build_parser", mock_built_parser)
    monkeypatch.setattr(ArgumentParser, "parse_args", mock_parser)
    monkeypatch.setattr(get_uniprot_proteins, "build_logger", mock_build_logger)
    monkeypatch.setattr(get_uniprot_proteins, "make_output_directory", mock_making_dir)
    monkeypatch.setattr(get_uniprot_proteins, "read_configuration", mock_configuration)

    get_uniprot_proteins.main()


def test_main(output_dir, null_logger, monkeypatch):
    """Test function 'main'."""

    def mock_built_parser(*args, **kwargs):
        parser_args = ArgumentParser(
            prog="get_uniprot_proteins.py",
            usage=None,
            description="Retrieve protein data from UniProtKB",
            conflict_handler="error",
            add_help=True,
        )
        return parser_args

    def mock_parser(*args, **kwargs):
        parser = Namespace(outdir=output_dir)
        return parser

    def mock_build_logger(*args, **kwargs):
        return null_logger

    def mock_making_dir(*args, **kwargs):
        return

    def mock_configuration(*args, **kwargs):
        return

    monkeypatch.setattr(get_uniprot_proteins, "build_parser", mock_built_parser)
    monkeypatch.setattr(ArgumentParser, "parse_args", mock_parser)
    monkeypatch.setattr(get_uniprot_proteins, "build_logger", mock_build_logger)
    monkeypatch.setattr(get_uniprot_proteins, "make_output_directory", mock_making_dir)
    monkeypatch.setattr(get_uniprot_proteins, "read_configuration", mock_configuration)

    get_uniprot_proteins.main()
