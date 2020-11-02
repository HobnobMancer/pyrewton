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

"""Tests for pyrewton genbank:get_ncbi_genomes submodule units.

These tests are inteded to be run from the root repository using:
pytest -v

To create new xml files to mock Entrez call modify the eutils URL
to perform the search, then save the resulting webpage
(this will automatically save it as a xml file)
https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?
dbfrom=Taxonomy&id=5061&db=Assembly&linkname=taxonomy_assembly
In the above URL separate out the search criteria used by '&'.
"""

import json
import pytest
import sys
import urllib.request

import pandas as pd

from argparse import Namespace, ArgumentParser
from urllib.error import URLError

from Bio import Entrez
from Bio.Entrez import Parser

from pyrewton.genbank.get_ncbi_genomes import get_ncbi_genomes

# Define dummy email for Entrez
Entrez.email = "my.email@my.domain"


@pytest.fixture
def entrez_dtd_dir(test_dir):
    """Fixes CircleCI failure to retrieve DTDs for Entrez.Parser from NCBI.
    
    See BioPython GitHub repository for demonstration of ability to pass
    Entrez.Parser a custom directory containting DTDs. Found in
    Tests.test_Entrez.py classCustomDirectoryTest:
    https://github.com/biopython/biopython/blob/master/Tests/test_Entrez.py
    """
    dir_path = test_dir / "test_inputs" / "Bio" / "Entrez" / "DTDs"
    handler = Parser.DataHandler(validate=False, escape=False)
    handler.directory = dir_path
    return


@pytest.fixture
def gt_ncbi_gnms_input_dir(test_dir):
    input_dir = test_dir / "test_inputs" / "gt_ncbi_gnms_test_inputs"
    return input_dir


@pytest.fixture
def gt_ncbi_gnms_test_inputs(gt_ncbi_gnms_input_dir):
    with (gt_ncbi_gnms_input_dir / "gt_ncbi_gnms_test_inputs.json").open("r") as ifh:
        test_inputs = json.load(ifh)
        inputs = []
        inputs.append(test_inputs["retries"])
        inputs.append(test_inputs["taxonomy_id"])
        inputs.append(test_inputs["genus_species_name"])
        inputs.append(test_inputs["line_number"])
        inputs.append(test_inputs["df_row_genus"])
        inputs.append(test_inputs["df_row_species"])
    return inputs


# create fixtures for testing parsing of input file


@pytest.fixture
def test_ncbi_species_file(gt_ncbi_gnms_input_dir):
    input_reading_path = gt_ncbi_gnms_input_dir / "gt_ncbi_gnms_reading_test_input.txt"
    return input_reading_path


@pytest.fixture
def tax_id_line():
    """Example line from input file."""
    return "NCBI:txid5061"


@pytest.fixture
def species_line():
    """Example line from onput file containing scientific name."""
    return "Aspergillus niger"


# create fixtures for testing scientific name and tax ID retrieval


@pytest.fixture
def efetch_result(gt_ncbi_gnms_input_dir):
    """Result when retrieving scientific name from NCBI."""
    results_file = gt_ncbi_gnms_input_dir / "efetch_result.xml"
    return results_file


@pytest.fixture
def efetch_result_empty(gt_ncbi_gnms_input_dir):
    """Result when retrieving scientific name and no name returned from NCBI."""
    results_file = gt_ncbi_gnms_input_dir / "efetch_result_empty.xml"
    return results_file


@pytest.fixture
def esearch_result(gt_ncbi_gnms_input_dir):
    """Result when searching for tax ID from NCBI."""
    results_file = gt_ncbi_gnms_input_dir / "esearch_result.xml"
    return results_file


@pytest.fixture
def esearch_result_empty(gt_ncbi_gnms_input_dir):
    """Result when no tax ID returned when call NCBI."""
    results_file = gt_ncbi_gnms_input_dir / "esearch_result_empty.xml"
    return results_file


# create fixtures for testing retrieval of accession numbers


@pytest.fixture
def na_df_row():
    """Create list to present pandas series containing only 'NA'."""
    mock_df_row = ["NA", "NA", "NA"]
    return mock_df_row


@pytest.fixture
def input_ncbi_df(gt_ncbi_gnms_test_inputs):
    row_data = []
    row_data.append(gt_ncbi_gnms_test_inputs[4])
    row_data.append(gt_ncbi_gnms_test_inputs[5])
    row_data.append(gt_ncbi_gnms_test_inputs[1])
    return row_data


@pytest.fixture
def accession_df(input_ncbi_df):
    row_data = [input_ncbi_df]
    df = pd.DataFrame(row_data, columns=["Genus", "Species", "NCBI Taxonomy ID"])
    return df


@pytest.fixture
def ncbi_args(test_dir, test_ncbi_species_file):
    argsdict = {
        "args": Namespace(
            genbank=True,
            retries=10,
            timeout=10,
            output=test_dir,
            input_file=test_ncbi_species_file,
        )
    }
    return argsdict


@pytest.fixture
def ncbi_args_file_exists(test_dir, test_ncbi_species_file):
    path = test_dir / "test_targets" / "gt_ncbi_gnms_test_targets2"
    path = path / "GCA_001599495_1_JCM_2005_assembly_v001_genomic.gbff.gz"
    argsdict = {
        "args": Namespace(
            genbank=True,
            retries=10,
            timeout=10,
            output=path,
            input_file=test_ncbi_species_file,
        )
    }
    return argsdict



@pytest.fixture
def ncbi_args_stdout(test_ncbi_species_file):
    argsdict = {
        "args": Namespace(
            genbank=True,
            retries=10,
            timeout=10,
            output=sys.stdout,
            input_file=test_ncbi_species_file,
        )
    }
    return argsdict


@pytest.fixture
def elink_result(gt_ncbi_gnms_input_dir):
    results_file = gt_ncbi_gnms_input_dir / "elink_result.xml"
    return results_file


@pytest.fixture
def elink_result_empty(gt_ncbi_gnms_input_dir):
    results_file = gt_ncbi_gnms_input_dir / "elink_result_empty.xml"
    return results_file


@pytest.fixture
def epost_result(gt_ncbi_gnms_input_dir):
    results_file = gt_ncbi_gnms_input_dir / "epost_result.xml"
    return results_file


@pytest.fixture
def efetch_accession_result(gt_ncbi_gnms_input_dir):
    results_file = gt_ncbi_gnms_input_dir / "efetch_accession_result.xml"
    return results_file


@pytest.fixture
def efetch_accession_result_empty(gt_ncbi_gnms_input_dir):
    results_file = gt_ncbi_gnms_input_dir / "efetch_accession_result_empty.xml"
    return results_file


@pytest.fixture
def mocked_webenv():
    webenv = ["123", "456"]
    return webenv


# create fixture to store test targets


@pytest.fixture
def gt_ncbi_gnms_targets(test_dir):
    target_dir = test_dir / "test_targets" / "gt_ncbi_gnms_test_targets"
    with (target_dir / "gt_ncbi_gnms_test_targets.json").open("r") as ifht:
        test_targets = json.load(ifht)
        target_sci_name = test_targets["target_genus_species_name"]
        target_tax_id = test_targets["target_taxonomy_id"]
    return target_sci_name, target_tax_id, target_dir


@pytest.fixture
def df_output(gt_ncbi_gnms_targets):
    target_dir = gt_ncbi_gnms_targets[2]
    output_dir = target_dir / "written_out_csv"
    return output_dir


@pytest.fixture
def coordination_args(df_output, test_ncbi_species_file):
    argsdict = {
        "args": Namespace(
            input_file=test_ncbi_species_file,
            retries=10,
            dataframe=df_output,
            force=True,
            nodelete=False,
        )
    }
    return argsdict


@pytest.fixture
def coordination_args_stdout(test_ncbi_species_file):
    argsdict = {
        "args": Namespace(
            input_file=test_ncbi_species_file,
            retries=10,
            dataframe=sys.stdout,
            force=True,
            nodelete=False,
        )
    }
    return argsdict


@pytest.fixture
def expected_web_env():
    expected_result = (
        "NCID_1_62192438_130.14.22.76_9001_1595074126_1125068541_0MetA0_S_MegaStore",
        "1",
    )
    return expected_result


@pytest.fixture
def expected_accesions():
    expected_result = "GCA_011316255.1"
    return expected_result


@pytest.fixture
def output_dir(test_dir):
    path = test_dir / "test_targets" / "gt_unprt_prtns_test_targets"
    return path


# Test parsing of input file


@pytest.mark.run(order=8)
def test_reading_input_file(
    test_ncbi_species_file, null_logger, gt_ncbi_gnms_test_inputs, monkeypatch
):
    def mock_line_parsing(*args, **kwargs):
        """Mock results from the function parse_link."""
        return ["genus", "species", "tax_id"]

    monkeypatch.setattr(get_ncbi_genomes, "parse_line", mock_line_parsing)

    """Tests script can open and read supplied input file."""
    get_ncbi_genomes.parse_input_file(
        test_ncbi_species_file, null_logger, gt_ncbi_gnms_test_inputs[1]
    )


@pytest.mark.run(order=9)
def test_input_file_check(test_dir, null_logger, gt_ncbi_gnms_test_inputs):
    """Test get_ncbi_genomes ability to detect when no input file available."""
    with pytest.raises(SystemExit) as pytest_wrapped_e:
        get_ncbi_genomes.parse_input_file(
            test_dir, null_logger, gt_ncbi_gnms_test_inputs
        )
    assert pytest_wrapped_e.type == SystemExit


@pytest.mark.run(order=10)
def test_line_parsing(tax_id_line, species_line, null_logger, monkeypatch):
    """Test function parse_line ability to parse lines of input file."""

    def mock_ncbi_data_retrieval(*args, **kwargs):
        """Mock results from functions with retrieve data from ncbi."""
        return "mock line"

    monkeypatch.setattr(
        get_ncbi_genomes, "get_genus_species_name", mock_ncbi_data_retrieval
    )
    monkeypatch.setattr(get_ncbi_genomes, "get_tax_id", mock_ncbi_data_retrieval)

    get_ncbi_genomes.parse_line(tax_id_line, null_logger, 1, 1)
    get_ncbi_genomes.parse_line(species_line, null_logger, 1, 1)


# Test called to Entrez to retrieve scientific name


@pytest.mark.run(order=11)
def test_scientific_name_retrieval(
    gt_ncbi_gnms_test_inputs,
    gt_ncbi_gnms_targets,
    null_logger,
    efetch_result,
    monkeypatch,
):
    """Tests Entrez call to NCBI to retrieve scientific name from taxonomy ID.
    Tests that correct output is returned from from get_genus_species_name()
    function in Extract_genomes_NCBI.py.
    """
    with open(efetch_result, "rb") as fh:
        result = fh

        def mock_entrez_sci_call(*args, **kwargs):
            """Mocks call to Entrez to retrieve scientific name."""
            return result

        monkeypatch.setattr(get_ncbi_genomes, "entrez_retry", mock_entrez_sci_call)

        assert gt_ncbi_gnms_targets[0] == get_ncbi_genomes.get_genus_species_name(
            gt_ncbi_gnms_test_inputs[1],
            null_logger,
            gt_ncbi_gnms_test_inputs[3],
            gt_ncbi_gnms_test_inputs[0],
        )


@pytest.mark.run(order=12)
def test_scientific_name_retrieval_indexerror_catch(
    gt_ncbi_gnms_test_inputs,
    null_logger,
    efetch_result_empty,
    monkeypatch,
):
    """Tests get_scientific name retrieval handling when no entry with the
    given taxonomy ID is found."""
    with open(efetch_result_empty, "rb") as fh:
        result = fh

        def mock_entrez_sci_call(*args, **kwargs):
            """Mocks call to Entrez to retrieve scientific name."""
            return result

        monkeypatch.setattr(get_ncbi_genomes, "entrez_retry", mock_entrez_sci_call)

        assert "NA" == get_ncbi_genomes.get_genus_species_name(
            gt_ncbi_gnms_test_inputs[1],
            null_logger,
            gt_ncbi_gnms_test_inputs[3],
            gt_ncbi_gnms_test_inputs[0],
        )


def test_scientific_name_retrieval_typeerror_catch(
    gt_ncbi_gnms_test_inputs,
    null_logger,
    monkeypatch,
):
    """Tests get_scientific name retrieval handling when no entry with the
    given taxonomy ID is found."""

    def mock_entrez_sci_call(*args, **kwargs):
        """Mocks call to Entrez to retrieve scientific name."""
        return

    monkeypatch.setattr(get_ncbi_genomes, "entrez_retry", mock_entrez_sci_call)

    assert "NA" == get_ncbi_genomes.get_genus_species_name(
        gt_ncbi_gnms_test_inputs[1],
        null_logger,
        gt_ncbi_gnms_test_inputs[3],
        gt_ncbi_gnms_test_inputs[0],
    )


# Test retrieval of taxonomy ID


@pytest.mark.run(order=13)
def test_taxonomy_id_retrieval(
    gt_ncbi_gnms_test_inputs,
    gt_ncbi_gnms_targets,
    null_logger,
    monkeypatch,
    esearch_result,
):
    """Tests Entrez call to NCBI to retrieve taxonomy ID from scientific name."""
    with open(esearch_result, "rb") as fh:
        result = fh

    def mock_entrez_txid_call(*args, **kwargs):
        """Mocks call to Entrez to retrieve taxonomy ID."""
        return result

    monkeypatch.setattr(get_ncbi_genomes, "entrez_retry", mock_entrez_txid_call)

    assert gt_ncbi_gnms_targets[1] == get_ncbi_genomes.get_tax_id(
        gt_ncbi_gnms_test_inputs[2],
        null_logger,
        1,
        gt_ncbi_gnms_test_inputs[0],
    )


@pytest.mark.run(order=14)
def test_tax_id_check(null_logger):
    """Tests searching of input scientific name for digits"""
    get_ncbi_genomes.get_tax_id("5061", null_logger, 1, 1)


@pytest.mark.run(order=15)
def test_tax_id_retrieval_indexerror_catch(
    gt_ncbi_gnms_test_inputs,
    null_logger,
    monkeypatch,
    esearch_result_empty,
):
    """Tests handling index Error when retrieving tax ID"""
    with open(esearch_result_empty, "rb") as fh:
        result = fh

    def mock_entrez_txid_call(*args, **kwargs):
        """Mocks call to Entrez to retrieve taxonomy ID."""
        return result

    monkeypatch.setattr(get_ncbi_genomes, "entrez_retry", mock_entrez_txid_call)

    assert "NA" == get_ncbi_genomes.get_genus_species_name(
        gt_ncbi_gnms_test_inputs[2],
        null_logger,
        gt_ncbi_gnms_test_inputs[3],
        gt_ncbi_gnms_test_inputs[0],
    )


def test_tax_id_retrieval_typeerror_catch(
    gt_ncbi_gnms_test_inputs,
    null_logger,
    monkeypatch,
):
    """Tests handling index Error when retrieving tax ID"""

    def mock_entrez_txid_call(*args, **kwargs):
        """Mocks call to Entrez to retrieve taxonomy ID."""
        return

    monkeypatch.setattr(get_ncbi_genomes, "entrez_retry", mock_entrez_txid_call)

    assert "NA" == get_ncbi_genomes.get_genus_species_name(
        gt_ncbi_gnms_test_inputs[2],
        null_logger,
        gt_ncbi_gnms_test_inputs[3],
        gt_ncbi_gnms_test_inputs[0],
    )

# Test retrieval of accession numbers from NCBI


@pytest.mark.run(order=16)
def test_df_cell_content_check(na_df_row, null_logger):
    """Test get_accession_numbers ability to catch "NA" content of cell."""
    get_ncbi_genomes.get_accession_numbers(na_df_row, null_logger, "args")


@pytest.mark.run(order=17)
def test_catch_failed_assembly_id(input_ncbi_df, null_logger, ncbi_args, monkeypatch):
    """Test catching failed retrieval of assembly IDs."""

    def mock_assembly_ids(*args, **kwargs):
        return "NA"

    monkeypatch.setattr(get_ncbi_genomes, "get_assembly_ids", mock_assembly_ids)

    assert "NA" == get_ncbi_genomes.get_accession_numbers(
        input_ncbi_df, null_logger, ncbi_args["args"]
    )


@pytest.mark.run(order=18)
def test_catch_failed_posting(input_ncbi_df, null_logger, ncbi_args, monkeypatch):
    """Test catching failed retrieval of posting assembly IDs."""

    def mock_assembly_ids(*args, **kwargs):
        ids = ["123", "456"]
        return ids

    def mock_posting_result(*args, **kwargs):
        return "NA"

    monkeypatch.setattr(get_ncbi_genomes, "get_assembly_ids", mock_assembly_ids)
    monkeypatch.setattr(get_ncbi_genomes, "post_assembly_ids", mock_posting_result)

    assert "NA" == get_ncbi_genomes.get_accession_numbers(
        input_ncbi_df, null_logger, ncbi_args["args"]
    )


@pytest.mark.run(order=19)
def test_catch_failed_accession_retrieval(
    input_ncbi_df, null_logger, ncbi_args, monkeypatch
):
    """Test catching failed retrieval of accession numbers."""

    def mock_assembly_ids(*args, **kwargs):
        ids = ["123", "456"]
        return ids

    def mock_posting_result(*args, **kwargs):
        webenv = "webenv"
        query_key = "QK"
        return webenv, query_key

    def mock_accession_retrieval(*args, **kwargs):
        result = "NA"
        return result

    monkeypatch.setattr(get_ncbi_genomes, "get_assembly_ids", mock_assembly_ids)
    monkeypatch.setattr(get_ncbi_genomes, "post_assembly_ids", mock_posting_result)
    monkeypatch.setattr(
        get_ncbi_genomes, "retrieve_accession_numbers", mock_accession_retrieval
    )

    assert "NA" == get_ncbi_genomes.get_accession_numbers(
        input_ncbi_df, null_logger, ncbi_args["args"]
    )


@pytest.mark.run(order=20)
def test_successful_accession_retrieval(
    input_ncbi_df, null_logger, ncbi_args, monkeypatch
):
    """Test successful retrieval of accession numbers."""

    def mock_assembly_ids(*args, **kwargs):
        ids = ["123", "456"]
        return ids

    def mock_posting_result(*args, **kwargs):
        webenv = "webenv"
        query_key = "QK"
        return webenv, query_key

    def mock_accession_retrieval(*args, **kwargs):
        accessions = ["abc", "def"]
        return accessions

    monkeypatch.setattr(get_ncbi_genomes, "get_assembly_ids", mock_assembly_ids)
    monkeypatch.setattr(get_ncbi_genomes, "post_assembly_ids", mock_posting_result)
    monkeypatch.setattr(
        get_ncbi_genomes, "retrieve_accession_numbers", mock_accession_retrieval
    )

    assert ["abc", "def"] == get_ncbi_genomes.get_accession_numbers(
        input_ncbi_df, null_logger, ncbi_args["args"]
    )


# Test retrieval of assembly IDs using Entrez.elink


@pytest.mark.run(order=21)
def test_failed_elink(
    input_ncbi_df, null_logger, ncbi_args, monkeypatch, elink_result_empty
):
    """Test catching of when no aseembly IDs retrieved from Entrez.elink"""
    with open(elink_result_empty, "rb") as fh:
        result = fh

        def mock_elink(*args, **kwargs):
            """mock Entre.elink when no result is returned"""
            return result

        monkeypatch.setattr(get_ncbi_genomes, "entrez_retry", mock_elink)

        assert "NA" == get_ncbi_genomes.get_assembly_ids(
            input_ncbi_df, null_logger, ncbi_args["args"]
        )


def test_no_elink(
    input_ncbi_df, null_logger, ncbi_args, monkeypatch, elink_result_empty
):
    """Test catching of when nothing returned from Entrez.elink"""

    def mock_elink(*args, **kwargs):
        """mock Entre.elink when no result is returned"""
        return

    monkeypatch.setattr(get_ncbi_genomes, "entrez_retry", mock_elink)

    assert "NA" == get_ncbi_genomes.get_assembly_ids(
        input_ncbi_df, null_logger, ncbi_args["args"]
    )


@pytest.mark.run(order=22)
def test_successful_elink(
    input_ncbi_df, null_logger, ncbi_args, elink_result, monkeypatch
):
    """Test processing of succesful elink retrieval of accession numbers."""
    with open(elink_result, "rb") as fh:
        result = fh

        def mock_elink(*args, **kwargs):
            """mock Entrez.elink"""
            return result

        monkeypatch.setattr(get_ncbi_genomes, "entrez_retry", mock_elink)

        get_ncbi_genomes.get_assembly_ids(input_ncbi_df, null_logger, ncbi_args["args"])


# Test posting of assembly IDs using Enrez.epost


@pytest.mark.run(order=23)
def test_failed_epost(input_ncbi_df, null_logger, ncbi_args, monkeypatch):
    """Test catching when nothing returned from Entrez.epost."""

    def mock_epost(*args, **kwargs):
        return

    monkeypatch.setattr(get_ncbi_genomes, "entrez_retry", mock_epost)

    assert "NA" == get_ncbi_genomes.post_assembly_ids(
        ["123", "456"], input_ncbi_df, null_logger, ncbi_args["args"]
    )


@pytest.mark.run(order=24)
def test_successful_epost(
    input_ncbi_df, null_logger, ncbi_args, monkeypatch, epost_result, expected_web_env
):
    """Test successful posting of assembly IDs using Entrez.epost."""
    with open(epost_result, "rb") as fh:
        result = fh

        def mock_epost(*args, **kwargst):
            return result

        monkeypatch.setattr(get_ncbi_genomes, "entrez_retry", mock_epost)

        assert expected_web_env == get_ncbi_genomes.post_assembly_ids(
            ["123", "456"], input_ncbi_df, null_logger, ncbi_args["args"]
        )


# Test retrieval of accession numbers using WebEnv data


@pytest.mark.run(order=25)
def test_failed_accession_retrieval(
    input_ncbi_df,
    null_logger,
    ncbi_args,
    monkeypatch,
    efetch_accession_result_empty,
    mocked_webenv,
):
    """Test handling data when no accession number contained in result from Entrez.efetch."""
    with open(efetch_accession_result_empty, "rb") as fh:
        result = fh

        def mock_efetch(*args, **kwargs):
            return result

        monkeypatch.setattr(get_ncbi_genomes, "entrez_retry", mock_efetch)

        "NA" == get_ncbi_genomes.retrieve_accession_numbers(
            mocked_webenv, input_ncbi_df, null_logger, ncbi_args["args"]
        )


def test_no_accession_retrieval(
    input_ncbi_df,
    null_logger,
    ncbi_args,
    monkeypatch,
    efetch_accession_result_empty,
    mocked_webenv,
):
    """Test handling data when nothing is returned from Entrez.efetch."""
    def mock_efetch(*args, **kwargs):
        return

    monkeypatch.setattr(get_ncbi_genomes, "entrez_retry", mock_efetch)

    "NA" == get_ncbi_genomes.retrieve_accession_numbers(
        mocked_webenv, input_ncbi_df, null_logger, ncbi_args["args"]
    )


@pytest.mark.run(order=26)
def test_successful_accession_number_retrieval(
    input_ncbi_df,
    null_logger,
    ncbi_args,
    monkeypatch,
    efetch_accession_result,
    mocked_webenv,
    expected_accesions,
):
    """Test processing of successful retrieval of accession numbers from Entrez.efetch."""
    with open(efetch_accession_result, "rb") as fh:
        result = fh

        def mock_efetch(*args, **kwargs):
            return result

        def mock_genbank_download(*args, **kwargs):
            return

        monkeypatch.setattr(get_ncbi_genomes, "entrez_retry", mock_efetch)
        monkeypatch.setattr(
            get_ncbi_genomes, "get_genbank_files", mock_genbank_download
        )

        assert expected_accesions == get_ncbi_genomes.retrieve_accession_numbers(
            mocked_webenv, input_ncbi_df, null_logger, ncbi_args["args"]
        )


# Test coordination of downloading GenBank files


@pytest.mark.run(order=27)
def test_coordinating_genbank_download_stdout(
    null_logger, ncbi_args_stdout, monkeypatch
):
    """Test coordination of GenBank file download, when output is stdout"""

    def mock_compile_url(*args, **kwargs):
        result = ["url", "filestem"]
        return result

    def mocK_download(*args, **kwargs):
        return

    monkeypatch.setattr(get_ncbi_genomes, "compile_url", mock_compile_url)
    monkeypatch.setattr(get_ncbi_genomes, "download_file", mocK_download)

    get_ncbi_genomes.get_genbank_files(
        "accession", "name", null_logger, ncbi_args_stdout["args"], suffix=".gbff"
    )


@pytest.mark.run(order=28)
def test_coordinating_genbank_download(null_logger, ncbi_args, monkeypatch):
    """Test coordination of GenBank file download, when output is not stdout"""

    def mock_compile_url(*args, **kwargs):
        result = ["url", "filestem"]
        return result

    def mocK_download(*args, **kwargs):
        return

    monkeypatch.setattr(get_ncbi_genomes, "compile_url", mock_compile_url)
    monkeypatch.setattr(get_ncbi_genomes, "download_file", mocK_download)

    get_ncbi_genomes.get_genbank_files(
        "accession", "name", null_logger, ncbi_args["args"], suffix=".gbff"
    )


# Test creaction of URL for GenBank download


@pytest.mark.run(order=29)
def test_compiling_url(null_logger):
    """Test generation of URL for downloading GenBank files."""
    get_ncbi_genomes.compile_url("test_accession", "test_name", null_logger, "suffix")


# Test downloading of a GenBank file


@pytest.mark.run(order=30)
def test_download(null_logger, ncbi_args, test_dir, monkeypatch):
    """Tests downloading of GenBank file"""

    def mock_url_open(*args, **kwargs):
        response = Namespace(info="info", get=1000)
        return response

    def mock_reading_response(*args, **kwargs):
        return "test test test"

    monkeypatch.setattr("urllib3.connectionpool.HTTPConnectionPool.urlopen", mock_url_open)
    monkeypatch.setattr("urllib3.response.HTTPResponse.read", mock_reading_response)

    get_ncbi_genomes.download_file(
        "http://www.google.com", ncbi_args["args"], test_dir, null_logger, "accession", "GenBank"
    )


def test_download_raise_urlerror(null_logger, ncbi_args, test_dir, monkeypatch):
    """Test file download when HTTP error is raised."""

    def mock_urlopen(*args, **kwargs):
        raise URLError("http://foo")

    monkeypatch.setattr("urllib3.connectionpool.HTTPConnectionPool.urlopen", mock_urlopen)

    get_ncbi_genomes.download_file(
        "http://foo", ncbi_args["args"], test_dir, null_logger, "accession", "GenBank"
    )


def test_download_file_exists(null_logger, ncbi_args_file_exists, test_dir, monkeypatch):
    """Tests downloading of GenBank file"""

    def mock_url_open(*args, **kwargs):
        return

    monkeypatch.setattr("urllib3.connectionpool.HTTPConnectionPool.urlopen", mock_url_open)

    get_ncbi_genomes.download_file(
        "http://foo", ncbi_args_file_exists["args"], test_dir, null_logger, "accession", "GenBank"
    )


# Test function main()


@pytest.mark.run(order=31)
def test_main(output_dir, accession_df, null_logger, monkeypatch):
    """Test function main()."""

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
        parser = Namespace(
            user="dummy_email@domain",
            output=output_dir,
            input_file="mock input",
            retries=10,
            dataframe="dataframe",
            force=False,
            nodelete=False,
        )
        return parser

    def mock_build_logger(*args, **kwargs):
        return null_logger

    def mock_parsed_inputfile(*args, **kwargs):
        # retrn dataframe when initially parsing input file
        df = accession_df
        return df

    def mock_making_dir(*args, **kwargs):
        return

    def mock_accession_retrieval(*args, **kwargs):
        accessions = [["123456"]]
        return accessions

    def mock_writing_out_df(*args, **kwargs):
        return

    monkeypatch.setattr(get_ncbi_genomes, "build_parser", mock_built_parser)
    monkeypatch.setattr(ArgumentParser, "parse_args", mock_parser)
    monkeypatch.setattr(get_ncbi_genomes, "build_logger", mock_build_logger)
    monkeypatch.setattr(get_ncbi_genomes, "make_output_directory", mock_making_dir)
    monkeypatch.setattr(get_ncbi_genomes, "parse_input_file", mock_parsed_inputfile)
    monkeypatch.setattr(
        get_ncbi_genomes, "get_accession_numbers", mock_accession_retrieval
    )
    monkeypatch.setattr(get_ncbi_genomes, "write_out_dataframe", mock_writing_out_df)

    get_ncbi_genomes.main()


@pytest.mark.run(order=32)
def test_main_no_email(output_dir, accession_df, null_logger, monkeypatch):
    """Test function main() when no email is given."""

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
        parser = Namespace(
            user=None, output=output_dir, input_file="mock input", retries=10
        )
        return parser

    def mock_build_logger(*args, **kwargs):
        return null_logger

    def mock_parsed_inputfile(*args, **kwargs):
        # retrn dataframe when initially parsing input file
        df = accession_df
        return df

    def mock_making_dir(*args, **kwargs):
        return

    def mock_accession_retrieval(*args, **kwargs):
        accessions = [["123456"]]
        return accessions

    def mock_writing_out_df(*args, **kwargs):
        return

    monkeypatch.setattr(get_ncbi_genomes, "build_parser", mock_built_parser)
    monkeypatch.setattr(ArgumentParser, "parse_args", mock_parser)
    monkeypatch.setattr(get_ncbi_genomes, "build_logger", mock_build_logger)
    monkeypatch.setattr(get_ncbi_genomes, "make_output_directory", mock_making_dir)
    monkeypatch.setattr(get_ncbi_genomes, "parse_input_file", mock_parsed_inputfile)
    monkeypatch.setattr(
        get_ncbi_genomes, "get_accession_numbers", mock_accession_retrieval
    )
    monkeypatch.setattr(get_ncbi_genomes, "write_out_dataframe", mock_writing_out_df)

    with pytest.raises(SystemExit) as pytest_wrapped_e:
        get_ncbi_genomes.main()
    assert pytest_wrapped_e.type == SystemExit


def test_main_argv(output_dir, accession_df, null_logger, monkeypatch):
    """Test function main() when argv is not None."""

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
        parser = Namespace(
            user="dummy_email@domain",
            output=output_dir,
            input_file="mock input",
            retries=10,
            dataframe="dataframe",
            force=False,
            nodelete=False,
        )
        return parser

    def mock_build_logger(*args, **kwargs):
        return null_logger

    def mock_parsed_inputfile(*args, **kwargs):
        # retrn dataframe when initially parsing input file
        df = accession_df
        return df

    def mock_making_dir(*args, **kwargs):
        return

    def mock_accession_retrieval(*args, **kwargs):
        accessions = [["123456"]]
        return accessions

    def mock_writing_out_df(*args, **kwargs):
        return

    monkeypatch.setattr(get_ncbi_genomes, "build_parser", mock_built_parser)
    monkeypatch.setattr(ArgumentParser, "parse_args", mock_parser)
    monkeypatch.setattr(get_ncbi_genomes, "build_logger", mock_build_logger)
    monkeypatch.setattr(get_ncbi_genomes, "make_output_directory", mock_making_dir)
    monkeypatch.setattr(get_ncbi_genomes, "parse_input_file", mock_parsed_inputfile)
    monkeypatch.setattr(
        get_ncbi_genomes, "get_accession_numbers", mock_accession_retrieval
    )
    monkeypatch.setattr(get_ncbi_genomes, "write_out_dataframe", mock_writing_out_df)

    get_ncbi_genomes.main(["argv"])


def test_main_stdout(accession_df, null_logger, monkeypatch):
    """Test function main() when writing to stdout."""

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
        parser = Namespace(
            user="dummy_email@domain",
            output=sys.stdout,
            input_file="mock input",
            retries=10,
            dataframe=sys.stdout,
            force=False,
            nodelete=False,
        )
        return parser

    def mock_build_logger(*args, **kwargs):
        return null_logger

    def mock_parsed_inputfile(*args, **kwargs):
        # retrn dataframe when initially parsing input file
        df = accession_df
        return df

    def mock_making_dir(*args, **kwargs):
        return

    def mock_accession_retrieval(*args, **kwargs):
        accessions = [["123456"]]
        return accessions

    def mock_writing_out_df(*args, **kwargs):
        return

    monkeypatch.setattr(get_ncbi_genomes, "build_parser", mock_built_parser)
    monkeypatch.setattr(ArgumentParser, "parse_args", mock_parser)
    monkeypatch.setattr(get_ncbi_genomes, "build_logger", mock_build_logger)
    monkeypatch.setattr(get_ncbi_genomes, "make_output_directory", mock_making_dir)
    monkeypatch.setattr(get_ncbi_genomes, "parse_input_file", mock_parsed_inputfile)
    monkeypatch.setattr(
        get_ncbi_genomes, "get_accession_numbers", mock_accession_retrieval
    )
    monkeypatch.setattr(get_ncbi_genomes, "write_out_dataframe", mock_writing_out_df)

    get_ncbi_genomes.main()


# test entrez_retry


def test_entry_retry(null_logger):
    """Test entrez_retry."""

    def mock_record(*args, **kwargs):
        return "test_record"

    assert "test_record" == get_ncbi_genomes.entrez_retry(null_logger, 1, mock_record)


def test_entrez_retry_none(null_logger):
    """Test entrez_retry when nothing is returned."""

    def mock_record(*args, **kwargs):
        return

    assert get_ncbi_genomes.entrez_retry(null_logger, 0, mock_record) is None
