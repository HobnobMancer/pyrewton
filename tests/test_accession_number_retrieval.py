#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import json
import logging
import unittest

from argparse import Namespace
from pathlib import Path

import pytest

from Bio import Entrez

from pyrewton.genbank.get_ncbi_genomes import get_ncbi_genomes

Entrez.email = "proteng.ext_gnm_ncbi@my.domain"


class Test_call_to_AssemblyDb(unittest.TestCase):

    """Class defining tests of Extract_genomes_NCBI.py assembly retrieval"""

    # Establish inputs for tests and expected outputs

    def setUp(self):
        """"Retrieve inputs and targets for tests."""

        # Define test directories
        self.test_dir = Path("tests")
        self.input_dir = self.test_dir / "test_inputs" / "gt_ncbi_gnms_test_inputs"

        # Null logger instance
        self.logger = logging.getLogger("Test_ac_number_retrieval logger")
        self.logger.addHandler(logging.NullHandler())

        # Parse file containing test inputs
        with (self.input_dir / "gt_ncbi_gnms_test_inputs.json").open("r") as ifh:
            test_inputs = json.load(ifh)

        # create dataframe for test
        # Genus / Species / Taxonomy ID
        self.row_data = []
        self.row_data.append(test_inputs["df_row_genus"])
        self.row_data.append(test_inputs["df_row_species"])
        self.row_data.append(test_inputs["taxonomy_id"])

        # Define Namespace and disable genbank download
        self.argsdict = {"args": Namespace(genbank=False, retries=10)}

    # Define function to test

    @pytest.mark.run(order=13)
    def test_accession_number_retrieval(self):
        """Tests multiplpe Entrez calls to NCBI to retrieve accession numbers."""
        get_ncbi_genomes.get_accession_numbers(
            self.row_data, self.logger, self.argsdict["args"]
        )

    @pytest.mark.tun(order=14)
    def test_compiling_url(self):
        """Test generation of URL for downloading GenBank files."""
        get_ncbi_genomes.compile_url(
            "test_accession", "test_name", self.logger, "suffix"
        )
