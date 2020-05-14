#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import logging
import unittest
from pathlib import Path

from Bio import Entrez
import pytest

from Section1_Extracting_Genomes import Extract_genomes_NCBI

Entrez.email = "proteng.ext_gnm_ncbi@my.domain"


class Test_call_to_AssemblyDb(unittest.TestCase):

    """Class defining tests of Extract_genomes_NCBI.py assembly retrieval"""

    # Establish inputs for tests and expected outputs

    def setUp(self):
        """"Retrieve inputs and targets for tests."""

        # Define test directories
        self.test_dir = Path("tests")
        self.input_dir = self.test_dir / "test_inputs" / "test_ext_gnm_ncbi"

        # Null logger instance
        self.logger = logging.getLogger("Test_name_and_ID_Retrieval logger")
        self.logger.addHandler(logging.NullHandler())

        # Parse file containing test inputs
        self.input_file_path = self.input_dir / "test_inputs.txt"

        with open(self.input_file_path) as file:
            input_list = file.read().splitlines()

        # Define test inputs
        for line in input_list:
            if line.startswith("retries"):
                self.retries = line[-1:]
            elif line.startswith("input_taxonomy_id:"):
                self.input_tax_id = line[18:]
            elif line.startswith("input_assembly_id_list:"):
                self.input_assembly_id_list = line[23:]

    # Define tests

    @pytest.mark.run(order=5)
    def test_assembly_id_retrieval(self):
        """Tests multiplpe Entrez calls to NCBI to retrieve accession numbers."""

        Extract_genomes_NCBI.get_accession_numbers(
            self.input_tax_id, self.logger, self.retries
        )
