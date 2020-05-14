#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import logging
import unittest
from pathlib import Path

from Bio import Entrez
import pytest

from Section1_Extracting_Genomes import Extract_genomes_NCBI

Entrez.email = "proteng.ext_gnm_ncbi@my.domain"


class Test_parsing_input_file(unittest.TestCase):

    """Class defining tests of Extract_genomes_NCBI.py ability to parse
    provided input"""

    # Establish inputs for tests and expected outputs

    def setUp(self):
        """"Set up null logger instance and input file path."""

        # Define test directories
        self.test_dir = Path("tests")
        self.input_dir = (
            self.test_dir / "test_inputs" / "test_ext_gnm_ncbi" / "input_reading"
        )

        # Null logger instance
        self.logger = logging.getLogger("Test_name_and_ID_Retrieval logger")
        self.logger.addHandler(logging.NullHandler())

        # Path to example input file provided with the programme
        self.input_reading_path = self.input_dir / "test_input_file.txt"

        # Retrieve 'retries' argument from test inputs
        self.input_file_path = self.input_dir / "test_inputs.txt"

        with open(self.input_file_path) as file:
            input_list = file.read().splitlines()

        for line in input_list:
            if line.startswith("retries"):
                self.retries = line[-1:]

    # Define tests

    @pytest.mark.run(order=1)
    def test_example_input_file_exists(self):
        """Tests that test input file, supplied with programme is present."""

        self.assertTrue(self.input_reading_path.is_file())

    @pytest.mark.run(order=2)
    def test_input_file(self):
        """Tests script can open and read supplied input file."""

        Extract_genomes_NCBI.parse_input_file(
            self.input_reading_path, self.logger, self.retries
        )
