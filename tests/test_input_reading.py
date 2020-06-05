#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import json
import logging
import unittest
from pathlib import Path

from Bio import Entrez
import pytest

from pyrewton.genbank.get_ncbi_genomes import get_ncbi_genomes

Entrez.email = "my_email@my.domain"


class Test_parsing_input_file(unittest.TestCase):

    """Class defining tests of Extract_genomes_NCBI.py ability to parse
    provided input"""

    # Establish inputs for tests and expected outputs

    def setUp(self):
        """"Set up null logger instance and input file path."""

        # Define test directories
        self.test_dir = Path("tests")
        self.input_dir = self.test_dir / "test_inputs" / "gt_ncbi_gnms_test_inputs"

        # Null logger instance
        self.logger = logging.getLogger("Test_name_and_ID_Retrieval logger")
        self.logger.addHandler(logging.NullHandler())

        # Path to example input file provided with the programme
        self.input_reading_path = self.input_dir / "gt_ncbi_gnms_reading_test_input.txt"

        # Retrieve 'retries' argument from test inputs
        with (self.input_dir / "gt_ncbi_gnms_test_inputs.json").open("r") as ifh:
            test_inputs = json.load(ifh)
        self.retries = test_inputs["retries"]

    # Define tests

    @pytest.mark.run(order=5)
    def test_example_input_file_exists(self):
        """Tests that test input file, supplied with programme is present."""
        self.assertTrue(self.input_reading_path.is_file())

    @pytest.mark.run(order=6)
    def test_input_file(self):
        """Tests script can open and read supplied input file."""
        get_ncbi_genomes.parse_input_file(
            self.input_reading_path, self.logger, self.retries
        )
