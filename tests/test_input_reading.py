#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import logging
import unittest
from pathlib import Path

from Bio import Entrez
import pytest

from Section1_Extracting_Genomes import Extract_genomes_NCBI

Entrez.email = "my.email@my.domain"


class Test_parsing_input_file(unittest.TestCase):

    """Class defining tests of Extract_genomes_NCBI.py ability to parse
    provided input"""

    # Establish inputs for tests and expected outputs

    def setUp(self):
        """"Set up null logger instance and input file path."""

        # Null logger instance
        self.logger = logging.getLogger("Test_name_and_ID_Retrieval logger")
        self.logger.addHandler(logging.NullHandler())

        # Path to example input file provided with the programme
        self.input_file_path = Path("test_inputs/input_reading/test_input_file.txt")

    # Define tests

    @pytest.mark.run(order=1)
    def test_example_input_file_exists(self):
        """Tests that test input file, supplied with programme is present."""
        self.assertTrue(self.input_file_path.is_file())

    @pytest.mark.run(order=2)
    def test_input_file(self):
        """Tests script open and read supplied input file."""
        Extract_genomes_NCBI.parse_input_file(self.input_file_path, self.logger)
