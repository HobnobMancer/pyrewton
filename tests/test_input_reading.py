#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import logging
import unittest
from pathlib import Path

from Bio import Entrez
import pytest

from Section1_Extracting_Genomes import Extract_genomes_NCBI

Entrez.email = "eemh1@standrews.ac.uk"


class TestName_and_IDRetrieval(unittest.TestCase):

    """Class defining tests of Extract_genomes_NCBI.py ability to parse
    provided input"""

    # Establish inputs for tests and expected outputs

    def setUp(self):
        """"Set up null logger instance and input file path."""

        # Null logger instance
        self.logger = logging.getLogger("Test_name_and_ID_Retrieval logger")
        self.logger.addHandler(logging.NullHandler())

        # Path to example input file provided with the programme
        self.input_file_path = (
            Path(__file__).resolve().parent / "test_input_files" / "test_input_file.txt"
        )

    # Define tests

    def test_example_input_file_exists(self):

        """Tests that example input file, supplied with programme is present."""
        assert self.input_file_path.is_file() == True

    def test_input_file(self):
        """Tests script open and read supplied input file."""

        Extract_genomes_NCBI.parse_input_file(self.input_file_path, self.logger)
