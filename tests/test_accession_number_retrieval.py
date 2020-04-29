#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import logging
import unittest

from Bio import Entrez
import pytest

from Section1_Extracting_Genomes import Extract_genomes_NCBI

Entrez.email = "eemh1@standrews.ac.uk"


class TestName_and_IDRetrieval(unittest.TestCase):

    """Class defining tests of Extract_genomes_NCBI.py assembly retrieval"""

    # Establish inputs for tests and expected outputs

    def setUp(self):
        """"Null logger instance"""
        self.logger = logging.getLogger("Test_name_and_ID_Retrieval logger")
        self.logger.addHandler(logging.NullHandler())

    # Define tests

    def test_assembly_id_retrieval(self):
        """Tests Entrez call to NCBI to retrieve all accession numbers associated
        with taxonomy ID. Specifically, tests that something is returned after all
        Entrez calls"""

        # Input for Entrez call: '5061'
        # Test something was returned, and not an empty query

        Extract_genomes_NCBI.get_accession_numbers("5061", self.logger)
