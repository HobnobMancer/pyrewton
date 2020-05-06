#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import logging
import unittest

from Bio import Entrez
import pytest

from Section1_Extracting_Genomes import Extract_genomes_NCBI

Entrez.email = "my.email@my.domain"


class Test_call_to_AssemblyDb(unittest.TestCase):

    """Class defining tests of Extract_genomes_NCBI.py assembly retrieval"""

    # Establish inputs for tests and expected outputs

    def setUp(self):
        """"Null logger instance"""
        self.logger = logging.getLogger("Test_name_and_ID_Retrieval logger")
        self.logger.addHandler(logging.NullHandler())

        # define assembly ID list
        self.assembly_id_list = "6106631,6106621"

    # Define tests

    @pytest.mark.run(order=7)
    def test_assembly_id_retrieval(self):
        """Tests is turned after get_accession_numbers() Entrez calls"""

        Extract_genomes_NCBI.get_accession_numbers("5061", self.logger)

    @pytest.mark.run(order=8)
    def test_assembly_id_retry(self):
        """Tests function to retry Entrez call after network error encountered.

        Input: '5061' taxonomy ID.
        Test no errors encountered.
        """

        Extract_genomes_NCBI.get_a_id_retry("5061", self.logger)

    @pytest.mark.run(order=9)
    def test_assembly_posting_and_accession_retry(self):
        """Tests function to retry Entrez call after network error encountered.

        Input: assembly IDs'6106631,6106621'; taxonomy ID '5061'
        Test no error encountered.
        """

        record = Extract_genomes_NCBI.post_a_ids_retry(
            "6106631,6106621", self.logger, "5061"
        )

        Extract_genomes_NCBI.get_a_n_retry(
            record["QueryKey"], record["WebEnv"], self.logger, "5061"
        )
