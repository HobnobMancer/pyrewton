#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import logging
import unittest

from Bio import Entrez
import pytest

from Section1_Extracting_Genomes import Extract_genomes_NCBI

Entrez.email = "eemh1@standrews.ac.uk"


class TestName_and_IDRetrieval(unittest.TestCase):

    """Class defining tests of Extract_genomes_NCBI.py name
  and taxonomy ID retrieval. """

    # Establish inputs for tests and expected outputs

    def setUp(self):
        """"Null logger instance"""
        self.logger = logging.getLogger("Test_name_and_ID_Retrieval logger")
        self.logger.addHandler(logging.NullHandler())

    # Define tests

    @pytest.mark.run(order=3)
    def test_species_name_retrieval(self):
        """Tests Entrez call to NCBI to retrieve scientific name from taxonomy ID.
        
        Input: '5061'. Expected output: 'Aspergillus niger'.
        '0' represent an arbitary line numbers passed to the function in the
        original srcipt."""

        self.assertEqual(
            "Aspergillus niger",
            Extract_genomes_NCBI.get_genus_species_name("5061", self.logger, 0),
        )

    @pytest.mark.run(order=4)
    def test_taxonomy_id_retrieval(self):
        """Tests Entrez call to NCBI to retrieve taxonomy ID from scientific name.
        
        Input: 'Aspergillus nidulans'. Expected output: '162425'.
        '0' represents an arbitary line numbers passed to the function in the
        original script."""

        self.assertEqual(
            "NCBI:txid5061",
            Extract_genomes_NCBI.get_tax_ID("Aspergillus niger", self.logger, 0),
        )
