#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import logging
import unittest
from pathlib import Path

from Bio import Entrez
import pytest

from Section1_Extracting_Genomes import Extract_genomes_NCBI

# Define dummy email for Entrez
Entrez.email = "my.email@my.domain"


class TestName_and_IDRetrieval(unittest.TestCase):

    """Class defining tests of Extract_genomes_NCBI.py name
  and taxonomy ID retrieval. """

    # Establish inputs for tests and expected outputs

    def setUp(self):
        """"Set attributes for tests."""

        # Define test directories
        self.test_dir = Path("tests")
        self.input_dir = self.test_dir / "test_inputs" / "test_ext_gnm_ncbi"
        self.target_dir = self.test_dir / "test_targets" / "test_ext_gnm_ncbi"

        # Null logger instance
        self.logger = logging.getLogger("Test_name_and_ID_Retrieval logger")
        self.logger.addHandler(logging.NullHandler())

        # Parse file containing test inputs
        self.input_file_path = self.input_dir / "test_inputs.txt"

        with open(self.input_file_path) as file:
            input_list = file.read().splitlines()

        # Define test inputs
        for line in input_list:
            if line.startswith("input_taxonomy_id:"):
                self.input_tax_id = line[18:]
            elif line.startswith("input_genus_species_name:"):
                self.input_genus_species_name = line[25:]
            elif line.startswith("input_line_number:"):
                self.input_line_number = 0

        # Parse file containing test targets
        self.target_file_path = self.target_dir / "test_targets.txt"

        with open(self.target_file_path) as file:
            target_list = file.read().splitlines()

        # Degine test targets
        for line in target_list:
            if line.startswith("target_genus_species_name:"):
                self.target_genus_species_name = line[26:]
            elif line.startswith("target_taxonomy_id:"):
                self.target_tax_id = line[19:]

    # Define tests

    @pytest.mark.run(order=3)
    def test_species_name_retrieval(self):
        """Tests Entrez call to NCBI to retrieve scientific name from taxonomy ID.
        
        Tests that correct output is returned from from get_genus_species_name()
        function in Extract_genomes_NCBI.py.
        """

        self.assertEqual(
            self.target_genus_species_name,
            Extract_genomes_NCBI.get_genus_species_name(
                self.input_tax_id, self.logger, self.input_line_number
            ),
        )

    @pytest.mark.run(order=4)
    def test_species_name_network_retry(self):
        """Tests function to retry Entrez call after network error encountered.

        Tests that correct output is returned from from 
        get_genus_species_name_retry() function in Extract_genomes_NCBI.py.
        """

        self.assertEqual(
            self.target_genus_species_name,
            Extract_genomes_NCBI.get_genus_species_name_retry(
                self.input_tax_id, self.logger, self.input_line_number
            )[0]["ScientificName"],
        )

    @pytest.mark.run(order=5)
    def test_taxonomy_id_retrieval(self):
        """Tests Entrez call to NCBI to retrieve taxonomy ID from scientific name."""

        self.assertEqual(
            "NCBI:txid" + self.target_tax_id,
            Extract_genomes_NCBI.get_tax_ID(
                self.input_genus_species_name, self.logger, self.input_line_number
            ),
        )

    @pytest.mark.run(order=6)
    def test_tax_id_network_retry(self):
        """Tests function to retry Entrez call after network error encountered."""

        self.assertEqual(
            self.target_tax_id,
            Extract_genomes_NCBI.get_tax_id_retry(
                self.input_genus_species_name, self.logger, self.input_line_number
            )["IdList"][0],
        )
