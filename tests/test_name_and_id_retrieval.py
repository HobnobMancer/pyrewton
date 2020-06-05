#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import json
import logging
import unittest
from pathlib import Path

from Bio import Entrez
import pytest

from pyrewton.genbank.get_ncbi_genomes import get_ncbi_genomes

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
        self.input_dir = self.test_dir / "test_inputs" / "gt_ncbi_gnms_test_inputs"
        self.target_dir = self.test_dir / "test_targets" / "gt_ncbi_gnms_test_targets"

        # Null logger instance
        self.logger = logging.getLogger("Test_name_and_ID_Retrieval logger")
        self.logger.addHandler(logging.NullHandler())

        # Parse file containing test inputs
        with (self.input_dir / "gt_ncbi_gnms_test_inputs.json").open("r") as ifh:
            test_inputs = json.load(ifh)

        # Define test inputs
        self.retries = test_inputs["retries"]
        self.input_tax_id = test_inputs["taxonomy_id"]
        self.input_genus_species_name = test_inputs["genus_species_name"]
        self.input_line_number = test_inputs["line_number"]

        # Parse file containing test targets
        with (self.target_dir / "gt_ncbi_gnms_test_targets.json").open("r") as ifh:
            test_targets = json.load(ifh)

        # Degine test targets
        self.target_genus_species_name = test_targets["target_genus_species_name"]
        self.target_tax_id = test_targets["target_taxonomy_id"]

    # Define tests

    @pytest.mark.run(order=6)
    def test_scientific_name_retrieval(self):
        """Tests Entrez call to NCBI to retrieve scientific name from taxonomy ID.

        Tests that correct output is returned from from get_genus_species_name()
        function in Extract_genomes_NCBI.py.
        """

        self.assertEqual(
            self.target_genus_species_name,
            get_ncbi_genomes.get_genus_species_name(
                self.input_tax_id, self.logger, self.input_line_number, self.retries
            ),
        )

    @pytest.mark.run(order=7)
    def test_taxonomy_id_retrieval(self):
        """Tests Entrez call to NCBI to retrieve taxonomy ID from scientific name."""

        self.assertEqual(
            "NCBI:txid" + self.target_tax_id,
            get_ncbi_genomes.get_tax_id(
                self.input_genus_species_name,
                self.logger,
                self.input_line_number,
                self.retries,
            ),
        )
