#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import logging
import unittest

from argparse import Namespace
from pathlib import Path

import pandas as pd
import pytest

from Bio import Entrez

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

        with open(self.input_file_path) as ifh:
            input_list = ifh.read().splitlines()

        # Define test inputs
        self.df_row_data = []
        for line in input_list:
            if line.startswith("input_taxonomy_id:"):
                self.input_tax_id = line[18:]
            elif line.startswith("input_df_row_genus:"):
                self.df_row_data.append(line)[19:]
            elif line.startswith("input_df_row_species:"):
                self.df_row_data.append(line)[21:]

        # create dataframe for test
        # Genus / Species / Taxonomy ID
        self.test_df = pd.DataFrame(self.df_row_data, columns=["G", "S", "T_ID"])
        print(self.test_df)

        # Define Namespace and disable genbank download
        self.argsdict = {"args": Namespace(genbank=False, retries=10)}

    # Define function to test

    @pytest.mark.run(order=5)
    def test_accession_number_retrieval(self):
        """Tests multiplpe Entrez calls to NCBI to retrieve accession numbers."""
        Extract_genomes_NCBI.get_accession_numbers(
            self.test_df, self.logger, self.argsdict["args"],
        )
