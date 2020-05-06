#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import logging
import unittest
from pathlib import Path

from Bio import Entrez
import pytest

from Section1_Extracting_Genomes import Extract_genomes_NCBI

Entrez.email = "my.email@my.domain"


class Test_call_to_AssemblyDb(unittest.TestCase):

    """Class defining tests of Extract_genomes_NCBI.py assembly retrieval"""

    # Establish inputs for tests and expected outputs

    def setUp(self):
        """"Retrieve inputs and targets for tests."""

        # Null logger instance
        self.logger = logging.getLogger("Test_name_and_ID_Retrieval logger")
        self.logger.addHandler(logging.NullHandler())

        # Parse file containing test inputs
        self.input_file_path = Path(
            "tests/test_inputs/test_ext_gnm_ncbi/test_inputs.txt"
        )

        with open(self.input_file_path) as file:
            input_list = file.read().splitlines()

        # Define test inputs
        for line in input_list:
            if line.startswith("input_taxonomy_id:"):
                self.input_tax_id = line[18:]
            elif line.startswith("input_assembly_id_list:"):
                self.input_assembly_id_list = line[23:]

    # Define tests

    @pytest.mark.run(order=7)
    def test_assembly_id_retrieval(self):
        """Tests multiplpe Entrez calls to NCBI to retrieve accession numbers.
        
        Tests all Entrez called made in Extract_genomes_NCBI.py
        get_accession_numbers() function. These are a call using a taxonomy ID
        to retrieve all assemmbly IDs, posting of all assembly IDs and retrieval
        of accession numbers from posting. Tests these can be performed without
        encountering errors. Does not test against expected/target outcome
        becuase output data is subject to potentially frequent updates.
        """

        Extract_genomes_NCBI.get_accession_numbers(self.input_tax_id, self.logger)

    @pytest.mark.run(order=8)
    def test_assembly_id_retry(self):
        """Tests function to retry Entrez call after network error encountered.
        
        Tests call to NCBI Assembly database using Entrez to retrieve all
        associated assembly IDs, for a given taxonomy ID. Tests function
        can be performed without encoutinring errors. Does not test
        against expected/target outcome becuase output data is subject to
        potentially frequent updates.
        """

        Extract_genomes_NCBI.get_a_id_retry(self.input_tax_id, self.logger)

    @pytest.mark.run(order=9)
    def test_assembly_posting_and_accession_retry(self):
        """Tests function to retry Entrez call after network error encountered.

        Tests call to NCBI using Entrez to post provided assembly IDs,
        and retrieve associated accession numbrs. Tests function
        can be performed without encoutinring errors. Does not test
        against expected/target outcome becuase output data is subject to
        potentially frequent updates.
        """
        
        record = Extract_genomes_NCBI.post_a_ids_retry(
            self.input_assembly_id_list, self.logger, self.input_tax_id
        )

        Extract_genomes_NCBI.get_a_n_retry(
            record["QueryKey"], record["WebEnv"], self.logger, self.input_tax_id
        )
