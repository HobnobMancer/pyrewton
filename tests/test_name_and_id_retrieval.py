#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import logging
import unittest

import pytest

class TestName_and_IDRetrieval(unittest.TestCase):

  """Class defining tests of Extract_genomes_NCBI.py name
  and taxonomy ID retrieval. """

  # Establish inputs for tests and expected outputs

  @pytest.fixture
  def tax_id():
    id = 5061
    return (id)
  
  @pytest.fixture
  def expected_name():
    e_name = "Aspergillus niger"
    return (e_name)  
  
  @pytest.fixture
  def species_name():
    name = "Aspergillus nidulans"
    return (name)
    
  @pytest.fixture
  def expected_id():
    e_id = 162425
    return (e_id)
  
  # Null logger instance
  self.logger = logging.getLogger("Test_name_and_ID_Retrieval logger")
  self.logger.addHandler(logging.NullHandler())
  
  # Define tests

  def test_species_name_retrieval(self, tax_id, expected_name):
    """Tests Entrez call to NCBI to retrieve scientific name from taxonomy ID."""

    # 0 represents an arbitary line numbers passed to the function
    assert expected_name == Section1_Extracting_Genomes.Extract_genomes_NCBI:get_genus_species_name(tax_id, self.logger, 0)

  def test_taxonomy_id_retrieval(self, species_name, expected_id):
    """Tests Entrez call to NCBI to retrieve taxonomy ID from scientific name."""

    # 0 represents an arbitary line numbers passed to the function
    assert expected_id == Section1_Extracting_Genomes.Extract_genomes_NCBI: get_tax_ID(species_name, self.logger, 0)

if __name__ == '__main__'
  unittest.main()
