#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pytest

# Establish inputs for tests and expected outputs

@pytest.fixture
def tax_id():
  id = 5061
  
@pytest.fixture
def expected_name():
  e_name = "Aspergillus niger"

@pytest.fixture
def species_name():
  name = "Aspergillus nidulans"
  
@pytest.fixture
def expected_id():
  e_id = 162425
  
# Define tests
def test_species_name_retrieval(tax_id, expected_name):
"""Tests Entrez call to NCBI to retrieve scientific name from taxonomy ID"""

def test_taxonomy_id_retrieval(species_name, expected_id):
"""Tests Entrez call to NCBI to retrieve taxonomy ID from scientific name"""

