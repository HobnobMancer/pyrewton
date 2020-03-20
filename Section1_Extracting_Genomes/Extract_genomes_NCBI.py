#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Script pulls down assembles from NCBI database with genus/species names and taxonomy ID input

from Bio import Entrez
import datetime
import re
import pandas as pd

# Fill out pull down form:
Entrez.email = "eemh1@st-andrews.ac.uk"  # enter email address
date = datetime.datetime.now()
date_of_pulldown = date.strftime("%Y-%m-%d")

# open wokring_species_list.txt and extract lines, without new line character
with open("working_species_list.txt") as file:  # use file to refer to the file object
    input_list = file.read().splitlines()

# create lists to store taxonomy IDs and genus/species names
input_tax_id_list = []  # store taxonomy IDs from input file
genus_species_names_list = []  # store genus/species names from input file
taxonomy_id_list = []  # store all taxonomic IDs

# separate out taxID and G/S names
for entry in input_list:
    if entry[0] is not "#":
        if entry.startswith("NCBI:txid", 0, 9) is True:
            input_tax_id_list.append(entry)
        else:
            genus_species_names_list.append(entry)

# the above code put into a single function
def extract_data(input_file)
"""docstring""" # to be added
    with open("working_species_list.txt") as file:
        input_list = file.read().splitlines()
        for line in input_list:
            if line[0] is not "#":
                if line.startswith("NCBI:txid", 0, 9) is True:
                    input_tax_id_list.append(line)
                else:
                    genus_species_names_list.append(line)

# input tax_id_list and genus_spcies_names_list would need to be defined globally in the main script body




# pull taxID from NCBI, using genus/species names query
for entry in genus_species_names_list:
    with Entrez.esearch(db="Taxonomy", term=entry) as handle:
        record = Entrez.read(handle)
    id = str(record["IdList"])
    id = id.replace("['", "")
    id = id.replace("']", "")
    id = "NCBI:txid" + id
    taxonomy_id_list.append(id)

# add taxID from input list to the total tax ID list
taxonomy_id_list = taxonomy_id_list + input_tax_id_list

# pull genus/species names from NCBI, using taxID query
for entry in input_tax_id_list:
    t_id = entry[9:]
    with Entrez.efetch(db="Taxonomy", id=t_id, retmode="xml") as handle:
        record = Entrez.read(handle)
genus_species_name = str(record[0]["ScientificName"])
genus_species_names_list.append(
    genus_species_name
)  # add genus/species name to the list

# Fetch accession numbers from NCBI
# create accession_numbers_list

# split genus and species names
# create genus_list
# create species_list

# create dataframe to store genus, species and taxID
species_table = pd.DataFrame(
    {
        "Genus": genus_list,
        "Species": species_list,
        "Taxonomy ID": taxonomy_id_list,
        "Accession numbers": accession_numbers_list,
    }
)

# .... Then do the assembly pull down
