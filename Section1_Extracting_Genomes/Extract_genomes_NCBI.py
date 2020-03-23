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


def main():
    """Read input file and pass genus/species names and taxonomy IDs to approriapte functions
    to pull down the accosicated taxonomy ID and genus/species name, respectively.

    Open input file, a plain text file containing comments (indicated by the first line
    character '#'), taxonomy IDs (starting with 'NCBI:txid'), and genus/species names.
    Pass over the lines in the plain text file. The line contains a taxonomy ID the associated
    genus/species name is returned by calling get_genus_species_name(), and if line contains
    a genus/species name the taxonomy ID is obtained by calling get_tax_ID().
    The genus and species names are separated, and then these names with the taxonomy ID are 
    passed to a list (called 'working_species').
    The list containing the genus name, species name and taxonomy ID from the current line
    being worked on from the input file, is passed to the tuple 'all_species_data'.

    Store genus name, species name and taxonomy ID tuplets stored in the 'all_species_data'
    tuple, which is used to create a dataframe, of three columnes: 'Genus', 'Species',
    and 'NCBI Taxonomy ID'.

    The 'NCBI Taxonomy ID' column is passed to the get_accession_numbers() function. Add
    the retrieved accession numbers to a new column in the dataframe ('NCBI Accession Numbers')
    """

    # create tuple to store genus name, species name and taxonomy ID tuplets
    all_species_data = []

    # open working_species_list.txt and extract lines, without new line character
    # then create genus name, species name and taxonomy ID tuplet
    with open("working_species_list.txt") as file:
        input_list = file.read().splitlines()
        for line in input_list:
            if line[0] is not "#":
                if line.startswith("NCBI:txid", 0, 9) is True:
                    gs_name = get_genus_species_name(line[9:])
                    line_data = gs_name.split()
                    line_data.append(line)
                    all_species_data.append(line_data)
                else:
                    tax_id = get_tax_ID(line)
                    line_data = line.split()
                    line_data.append(tax_id)
                    all_species_data.append(line_data)

    # create dataframe containg three columns: 'Genus', 'Species', 'NCBI Taxonomy ID'
    species_table = pd.DataFrame(
        all_species_data, columns=["Genus", "Species", "NCBI Taxonomy ID"]
    )

    # pull down all accession numbers associated with each Taxonomy ID, from NCBI
    all_accession_numbers = []
    for NCBI_taxonomy_id in species_table["NCBI Taxonomy ID"]:
        working_tax_id = NCBI_taxonomy_id[9:]
        working_accession_numbers = get_accession_numbers(working_tax_id)
        working_accession_numbers = str(working_accession_numbers)
        working_accession_numbers = working_accession_numbers.replace("[", "")
        working_accession_numbers = working_accession_numbers.replace("]", "")
        working_accession_numbers = working_accession_numbers.replace("'", "")
        print(working_accession_numbers)
        all_accession_numbers.append(working_accession_numbers)

    # add accession numbers to the dataframe
    species_table["NCBI Accession Numbers"] = all_accession_numbers
    print("\n", species_table, "\n")


def get_genus_species_name(taxonomy_id):
    """Pull fetch scientific name of species from NCBI,
    using its NCBI taxonomy ID as the query.

    Use Entrez efetch function to pull down the scientific name in the 
    NCBI Taxonomy database, associated with the taxonomy ID passed to the
    function.
    """

    with Entrez.efetch(db="Taxonomy", id=taxonomy_id, retmode="xml") as handle:
        record = Entrez.read(handle)
    genus_species_name = str(record[0]["ScientificName"])
    return genus_species_name


def get_tax_ID(genus_species):
    """Pull down taxonomy ID from NCBI, using genus/species name as query.

    Use Entrez esearch function to pull down the NCBI taxonomy ID of the
    species name passed to the function. Then the taxonomy ID is formatted
    and returned in the recommended referencing formate of 'NCBI:txidXXX',
    where 'XXX' is taxonomy ID.
    """

    with Entrez.esearch(db="Taxonomy", term=genus_species) as handle:
        record = Entrez.read(handle)
    id = str(record["IdList"])
    id = id.replace("['", "")
    id = id.replace("']", "")
    id = "NCBI:txid" + id
    return id


def get_accession_numbers(taxonomy_id_column):
    """Pulls down NCBI accession numbers using their taxonomy IDs as the
    query.

    Passes over a Pandas dataframe column containing NCBI taxonomy IDs. Pass
    Taxonomy IDs to Entrez eLink to pull down associated accession numbers.
    Accession numbers stored in list 'accession_numbers_list'.
    """

    taxonomy_id = taxonomy_id_column[9:]

    with Entrez.elink(
        dbfrom="Taxonomy",
        id=taxonomy_id_column,
        db="Assembly",
        linkname="taxonomy_assembly",
    ) as handle_one:
        record_one = Entrez.read(handle_one)
        id_list = []
        for id in record_one[0]["LinkSetDb"][0]["Link"]:
            id_string = str(id)
            id_string = id_string.replace("{'Id': '", "")
            id_string = id_string.replace("'}", "")
            id_list.append(id_string)

    NCBI_accession_numbers_list = []
    for assembly_db_id in id_list:
        with Entrez.efetch(
            db="Assembly", id=assembly_db_id, rettype="docsum", retmode="xml"
        ) as handle_two:
            record_two = Entrez.read(handle_two)
            NCBI_accession_number = record_two["DocumentSummarySet"]["DocumentSummary"][
                0
            ]["AssemblyAccession"]
            NCBI_accession_numbers_list.append(NCBI_accession_number)
    return NCBI_accession_numbers_list
