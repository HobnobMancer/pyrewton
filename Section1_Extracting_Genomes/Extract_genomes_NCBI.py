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
    """Generates datafame containing genus/species name with NCBI taxonomy and accession numbers

    Pass input file (containing unique species on each line, idenfitied by their genus/species
    name or NCBI Taxnomy ID) tp parse_input_file() function, which aquires missing genus/species
    names and NCBI taxonomy IDs as approprirate. Genus/species names and associated NCBI
    Taxonomy ID are stored in a generated dataframe with three columns: 'Genus', 'Species',
    and 'NCBI Taxonomy ID'. Parse_input_file() returns the generated dataframe which is stored in
    the variable 'species_table'.

    Pass the 'species_table' dataframe to collate_accession_numbers(), which aquires all associated
    aquires all associated NCBI accession numbers associated for each Taxonomy ID. Aquired accession
    numbers are stored in a new column in the dataframe, titled 'NCBI Accession Numbers'.
    The modified dataframe is returned from collate_accession_numbers() and stored in the variable
    'species_table'.
    """

    # create dataframe storing genus, species and NCBI Taxonomy ID, called 'species_table'
    species_table = parse_input_file(
        "working_species_list.txt"
    )  # pass name of input file with extension

    # pull down all accession numbers associated with each Taxonomy ID, from NCBI
    # add accession numbers to 'species_table' dataframe
    species_table = collate_accession_numbers(species_table)
    print("\n", species_table)


def parse_input_file(input_filename):
    """Parse input file, returning dataframe of species names and NCBI Taxonomy IDs.

    Read input file passed to the function, performing the appropriate action depending
    on each line's content.

    Comments: Indicated with the first character of '#', perform no action.
    NCBI Taxonomy ID: Indicated with the first nine characters of 'NCBI:txid', pass 
    the Taxonomy ID (excluding the 'NCBI:txid' prefix) to get_genus_species_name() to 
    return associated scientific name.
    Genus/species name: Indicated by lack of '#' and 'NCBI:txid', pass genus/species
    name to get get_tax_ID() to return associated NCBI Taxonomy ID with 'NCBI:txid'
    prefix.

    Split genus and species name into separate variables. Store genus, species and
    associated Taxonomy ID in the list 'line_data'. Append 'line_data' to
    'all_species_data' list. Repeat for each line.

    Generate a dataframe with three columns: 'Genus', 'Species' and 'NCBI Taxonomy ID'.
    Use data from 'all_species_data' to fill out dataframe.

    Return dataframe.
    """

    # open working_species_list.txt and extract lines, without new line character
    # then create genus name, species name and taxonomy ID tuplet
    all_species_data = []
    with open(input_filename) as file:
        input_list = file.read().splitlines()
        print("Reading input file, and acquiring Tax IDs and Genus/Species names")
        lines = len(input_list)
        line_count = 1
        for line in input_list:
            if line[0] != "#":
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
            print("line", line_count, "/", lines)
            line_count += 1

    # create dataframe containg three columns: 'Genus', 'Species', 'NCBI Taxonomy ID'
    print("Generating genus, species and taxonomy ID dataframe")
    species_table = pd.DataFrame(
        all_species_data, columns=["Genus", "Species", "NCBI Taxonomy ID"]
    )
    print("Dataframe completed")
    return species_table


def get_genus_species_name(taxonomy_id):
    """Fetch scientfic name associated with the NCBI Taxonomy ID.

    Use Entrez efetch function to pull down the scientific name (genus/species name)
    in the     NCBI Taxonomy database, associated with the taxonomy ID passed to the
    function.
    """

    with Entrez.efetch(db="Taxonomy", id=taxonomy_id, retmode="xml") as handle:
        record = Entrez.read(handle)
    genus_species_name = str(record[0]["ScientificName"])
    return genus_species_name


def get_tax_ID(genus_species):
    """Pull down taxonomy ID from NCBI, using genus/species name as query.

    Use Entrez esearch function to pull down the NCBI Taxonomy ID of the
    species name passed to the function. Return the NCBI Taxonomy ID with
    the prefix 'NCBI:txid'.
    """

    with Entrez.esearch(db="Taxonomy", term=genus_species) as handle:
        record = Entrez.read(handle)
    fetched_taxonomy_id = "NCBI:txid" + "".join(
        re.findall("\d.*?\d", str(record["IdList"]))
    )
    return fetched_taxonomy_id


def collate_accession_numbers(species_table):
    """Return dataframe with column containing all associated NCBI accession numbers.

    Pass each Taxonomy ID in the 'NCBI Taxonomy ID' column of the 'species_table'
    dataframe, find all associated NCBI accession numbers by calling get_accession_numbers().

    Format object returned from get_accession_numbers into a human readable list. Append list
    of accession numbers to 'all_accession_numbers' list. Create fourth column, called
    'NCBI Accession Numbers' in the 'species_table' dataframe and populate with data in
    'all_accession_numbers'. 
    
    Return modified dataframe.
    """

    all_accession_numbers = []
    tax_id_total_count = len(species_table["NCBI Taxonomy ID"])
    tax_id_counter = 1
    print("Aquiring accession numbers from NCBI Assembly database")
    for NCBI_taxonomy_id in species_table["NCBI Taxonomy ID"]:
        working_tax_id = NCBI_taxonomy_id[9:]
        working_accession_numbers = ", ".join(
            re.findall("GCA_\d*.\d", get_accession_numbers(working_tax_id))
        )
        all_accession_numbers.append(working_accession_numbers)
        print("Species", tax_id_counter, "/", tax_id_total_count)
        tax_id_counter += 1

    # add accession numbers to the dataframe
    print("Adding accession numbers to dataframe")
    species_table["NCBI Accession Numbers"] = all_accession_numbers
    return species_table


def get_accession_numbers(taxonomy_id_column):
    """Return all NCBI accession numbers associated with NCBI Taxonomy Id.

    Use Entrez elink function to pull down the assembly IDs of all genomic assemblies
    stored in the NCBI Assembly database that are associated with each NCBI Taxonomy ID
    stored in the dataframe passed to the function.

    Use Entrez efetch function to pull down the NCBI accession number associated with
    each NCBI assembly ID. Return the dataframe with additional column titled
    'NCBI Assembly Numbers' containing all assembly numbers associated with each NCBI
    Taxonomy ID. The dataframe returned from get_accession_numbers() is stored in the
    variable 'species_table'.
    """

    with Entrez.elink(
        dbfrom="Taxonomy",
        id=taxonomy_id_column,
        db="Assembly",
        linkname="taxonomy_assembly",
    ) as handle_one:
        record_one = Entrez.read(handle_one)
        assembly_id_list = re.findall(
            "'\d*'", str(record_one[0]["LinkSetDb"][0]["Link"])
        )

    NCBI_accession_numbers_list = []
    for assembly_id in assembly_id_list:
        with Entrez.efetch(
            db="Assembly", id=assembly_id, rettype="docsum", retmode="xml"
        ) as handle_two:
            record_two = Entrez.read(handle_two, validate=False)
            NCBI_accession_number = record_two["DocumentSummarySet"]["DocumentSummary"][
                0
            ]["AssemblyAccession"]
            NCBI_accession_numbers_list.append(NCBI_accession_number)
    return str(NCBI_accession_numbers_list)


if __name__ == "__main__":
    main()
