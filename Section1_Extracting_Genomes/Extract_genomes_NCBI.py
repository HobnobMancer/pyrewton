#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Script pulls down assembles from NCBI database with genus/species names and taxonomy ID input

import argparse
import datetime
import logging
import re

from Bio import Entrez
import pandas as pd


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

    # Create parser object
    parser = argparse.ArgumentParser(
        prog="Extract_genomes_NCBI.py",
        description="Programme to pull down genomic assembleis from NCBI",
    )

    # Add arguments to parser
    # Add user email address
    parser.add_argument(
        "user_email",
        type=str,
        metavar="user email address",
        default=None,
        help="Email address of user, this must be provided",
    )
    # Add input file name option
    parser.add_argument(
        "-i",
        "--input_file",
        type=str,
        metavar="input file name",
        default="Extract_genomes_NCBI_input_file.txt",
        help="Name of input file (including extension) containing list of species, if not using input file provided (Default: Extract_genomes_NCBI_input_file.txt)",
    )
    # Add output file name option
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        metavar="output file name",
        default=None,
        help="Name of output file (including extension) if not using input file provided (Default: None)",
    )
    # Add log file name option
    parser.add_argument(
        "-l",
        "--log",
        type=str,
        metavar="log file name",
        default=None,
        help="Additional string added to log file name (Default: None)",
    )

    # Parse arguments into args variable
    args = parser.parse_args()

    # Add users email address from parser
    Entrez.email = args.user_email

    # Capture date and time script is executed
    date = datetime.datetime.now()
    date_of_pulldown = date.strftime("%Y-%m-%d")
    time_of_pulldown = date.strftime("%H:%M")

    # Initiate logger
    build_logger("Extract_genomes_NCBI", args.log, date_of_pulldown, time_of_pulldown)
    logger = logging.getLogger("Extract_genomes_NCBI")
    logger.info("Run initated")

    # create dataframe storing genus, species and NCBI Taxonomy ID, called 'species_table'
    # default: "Extract_genomes_NCBI_input_file.txt"
    species_table = parse_input_file(args.input_file)

    # pull down all accession numbers associated with each Taxonomy ID, from NCBI
    # add accession numbers to 'species_table' dataframe
    species_table = collate_accession_numbers(species_table)
    logger.info("Generated species table")
    print("\nSpecies table:\n", species_table)


def build_logger(
    script_name, custom_string, date_of_pulldown, time_of_pulldown
) -> logging.Logger:
    """"Return a logger for this script.
    
    script_name: Name of script
    custom_string: Additional string parsed from cmdline by user
    date_of_pulldown: Data run was initated
    time_of_pulldown: Time run was initated

    Enables logger for script, sets parameters and creates new file to store log.
    """

    logger = logging.getLogger(script_name)
    logger.setLevel(logging.DEBUG)

    # Set format of loglines
    log_formatter = logging.Formatter(
        script_name + ": {} - {}".format("%(asctime)s", "%(message)s")
    )

    # Setup console handler to log to terminal
    console_log_handler = logging.StreamHandler()
    console_log_handler.setLevel(logging.DEBUG)
    console_log_handler.setFormatter(log_formatter)
    logger.addHandler(console_log_handler)

    # Setup file handler to log to a file
    if custom_string is not None:
        file_log_handler = logging.FileHandler(
            "LOG:"
            + script_name
            + "__"
            + custom_string
            + "_DATE_{}_TIME_{}.log".format(date_of_pulldown, time_of_pulldown)
        )
    else:
        file_log_handler = logging.FileHandler(
            "LOG:"
            + script_name
            + "_DATE_{}_TIME_{}.log".format(date_of_pulldown, time_of_pulldown)
        )

    file_log_handler.setLevel(logging.DEBUG)
    file_log_handler.setFormatter(log_formatter)
    logger.addHandler(file_log_handler)

    return logger


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
    logger = logging.getLogger("Extract_genomes_NCBI")

    # open working_species_list.txt and extract lines, without new line character
    # then create genus name, species name and taxonomy ID tuplet
    all_species_data = []
    with open(input_filename) as file:
        input_list = file.read().splitlines()
        number_of_lines = len(input_list)
        line_count = 1
        logger.info(
            "Reading input file, and acquiring Tax IDs and Genus/Species names."
        )
        for line in input_list:
            logger.info("Processing line {} of {}".format(line_count, number_of_lines))
            if line[0] != "#":
                if line.startswith("NCBI:txid"):
                    gs_name = get_genus_species_name(line[9:])
                    line_data = gs_name.split()
                    line_data.append(line)
                    all_species_data.append(line_data)
                else:
                    tax_id = get_tax_ID(line)
                    line_data = line.split()
                    line_data.append(tax_id)
                    all_species_data.append(line_data)
            logger.info(
                "Finished processing line {} of {}".format(line_count, number_of_lines)
            )
            line_count += 1
    logger.info("Finished reading and closed input file")

    # create dataframe containg three columns: 'Genus', 'Species', 'NCBI Taxonomy ID'
    logger.info("Generating genus, species and taxonomy ID dataframe")
    species_table = pd.DataFrame(
        all_species_data, columns=["Genus", "Species", "NCBI Taxonomy ID"]
    )
    logger.info("Dataframe completed")
    return species_table


def get_genus_species_name(taxonomy_id):
    """Fetch scientfic name associated with the NCBI Taxonomy ID.

    Use Entrez efetch function to pull down the scientific name (genus/species name)
    in the     NCBI Taxonomy database, associated with the taxonomy ID passed to the
    function.
    """
    logger = logging.getLogger("Extract_genomes_NCBI")
    logger.info("(Successfully retrieved genus/species name using Taxonomy ID)")

    with Entrez.efetch(db="Taxonomy", id=taxonomy_id, retmode="xml") as handle:
        record = Entrez.read(handle)

    return record[0]["ScientificName"]


def get_tax_ID(genus_species):
    """Pull down taxonomy ID from NCBI, using genus/species name as query.

    Use Entrez esearch function to pull down the NCBI Taxonomy ID of the
    species name passed to the function. Return the NCBI Taxonomy ID with
    the prefix 'NCBI:txid'.
    """
    logger = logging.getLogger("Extract_genomes_NCBI")
    logger.info("(Retrieving taxonomy ID using genus/species name)")

    with Entrez.esearch(db="Taxonomy", term=genus_species) as handle:
        record = Entrez.read(handle)

    return "NCBI:txid" + record["IdList"][0]


def collate_accession_numbers(species_table):
    """Return dataframe with column containing all associated NCBI accession numbers.

    Pass each Taxonomy ID in the 'NCBI Taxonomy ID' column of the 'species_table'
    dataframe to get_accession_numbers() to find all associated NCBI accession numbers.
    Format object returned from get_accession_numbers into a human readable list. Append list
    of accession numbers to 'all_accession_numbers' list. Create fourth column, called
    'NCBI Accession Numbers' in the 'species_table' dataframe and populate with data in
    'all_accession_numbers'. 
    
    Return modified dataframe, with four columns.
    """
    logger = logging.getLogger("Extract_genomes_NCBI")

    all_accession_numbers = []
    tax_id_total_count = len(species_table["NCBI Taxonomy ID"])
    tax_id_counter = 1
    logger.info(
        "Aquiring accession numbers from NCBI Assembly database for NCBI Taxonomy IDs"
    )
    for NCBI_taxonomy_id in species_table["NCBI Taxonomy ID"]:
        logger.info(
            "Acquiring accession numbers for NCBI Taxonomy ID {} of {}".format(
                tax_id_counter, tax_id_total_count
            )
        )
        working_tax_id = NCBI_taxonomy_id[9:]
        all_accession_numbers.append(get_accession_numbers(working_tax_id))
        logger.info(
            "Completed retrieving accession numbers of Taxonomy ID {} of {}".format(
                tax_id_counter, tax_id_total_count
            )
        )
        tax_id_counter += 1

    # add accession numbers to the dataframe
    logger.info("Adding accession numbers to dataframe")
    species_table["NCBI Accession Numbers"] = all_accession_numbers
    return species_table


def get_accession_numbers(taxonomy_id_column):
    """Return all NCBI accession numbers associated with NCBI Taxonomy Id.

    Use Entrez elink function to pull down the assembly IDs of all genomic assemblies
    stored in the NCBI Assembly database that are associated with each NCBI Taxonomy ID
    stored in the dataframe passed to the function.
    Use Entrez epost to post all assembly IDs to NCBI as a single query for subsequent
    Entrez efetch of all associated accession numbers.
    Accession numbers are returned as a string 'NCBI_accession_numbers'.
    """
    logger = logging.getLogger("Extract_genomes_NCBI")

    with Entrez.elink(
        dbfrom="Taxonomy",
        id=taxonomy_id_column,
        db="Assembly",
        linkname="taxonomy_assembly",
    ) as assembly_number_handle:
        assembly_number_record = Entrez.read(assembly_number_handle)
        assembly_id_list = [
            dict["Id"] for dict in assembly_number_record[0]["LinkSetDb"][0]["Link"]
        ]

    logger.info(
        "(Finished processing retrieval of associated assembly IDs for Taxonomy ID)"
    )

    epost_search_results = Entrez.read(
        Entrez.epost("Assembly", id=str(",".join(assembly_id_list)))
    )
    epost_webenv = epost_search_results["WebEnv"]
    epost_query_key = epost_search_results["QueryKey"]
    logger.info(
        "(Finished processing posting of assembly IDs for accession number fetch for Taxonomy ID)"
    )

    NCBI_accession_numbers_list = []

    with Entrez.efetch(
        db="Assembly",
        query_key=epost_query_key,
        WebEnv=epost_webenv,
        rettype="docsum",
        retmode="xml",
    ) as accession_handle:
        accession_record = Entrez.read(accession_handle, validate=False)
        # Find total number of accession numbers
        number_of_accession_numbers = len(
            accession_record["DocumentSummarySet"]["DocumentSummary"]
        )
        # Add each accession number in the accession record to the accession number list
        for index_number in range(number_of_accession_numbers):
            new_accession_number = accession_record["DocumentSummarySet"][
                "DocumentSummary"
            ][index_number]["AssemblyAccession"]
            NCBI_accession_numbers_list.append(new_accession_number)
            index_number += 1

    logger.info(
        "(Finished processing retrieval associated accession numbers of Taxonomy ID)"
    )

    # Process accession numbers into human readable list for dataframe
    NCBI_accession_numbers = ", ".join(NCBI_accession_numbers_list)

    return NCBI_accession_numbers


if __name__ == "__main__":
    main()
