#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Script pulls down assembles from NCBI database with genus/species names and taxonomy ID input

import argparse
import datetime
import logging
import re
from pathlib import Path

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

    # Capture date and time script is executed
    date = datetime.datetime.now()
    date_of_pulldown = date.strftime("%Y-%m-%d")
    time_of_pulldown = date.strftime("%H:%M")

    # Create parser object
    parser = argparse.ArgumentParser(
        prog="Extract_genomes_NCBI.py",
        description="Programme to pull down genomic assembleis from NCBI",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
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
    # If not given use standard input file ("Extract_genomes_NCBI_input_file.txt")
    parser.add_argument(
        "-i",
        "--input_file",
        type=Path,
        metavar="input file name",
        default=Path(__file__)
        .resolve()
        .parent.joinpath("Extract_genomes_NCBI_input_file.txt"),
        help="input filename",
    )
    # Add output file name option
    # If not given, file will be written to CWD
    # Development note - need to add file exte
    # nsion for default output file
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        metavar="output file name",
        default=Path.cwd().joinpath(
            "Extract_genomes_NCBI_{}_{}".format(date_of_pulldown, time_of_pulldown)
        ),
        help="output filename",
    )
    # Add log file name option
    # If not given not log file will be written out
    parser.add_argument(
        "-l",
        "--log",
        type=Path,
        metavar="log file name",
        default=None,
        help="Defines log file name and/or path",
    )

    # Parse arguments into args variable
    args = parser.parse_args()

    # Add users email address from parser
    Entrez.email = args.user_email

    # Initiate logger
    # Note: log file only created if specified at cmdline
    build_logger("Extract_genomes_NCBI", args.log, date_of_pulldown, time_of_pulldown)
    logger = logging.getLogger("Extract_genomes_NCBI")
    logger.info("Run initated")

    # Invoke main usage of script
    try:
        # Create dataframe storing genus, species and NCBI Taxonomy ID, called 'species_table'
        species_table = parse_input_file(args.input_file, logger)
    except Exception:
        logger.error(
            "::ERROR:: Error encounted during name and tax ID retrieval, see stack info:",
            exc_info=1,
            stack_info=1,
        )

    try:
        # pull down all accession numbers associated with each Taxonomy ID, from NCBI
        # add accession numbers to 'species_table' dataframe
        species_table = collate_accession_numbers(species_table, logger)
        logger.info("Generated species table")
        print("\nSpecies table:\n", species_table)
    except Exception:
        logger.error(
            "::ERROR:: Error encounted during accession number acquisition, see stack info::",
            exc_info=1,
            stack_info=1,
        )


def build_logger(
    script_name, log_file, date_of_pulldown, time_of_pulldown
) -> logging.Logger:
    """"Return a logger for this script.
    
    script_name: Name of script
    custom_string: Additional string parsed from cmdline by user - required for log
                    file to be written out
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
    if log_file is not None:
        file_log_handler = logging.FileHandler(log_file)
        file_log_handler.setLevel(logging.DEBUG)
        file_log_handler.setFormatter(log_formatter)
        logger.addHandler(file_log_handler)

    return logger


def parse_input_file(input_filename, logger):
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
    try:
        with open(input_filename) as file:
            input_list = file.read().splitlines()
    except Exception:
        logger.error("ERROR: Failed to read input file", exc_info=1, stack_info=1)

    number_of_lines = len(input_list)
    line_count = 1

    # Parse input, retrieving tax ID or scientific name as appropriate
    for line in input_list:
        logger.info("Processing line {} of {}".format(line_count, number_of_lines))
        if line[0] != "#":
            if line.startswith("NCBI:txid"):
                gs_name = get_genus_species_name(line[9:], logger)
                line_data = gs_name.split()
                line_data.append(line)
                all_species_data.append(line_data)
            else:
                tax_id = get_tax_ID(line, logger)
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


def get_genus_species_name(taxonomy_id, logger):
    """Fetch scientfic name associated with the NCBI Taxonomy ID.

    Use Entrez efetch function to pull down the scientific name (genus/species name)
    in the     NCBI Taxonomy database, associated with the taxonomy ID passed to the
    function.
    """

    logger.info("(Retrieving genus/species name using Taxonomy ID)")

    try:
        with Entrez.efetch(db="Taxonomy", id=taxonomy_id, retmode="xml") as handle:
            record = Entrez.read(handle)
    except Exception:
        logger.error("ERROR: Entrez failed to retrieve scientific name", exc_info=1)

    return record[0]["ScientificName"]


def get_tax_ID(genus_species, logger):
    """Pull down taxonomy ID from NCBI, using genus/species name as query.

    Use Entrez esearch function to pull down the NCBI Taxonomy ID of the
    species name passed to the function. Return the NCBI Taxonomy ID with
    the prefix 'NCBI:txid'.
    """

    logger.info("(Retrieving taxonomy ID using genus/species name)")
    try:
        with Entrez.esearch(db="Taxonomy", term=genus_species) as handle:
            record = Entrez.read(handle)
    except Exception:
        logger.error("ERROR: Entrez failed to retrieve taxonomy ID", exc_info=1)

    return "NCBI:txid" + record["IdList"][0]


def collate_accession_numbers(species_table, logger):
    """Return dataframe with column containing all associated NCBI accession numbers.

    Pass each Taxonomy ID in the 'NCBI Taxonomy ID' column of the 'species_table'
    dataframe to get_accession_numbers() to find all associated NCBI accession numbers.
    Format object returned from get_accession_numbers into a human readable list. Append list
    of accession numbers to 'all_accession_numbers' list. Create fourth column, called
    'NCBI Accession Numbers' in the 'species_table' dataframe and populate with data in
    'all_accession_numbers'. 
    
    Return modified dataframe, with four columns.
    """

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
        all_accession_numbers.append(get_accession_numbers(working_tax_id, logger))
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


def get_accession_numbers(taxonomy_id_column, logger):
    """Return all NCBI accession numbers associated with NCBI Taxonomy Id.

    Use Entrez elink function to pull down the assembly IDs of all genomic assemblies
    stored in the NCBI Assembly database that are associated with each NCBI Taxonomy ID
    stored in the dataframe passed to the function.
    Use Entrez epost to post all assembly IDs to NCBI as a single query for subsequent
    Entrez efetch of all associated accession numbers.
    Accession numbers are returned as a string 'NCBI_accession_numbers'.
    """

    try:
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
    except Exception:
        logger.error("ERROR: Entrez failed to retrieve assembly ID", exc_info=1)
    finally:
        logger.info("(Finished processing retrieval of assembly ID)", exc_info=1)

    try:
        epost_search_results = Entrez.read(
            Entrez.epost("Assembly", id=str(",".join(assembly_id_list)))
        )
        epost_webenv = epost_search_results["WebEnv"]
        epost_query_key = epost_search_results["QueryKey"]
    except Exception:
        logger.error("ERROR: Entrez failed to post assembly IDs", exc_info=1)
    finally:
        logger.info(
            "(Finihsed processing positing of assembly IDs for accession number fetch"
        )

    NCBI_accession_numbers_list = []

    try:
        with Entrez.efetch(
            db="Assembly",
            query_key=epost_query_key,
            WebEnv=epost_webenv,
            rettype="docsum",
            retmode="xml",
        ) as accession_handle:
            accession_record = Entrez.read(accession_handle, validate=False)
    except Exception:
        logger.error("ERROR: Entrez failed to retrieve accession numbers", exc_info=1)

    # Find total number of accession numbers
    number_of_accession_numbers = len(
        accession_record["DocumentSummarySet"]["DocumentSummary"]
    )

    # Retrieve accession numbers
    try:
        for index_number in range(number_of_accession_numbers):
            new_accession_number = accession_record["DocumentSummarySet"][
                "DocumentSummary"
            ][index_number]["AssemblyAccession"]
            NCBI_accession_numbers_list.append(new_accession_number)
            index_number += 1
    except Exception:
        logger.error(
            "ERROR: Error encounted when fetching accession numbers", exc_info=1
        )

    logger.info(
        "(Finished processing retrieval associated accession numbers of Taxonomy ID)"
    )

    # Process accession numbers into human readable list for dataframe
    NCBI_accession_numbers = ", ".join(NCBI_accession_numbers_list)

    return NCBI_accession_numbers


if __name__ == "__main__":
    main()
