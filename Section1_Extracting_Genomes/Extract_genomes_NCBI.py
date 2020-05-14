#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Pull down genomic assemblies from NCBI database.

:func main: generate a dataframe of scientific names, taxonomy IDs and accession numbers
:func build_logger: creates logger object
:func parse_input_file: parse input file
:func get_genus_species_name: retrieve scientific name from taxonomy ID
:func get_tax_id: retrieve NCBI taxonomy ID from scientific name
:func collate_accession_numbers: parse taxonomy ID column in dataframe
:func get_accession_numbers: retrieves all accessions associated to given taxonomy ID

Generates dataframe containing scientific names, taxonomy IDs and accession numbers.
Pulls down and stores genomic assemblies from NCBI Assebmly database.
"""
# Script pulls down assembles from NCBI database with genus/species names and taxonomy ID input

import argparse
import datetime
import logging
import re
import sys
import time
from pathlib import Path

from Bio import Entrez
import pandas as pd


def main():
    """Generate datafame containing scientific names, taxonomy IDs and accession numbers.
    
    Pass input file (containing unique species on each line, idenfitied by their genus/species
    name or NCBI Taxonomy ID) tp parse_input_file() function, which acquires missing genus/species
    names and NCBI taxonomy IDs as appropriate. Genus/species names and associated NCBI
    Taxonomy ID are stored in a generated dataframe with three columns: 'Genus', 'Species',
    and 'NCBI Taxonomy ID'. Parse_input_file() returns the generated dataframe which is stored in
    the variable 'species_table'.
    Pass the 'species_table' dataframe to collate_accession_numbers(), which acquires all associated
    NCBI accession numbers for each Taxonomy ID. Acquired accession numbers are stored in a new
    column in the dataframe, titled 'NCBI Accession Numbers'. The modified dataframe is returned
    from collate_accession_numbers() and stored in the variable 'species_table'.

    Return 'species_table' dataframe.
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
    # If not given input will be taken from STDIN
    parser.add_argument(
        "-i",
        "--input_file",
        type=Path,
        metavar="input file name",
        default=sys.stdin,
        help="input filename",
    )
    # Add output file name option
    # Must include file extension
    # If not given, output will be written to STDOUT
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        metavar="output file name",
        default=sys.stdout,
        help="output filename",
    )
    # Add log file name option
    # If not given, no log file will be written out
    parser.add_argument(
        "-l",
        "--log",
        type=Path,
        metavar="log file name",
        default=None,
        help="Defines log file name and/or path",
    )
    # Add custom maximum number of retries if network error is encountered
    parser.add_argument(
        "-r",
        "--retries",
        type=int,
        metavar="maximum number of retries",
        default=10,
        help="Defines the maximum number of retries if network errors are encountered",
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
    # Create dataframe storing genus, species and NCBI Taxonomy ID, called 'species_table'
    species_table = parse_input_file(args.input_file, logger, args.retries)

    species_table = collate_accession_numbers(species_table, logger, args.retries)
    logger.info("Generated species table")
    print("\nSpecies table:\n", species_table)


def build_logger(
    script_name, log_file, date_of_pulldown, time_of_pulldown
) -> logging.Logger:
    """Return a logger for this script.

    Enables logger for script, sets parameters and creates new file to store log.
    
    :param script_name: Name of script
    :param custom_string: Additional string parsed from cmdline by user - required for log
                    file to be written out
    :param date_of_pulldown: Data run was initiated
    :param time_of_pulldown: Time run was initiated

    Return logger object.
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


def parse_input_file(input_filename, logger, retries):
    """Parse input file, returning dataframe of species names and NCBI Taxonomy IDs.

    Read input file passed to the function, performing the appropriate action depending
    on each line's content.
    Comments: Indicated with the first character of '#', perform no action.
    NCBI Taxonomy ID: Indicated with the first nine characters of 'NCBI:txid', pass
    the Taxonomy ID (excluding the 'NCBI:txid' prefix) to get_genus_species_name() to
    return associated scientific name.
    Genus/species name: Indicated by lack of '#' and 'NCBI:txid', pass genus/species
    name to get get_tax_id() to return associated NCBI Taxonomy ID with 'NCBI:txid'
    prefix.
    Split genus and species name into separate variables. Store genus, species and
    associated Taxonomy ID in the list 'line_data'. Append 'line_data' to
    'all_species_data' list. Repeat for each line.
    Generate a dataframe with three columns: 'Genus', 'Species' and 'NCBI Taxonomy ID'.
    Use data from 'all_species_data' to fill out dataframe.
    
    :param input_filename: args, if specific name of input file, otherwise input taken from STDIN
    :param logger: logger object
    :param retries: parser argument, maximum number of retries excepted if network error encountered

    Return dataframe.
    """
    # open working_species_list.txt and extract lines, without newline character
    # then create genus name, species name and taxonomy ID tuplet
    all_species_data = []

    # test path to input file exists, if not exit programme
    if input_filename.is_file() is False:
        # report to user and exit programme
        logger.error(
            (
                "Input file not found. Check filename, extension and directory is correct."
                "Terminating program."
            ),
            exc_info=1,
        )
        sys.exit(1)

    # if path to input file exists proceed
    # parse input file
    try:
        with open(input_filename) as file:
            input_list = file.read().splitlines()
    except IOError:
        # if input file cannot be read, programme terminates
        logger.error(
            "Input file was found but could not be read\nTerminating programming",
            exc_info=1,
        )
        sys.exit(1)

    number_of_lines = len(input_list)
    line_count = 1

    # Parse input, retrieving tax ID or scientific name as appropriate
    for line in input_list:
        logger.info("Processing line {} of {}".format(line_count, number_of_lines))
        if line[0] != "#":
            if line.startswith("NCBI:txid"):
                gs_name = get_genus_species_name(line[9:], logger, line_count, retries)
                line_data = gs_name.split()
                line_data.append(line)
                all_species_data.append(line_data)
            else:
                tax_id = get_tax_id(line, logger, line_count, retries)
                line_data = line.split()
                line_data.append(tax_id)
                all_species_data.append(line_data)
        logger.info(
            "Finished processing line {} of {}".format(line_count, number_of_lines)
        )
        line_count += 1

    logger.info("Finished reading and closed input file")

    # create dataframe containing three columns: 'Genus', 'Species', 'NCBI Taxonomy ID'
    logger.info("Generating genus, species and taxonomy ID dataframe")
    species_table = pd.DataFrame(
        all_species_data, columns=["Genus", "Species", "NCBI Taxonomy ID"]
    )
    logger.info("Dataframe completed")
    return species_table


def get_genus_species_name(taxonomy_id, logger, line_number, retries):
    """Fetch scientific name associated with the NCBI Taxonomy ID.

    Use Entrez efetch function to pull down the scientific name (genus/species name)
    in the NCBI Taxonomy database, associated with the taxonomy ID passed to the
    function.
    
    :param taxonomy_id: str, NCBI taxonomy ID
    :param logger: logger object
    :param line_number: int, line number in input file containing taxonomy ID
    :param retries: parser argument, maximum number of retries excepted if network error encountered
    
    Return scientific name.
    """
    logger.info("(Retrieving scientific name for NCBI:txid{})".format(taxonomy_id))

    with entrez_retry(
        logger, retries, Entrez.efetch, db="Taxonomy", id=taxonomy_id, retmode="xml"
    ) as handle:
        record = Entrez.read(handle)

    # extract scientific name from record
    try:
        return record[0]["ScientificName"]

    except IndexError:
        logger.error(
            (
                "Entrez failed to retrieve scientific name, for species in line {} of input file."
                "Potential typo in taxonomy ID, check input."
            ).format(line_number),
            exc_info=1,
        )
        return "NA"


def get_tax_id(genus_species, logger, line_number, retries):
    """Pull down taxonomy ID from NCBI, using genus/species name as query.

    Use Entrez esearch function to pull down the NCBI Taxonomy ID of the
    species name passed to the function. Return the NCBI Taxonomy ID with
    the prefix 'NCBI:txid'.
    
    :param genus_species: str, scientific name of species
    :param logger: logger object
    :param line_number: int, number of line containing the species name in the input file.
    :param retries: parser argument, maximum number of retries excepted if network error encountered
    
    Return NCBI taxonomy ID.
    """
    # check for potential mistake in taxonomy ID prefix
    if bool(re.search(r"\d", genus_species)) is True:
        logger.error(
            (
                "Warning: Number found in genus-species name, when trying to retrieve taxonomy ID"
                "for species on line {}. Potential typo in taxonomy ID, so taxonomy ID "
                "misinterpreted as genus-species name."
            ).format(line_number),
            exc_info=1,
        )
        sys.exit(1)

    else:
        logger.info("(Retrieving NCBI taxonomy ID for {})".format(genus_species))

        with entrez_retry(
            logger, retries, Entrez.esearch, db="Taxonomy", term=genus_species
        ) as handle:
            record = Entrez.read(handle)

    # extract taxonomy ID from record
    try:
        return "NCBI:txid" + record["IdList"][0]

    except IndexError:
        logger.error(
            (
                "Entrez failed to retrieve taxonomy ID, for species in line {} of input file."
                "Potential typo in species name."
            ).format(line_number),
            exc_info=1,
        )
        return "NA"


def collate_accession_numbers(species_table, logger, retries):
    """Return dataframe with column containing all associated NCBI accession numbers.

    Pass each Taxonomy ID in the 'NCBI Taxonomy ID' column of the 'species_table'
    dataframe to get_accession_numbers() to find all associated NCBI accession numbers.
    Format object returned from get_accession_numbers into a human-readable list. Append list
    of accession numbers to 'all_accession_numbers' list. Create fourth column, called
    'NCBI Accession Numbers' in the 'species_table' dataframe and populate with data in
    'all_accession_numbers'. 
    
    :param species_table: dataframe, dataframe containing scientific names and taxonomy IDs
    :param logger: logger object
    :param retries: parser argument, maximum number of retries excepted if network error encountered
    
    Return modified dataframe, with four columns.
    """
    logger.info(
        "Aquiring accession numbers from NCBI Assembly database for NCBI Taxonomy IDs"
    )
    species_table["NCBI Accession Numbers"] = species_table.apply(
        lambda column: get_accession_numbers(
            column["NCBI Taxonomy ID"][9:], logger, retries
        ),
        axis=1,
    )

    logger.info("Finished adding accession numbers to dataframe")
    return species_table


def get_accession_numbers(taxonomy_id, logger, retries):
    """Return all NCBI accession numbers associated with NCBI Taxonomy Id.

    Use Entrez elink function to pull down the assembly IDs of all genomic
    assemblies stored in the NCBI Assembly database that are associated
    with each NCBI Taxonomy ID stored in the dataframe passed to the function.
    Use Entrez epost to post all assembly IDs to NCBI as a single query for
    subsequent Entrez efetch of all associated accession numbers.
    Accession numbers are returned as a string 'NCBI_accession_numbers'.
    
    :param taxonomy_id: str, NCBI taxonomy ID
    :param logger: logger object
    :param retries: parser argument, maximum number of retries excepted if network error encountered

    
    Return NCBI accession numbers.
    """
    # If previously failed to retrieve the taxonomy ID cancel retrieval of accession numbers
    if taxonomy_id == "NA":
        return "NA"

    # Retrieve all IDs of genomic assemblies for taxonomy ID

    logger.info("(Retrieving assembly IDs for NCBI:txid{}".format(taxonomy_id))

    with entrez_retry(
        logger,
        retries,
        Entrez.elink,
        dbfrom="Taxonomy",
        id=taxonomy_id,
        db="Assembly",
        linkname="taxonomy_assembly",
    ) as assembly_number_handle:
        assembly_number_record = Entrez.read(assembly_number_handle)

    # test record was returned, if failed to return exit retrieval of assembly IDs
    if assembly_number_record == "NA":
        logger.error(
            (
                "Entrez failed to retrieve assembly IDs, for NCBI:txid{}."
                "Exiting retrieval of accession numbers for {}"
            ).format(taxonomy_id, taxonomy_id),
            exc_info=1,
        )
        return assembly_number_record

    # extract assembly IDs from record
    try:
        assembly_id_list = [
            dict["Id"] for dict in assembly_number_record[0]["LinkSetDb"][0]["Link"]
        ]

    except IndexError:
        logger.error(
            (
                "Entrez failed to retrieve assembly IDs, for NCBI:txid{}."
                "Exiting retrieval of accession numbers for {}"
            ).format(taxonomy_id, taxonomy_id),
            exc_info=1,
        )
        return "NA"

    logger.info("(Finished processing retrieval of assembly ID)")

    # Post all assembly IDs to Entrez-NCBI for downstream pulldown of accession number
    # assocated with each assembly ID

    # compile list of ids in suitable format for epost
    id_post_list = str(",".join(assembly_id_list))
    epost_search_results = Entrez.read(
        entrez_retry(logger, retries, Entrez.epost, "Assembly", id=id_post_list)
    )

    # test record was returned, if failed to return exit retrieval of assembly IDs
    if epost_search_results == "NA":
        logger.error(
            (
                "Entrez failed to post assembly IDs, for NCBI:txid{}."
                "Exiting retrieval of accession numbers for {}"
            ).format(taxonomy_id, taxonomy_id),
            exc_info=1,
        )
        return epost_search_results

    # Retrieve web environment and query key from Entrez epost
    epost_webenv = epost_search_results["WebEnv"]
    epost_query_key = epost_search_results["QueryKey"]

    logger.info(
        "(Finished processing positing of assembly IDs for accession number fetch)"
    )

    # Pull down all accession numbers

    ncbi_accession_numbers_list = []

    with entrez_retry(
        logger,
        retries,
        Entrez.efetch,
        db="Assembly",
        query_key=epost_query_key,
        WebEnv=epost_webenv,
        rettype="docsum",
        retmode="xml",
    ) as accession_handle:
        accession_record = Entrez.read(accession_handle, validate=False)

    if accession_record == "NA":
        logger.error(
            (
                "Entrez failed to post assembly IDs, for NCBI:txid{}."
                "Exiting retrieval of accession numbers for {}"
            ).format(taxonomy_id, taxonomy_id),
            exc_info=1,
        )
        return accession_record

    # Extract accession numbers from document summary

    try:
        # Find total number of accession numbers
        number_of_accession_numbers = len(
            accession_record["DocumentSummarySet"]["DocumentSummary"]
        )

        for index_number in range(number_of_accession_numbers):
            new_accession_number = accession_record["DocumentSummarySet"][
                "DocumentSummary"
            ][index_number]["AssemblyAccession"]
            ncbi_accession_numbers_list.append(new_accession_number)
            index_number += 1

    except IndexError:
        logger.error(
            "No accession number retrieved from NCBI\nfor assembly ID {} of {}".format(
                index_number, number_of_accession_numbers
            ),
            exc_info=1,
        )
        return "NA"

    logger.info(
        "(Finished processing retrieval associated accession numbers of Taxonomy ID)"
    )

    # Process accession numbers into human readable list for dataframe
    ncbi_accession_numbers = ", ".join(ncbi_accession_numbers_list)

    return ncbi_accession_numbers


def entrez_retry(logger, retries, entrez_func, *func_args, **func_kwargs):
    """Calls to NCBI using Entrez.

    Maximum number of retries is 10, retry initated when network error encountered.

    :param logger: logger object
    :param retries: parser argument, maximum number of retries excepted if network error encountered
    :param entrez_func: function, call method to NCBI
    :param *func_args: tuple, arguments passed to Entrez function
    :param ** func_kwargs: dictionary, keyword arguments passed to Entrez function

    Returns record.
    """
    record = None
    retries = 10
    tries = 0

    while record is None and tries < retries:
        try:
            record = entrez_func(*func_args, **func_kwargs)

        except IOError:
            # log retry attempt
            if tries < retries:
                logger.warning(
                    "Network error encountered during try no.{}.\nRetrying in 10s".format(
                        tries
                    ),
                    exc_info=1,
                )
                time.sleep(10)
            tries += 1

    if record is None:
        logger.error(
            "Network error encountered too many times. Exiting attempt to call to NCBI"
        )
        return "NA"

    return record


if __name__ == "__main__":
    main()
