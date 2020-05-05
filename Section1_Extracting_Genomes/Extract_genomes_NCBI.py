#!/usr/bin/env python
# -*- coding: utf-8 -*-

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
    species_table = parse_input_file(args.input_file, logger)

    species_table = collate_accession_numbers(species_table, logger)
    logger.info("Generated species table")
    print("\nSpecies table:\n", species_table)


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

    # test path to input file exists, if not exit programme
    if input_filename.is_file() == False:
        # report to user and exit programme
        logger.error(
            "Input file not found.\nCheck filename, extension and directory is correct.\nTerminating program.",
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
                gs_name = get_genus_species_name(line[9:], logger, line_count)
                line_data = gs_name.split()
                line_data.append(line)
                all_species_data.append(line_data)
            else:
                tax_id = get_tax_ID(line, logger, line_count)
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


def get_genus_species_name(taxonomy_id, logger, line_number):
    """Fetch scientfic name associated with the NCBI Taxonomy ID.

    Use Entrez efetch function to pull down the scientific name (genus/species name)
    in the     NCBI Taxonomy database, associated with the taxonomy ID passed to the
    function.
    """

    logger.info("(Retrieving genus/species name using Taxonomy ID)")

    try:
        with Entrez.efetch(db="Taxonomy", id=taxonomy_id, retmode="xml") as handle:
            record = Entrez.read(handle)

    except IndexError:
        logger.error(
            "Entrez failed to retrieve scientific name, for speices in line {} of input file.\nPotential tpyo in taxonomy ID, check input.\nTerminating programming to avoid downstream processing errors".format(
                line_number
            ),
            exc_info=1,
        )
        sys.exit(1)

    except IOError:
        # log error
        logger.error(
            "Network error encountered during retrieval of scientific\nname for species in line {} of inputfile,\nretrying in 10s".format(
                line_number
            )
        )
        time.sleep(10)
        # initiate retries
        record = get_g_s_name_retry(taxonomy_id, logger, line_number)

    return record[0]["ScientificName"]


def get_tax_ID(genus_species, logger, line_number):
    """Pull down taxonomy ID from NCBI, using genus/species name as query.

    Use Entrez esearch function to pull down the NCBI Taxonomy ID of the
    species name passed to the function. Return the NCBI Taxonomy ID with
    the prefix 'NCBI:txid'.
    """
    # check for potential mistake in taxonomy ID prefix
    if bool(re.search(r"\d", genus_species)) == True:
        logger.error(
            "Warning: Number found in genus-species name,\n when trying to retrieve taxonomy ID for species on line {}.\nPotential typo in taxonomy ID, so taxonomy ID misinturpretted as genus-species name.\nTerminating program to avoid downstream processing errors".format(
                line_number
            ),
            exc_info=1,
        )
        sys.exit(1)

    else:
        logger.info("(Retrieving taxonomy ID using genus/species name)")

        try:
            with Entrez.esearch(db="Taxonomy", term=genus_species) as handle:
                record = Entrez.read(handle)

        except IndexError:
            logger.error(
                "Entrez failed to retrieve taxonomy ID, for species in line {} of input file.\nPotential typo in species name.\nTerminating program to avoid downstream processing errors".format(
                    line_number
                ),
                exc_info=1,
            )
            sys.exit(1)

        except IOError:
            # log error
            logger.error("Network error ecountered, retrying in 10s")
            time.sleep(10)
            # initate retries
            record = get_t_id_retry(genus_species, logger, line_number)

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


def get_accession_numbers(taxonomy_id, logger):
    """Return all NCBI accession numbers associated with NCBI Taxonomy Id.

    Use Entrez elink function to pull down the assembly IDs of all genomic assemblies
    stored in the NCBI Assembly database that are associated with each NCBI Taxonomy ID
    stored in the dataframe passed to the function.
    Use Entrez epost to post all assembly IDs to NCBI as a single query for subsequent
    Entrez efetch of all associated accession numbers.
    Accession numbers are returned as a string 'NCBI_accession_numbers'.
    """

    # Retrieve all IDs of genomic assemblies for taxonomy ID

    try:
        with Entrez.elink(
            dbfrom="Taxonomy",
            id=taxonomy_id,
            db="Assembly",
            linkname="taxonomy_assembly",
        ) as assembly_number_handle:
            assembly_number_record = Entrez.read(assembly_number_handle)
            assembly_id_list = [
                dict["Id"] for dict in assembly_number_record[0]["LinkSetDb"][0]["Link"]
            ]

    except IndexError:
        logger.error(
            "Entrez failed to retrieve assembly IDs, for NCBI:txid{}.\nExiting retrieval of accession numbers for {}".format(
                taxonomy_id, taxonomy_id
            ),
            exc_info=1,
        )
        return ()

    except IOError:
        # log error
        logger.error(
            "Network error encountered during retrieval of assembly IDs, for NCBI:txid{},\n retrying in 10s".format(
                taxonomy_id
            )
        )
        time.sleep(10)
        # initate retries
        assembly_id_list = get_a_id_retry(taxonomy_id, logger)
        if assembly_id_list == None:
            return ()

    logger.info("(Finished processing retrieval of assembly ID)")

    # Post all assembly IDs to Entrez-NCBI for downstream pulldown of accession number
    # assocated with each assembly ID

    try:
        epost_search_results = Entrez.read(
            Entrez.epost("Assembly", id=str(",".join(assembly_id_list)))
        )

    except RuntimeError:
        logger.error(
            "Entrez failed to post assembly IDs for NCBI:txid{}.\nEntrez could not retrieve document summary.\nPotential incorrect formatting of assembly IDs,\nor query too large.\nExiting retrieval of accession numbers for NCBI:txid{}".format(
                taxonomy_id, taxonomy_id
            ),
            exc_info=1,
        )
        return ()

    except IOError:
        # log error
        logger.error(
            "Network error encountered during assembly id posting for NCBI:txid{},\n retrying in 10s".format(
                taxonomy_id
            )
        )
        time.sleep(10)
        # initiate retries
        epost_search_results = post_a_ids_retry(assembly_id_list, logger, taxonomy_id)
        if epost_search_results == None:
            return ()

    # Retrieve web environment and query key from Entrez epost
    epost_webenv = epost_search_results["WebEnv"]
    epost_query_key = epost_search_results["QueryKey"]

    logger.info(
        "(Finihsed processing positing of assembly IDs for accession number fetch)"
    )

    # Pull down all accession numbers

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

    except RuntimeError:
        logger.error(
            "Entrez could not retrieve document summary when retrieving accession numbers.\nPotential incorrect formatting of assembly IDs,\nor query too large.\nExiting retrieval of accession numbers for NCBI:txid{}".format(
                taxonomy_id
            ),
            exc_info=1,
        )
        return ()

    except IOError:
        # log error
        logger.error(
            "Network error encountered during accession number retrieval\nfor NCBI:txid{},retrying in 10s".format(
                taxonomy_id
            )
        )
        time.sleep(10)
        # initate retries
        accession_record = get_a_n_retry(
            epost_query_key, epost_webenv, logger, taxonomy_id
        )
        if accession_record == None:
            return ()

    # Find total number of accession numbers
    number_of_accession_numbers = len(
        accession_record["DocumentSummarySet"]["DocumentSummary"]
    )

    # Extract accession numbers from document summary

    try:
        for index_number in range(number_of_accession_numbers):
            new_accession_number = accession_record["DocumentSummarySet"][
                "DocumentSummary"
            ][index_number]["AssemblyAccession"]
            NCBI_accession_numbers_list.append(new_accession_number)
            index_number += 1

    except IndexError:
        logger.error(
            "No accession number retrieved from NCBI\nfor assembly ID {} of {}".format(
                index_number, number_of_accession_numbers
            ),
            exc_info=1,
        )

    logger.info(
        "(Finished processing retrieval associated accession numbers of Taxonomy ID)"
    )

    # Process accession numbers into human readable list for dataframe
    NCBI_accession_numbers = ", ".join(NCBI_accession_numbers_list)

    return NCBI_accession_numbers


# Functions for retrying call to NCBI with network error is encountered

# If network error encountered during retrieval of scientific name
def get_g_s_name_retry(taxonomy_id, logger, line_number, retries=10):
    """Retries call to NCBI Taxnomy database to retrieve scientific name.

    Maxinum number of retries is 10. Retry initated when network error
    encountered. If no recorded is returned for non-network error issue,
    such as record does not exist, or typo in given taxonomy ID,
    programme terminates.

    Function only invoked if network  error encountered during call to NCBI
    using Entrez in another function.

    Return record.
    """

    record = None
    tries = 0

    while record is None and tries < retries:
        try:
            with Entrez.efetch(db="Taxonomy", id=taxonomy_id, retmode="xml") as handle:
                record = Entrez.read(handle)

        # If network error not encountered but no record returned, terminate programme
        except IndexError:
            logger.error(
                "Entrez failed to retrieve scientific name, for speices in line {} of input file.\nPotential tpyo in taxonomy ID, check input.\nTerminating programming to avoid downstream processing errors".format(
                    line_number
                ),
                exc_info=1,
            )
            sys.exit(1)

        except IOError:
            # log retry attempt
            if tries < tries:
                logger.error(
                    "Network error encountered during try no.{}.\nRetrying in 10s".format(
                        tries
                    ),
                    exc_info=1,
                )
                time.sleep(10)
            tries += 1

    if record is None:
        logger.error(
            "Network error encountered after 10 attempts.\nTerminating programme",
            exc_info=1,
        )
        sys.exit(1)

    return record


# If network error encountered during retrieval of taxonomy ID
def get_t_id_retry(genus_species, logger, line_number, retries=10):
    """Retries call to NCBI Taxonomy database to retrieve taxonomy ID.

    Maxinum number of retries is 10. Retry initated when network error
    encountered. If no recorded is returned for non-network error issue,
    such as record does not exist, or type in the given scientific name,
    programme terminates.

    Function all invoked if network error encountered during call to NCBI
    using Entrez in another function.

    Return record.
    """

    record = None
    tries = 0

    while record is None and tries < retries:
        try:
            with Entrez.esearch(db="Taxonomy", term=genus_species) as handle:
                record = Entrez.read(handle)

        # If network error not encountered but no recorded returned, terminate programme
        except IndexError:
            logger.error(
                "Entrez failed to retrieve taxonomy ID, for species in line {} of input file.\nPotential typo in species name.\nTerminating program to avoid downstream processing errors".format(
                    line_number
                ),
                exc_info=1,
            )
            sys.exit(1)

        except IOError:
            # log retry attempt
            if tries < tries:
                logger.error(
                    "Network error encountered during try no.{}.\nRetrying in 10s".format(
                        tries
                    ),
                    exc_info=1,
                )
                time.sleep(10)
            tries += 1

    if record is None:
        logger.error(
            "Network error encountered after 10 attempts.\nTerminating programme",
            exc_info=1,
        )
        sys.exit(1)

    return record


# If network error encountered during retrieval of assembly IDs
def get_a_id_retry(taxonomy_id, logger, retries=10):
    """Retries call to NCBI Assembly database to retrieve accession numbers.

    Maxinum number of retries is 10. Retry initated when network error
    encountered. If no record is returned for non-network error issue,
    such as record does not exist, or type is given taxonomy ID,
    exits function.

    Function only invoked if network error encountered during call to
    NCBI using Entrez in another function.

    Return list of assembly IDs for given taxonomy ID.
    """

    record = None
    tries = 0

    while record is None and tries < retries:
        try:
            with Entrez.elink(
                dbfrom="Taxonomy",
                id=taxonomy_id,
                db="Assembly",
                linkname="taxonomy_assembly",
            ) as handle:
                record = Entrez.read(handle)
                assembly_id_list = [
                    dict["Id"] for dict in record[0]["LinkSetDb"][0]["Link"]
                ]

        # If network error not encountered but no record returned
        # terminate collection of accession numbers
        except IndexError:
            logger.error(
                "Entrez failed to retrieve assembly IDs for NCBI:txid{}.\nExiting retrieval of accession numbers for {}".format(
                    taxonomy_id, taxonomy_id, exc_info=1
                )
            )
            return ()

        except IOError:
            # log retry attempt
            if tries < tries:
                logger.error(
                    "Network error encountered during try no.{}.\nRetrying in 10s".format(
                        tries
                    ),
                    exc_info=1,
                )
                time.sleep(10)
            tries += 1

    if record is None:
        logger.error(
            "Network error encountered after 10 attempts.\nExiting retrieval of accession numbers for NCBI:txid{}".format(
                taxonomy_id
            ),
            exc_info=1,
        )
        return ()

    return assembly_id_list


# If network error encountered during posting of assembly IDs
def post_a_ids_retry(assembly_id_list, logger, taxonomy_id, retries=10):
    """Retries call to NCBI Assembly database to post assembly IDs.

    Maxinum number of retries is 10. Retry initated when network error
    encountered. If no record is returned for non-network error issue,
    such as post is too large, exits function.

    Function only invoked if network error encountered during call to
    NCBI using Entrez in another function.

    Return record.
    """

    record = None
    tries = 0

    while record is None and tries < retries:
        try:
            record = Entrez.read(
                Entrez.epost("Assembly", id=str(",".join(assembly_id_list)))
            )

        # If non-network work error encountered
        # terminate posting of assembly IDs
        except RuntimeError:
            logger.error(
                "Entrez failed to post assembly IDs for NCBI:txid{}.\nEntrez could not retrieve document summary.\nPotentially incorrect formatting of assembly IDS,\nor query too large.\nExiting retrieval of accession numbers for NCBI:txid{}".format(
                    taxonomy_id, taxonomy_id
                ),
                exc_info=1,
            )
            return ()

        except IOError:
            # log retry attempt
            if tries < tries:
                logger.error(
                    "Network error encountered during try no.{}.\nRetrying in 10s".format(
                        tries
                    ),
                    exc_info=1,
                )
                time.sleep(10)
            tries += 1

    if record is None:
        logger.error(
            "Network error encountered after 10 attempts.\nExiting retrieval of accession numbers for NCBI:txid{}".format(
                taxonomy_id
            ),
            exc_info=1,
        )
        return ()

    return record


# If network error encountered during retrieval of accession numbers
def get_a_n_retry(epost_query_key, epost_webenv, logger, taxonomy_id, retries=10):
    """Retries call to NCBI Assembly database to retrieve accession numbers.

    Maxinum number of retries is 10. Retry initated when network error
    encountered. If no record is returned for non-network error issue,
    exits retrieval of accession numbers for given taxonomy ID.

    Function only invoked if network error encountered during call to
    NCBI using Entrez in another function.

    Return record.
    """

    record = None
    tries = 0

    while record is None and tries < retries:
        try:
            with Entrez.efetch(
                db="Assembly",
                query_key=epost_query_key,
                WebEnv=epost_webenv,
                rettype="docsum",
                retmode="xml",
            ) as handle:
                record = Entrez.read(handle, validate=False)

        except RuntimeError:
            logger.error(
                "Entrez could not retrieve document summary when retrieving accession numbers.\nPotential incorrect formatting of assembly IDs,\nor query too large.\nExiting retrieval of accession numbers for NCBI:txid{}".format(
                    taxonomy_id
                ),
                exc_info=1,
            )
            return ()

        except IOError:
            # log retry attempt
            if tries < tries:
                logger.error(
                    "Network error encountered during try no.{}.\nRetrying in 10s".format(
                        tries
                    ),
                    exc_info=1,
                )
                time.sleep(10)
            tries += 1

    if record is None:
        logger.error(
            "Network error encountered after 10 attempts.\nExiting retrieval of accession numbers for NCBI:txid{}".format(
                taxonomy_id
            ),
            exc_info=1,
        )
        return ()

    return record


if __name__ == "__main__":
    main()
