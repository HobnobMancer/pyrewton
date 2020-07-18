#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author:
# Emma E. M. Hobbs

# Contact
# eemh1@st-andrews.ac.uk

# Emma E. M. Hobbs,
# Biomolecular Sciences Building,
# University of St Andrews,
# North Haugh Campus,
# St Andrews,
# KY16 9ST
# Scotland,
# UK

# The MIT License
"""Pull down GenBank files (.gbff) from NCBI database.

:cmd_args --user: email address of user
:cmd_args --dataframe: output directory for dataframe
:cmd_args --force: force writing in output directory
:cmd_args --genbank: enable/disable GenBank file download
:cmd_args --input_file: path of input file
:cmd_args --log: enable log file writing
:cmd_args --nodelete: enable/disabling removal of exiting files in output
:cmd_args --output: output directory for downloaded files
:cmd_args --retries: maximum number of retries if network error occurs
:cmd_args --timeout: timeout limit of URL connection
:cmd_args --verbose: set logging level to 'INFO'

:func build_parser: create parser object
:func main: generate a dataframe of scientific names, taxonomy IDs and accession numbers
:func build_logger: creates logger object
:func make_output_directory: create directory for genomic files to be written to
:func parse_input_file: parse input file
:func parse_line: coordinate retrieval of scientific names and taxonomy IDs
:func get_genus_species_name: retrieve scientific name from taxonomy ID
:func get_tax_id: retrieve NCBI taxonomy ID from scientific name
:func get_accession_numbers: retrieves all accessions associated to given taxonomy ID
:func get_genbank_files: organise download of genbank files
:func compile_URL: create URL for downloading file
:func download_file: download file using provided URL
:func entrez_retry: perform call to NCBI using Entrez
:func write_out_dataframe: write out species table as .csv file

Generates dataframe containing scientific names, taxonomy IDs and accession numbers.
Pulls down and stores genomic assemblies and GenBank filesfrom NCBI Assembly database.
"""

import logging
import re
import sys
import time

from socket import timeout
from typing import List, Optional
from urllib.error import HTTPError, URLError
from urllib.request import urlopen

import pandas as pd

from Bio import Entrez
from tqdm import tqdm

from pyrewton.loggers import build_logger
from pyrewton.parsers.parser_get_ncbi_genomes import build_parser
from pyrewton.file_io import make_output_directory, write_out_dataframe


def main(argv: Optional[List[str]] = None, logger: Optional[logging.Logger] = None):
    """Generate dataframe containing scientific names, taxonomy IDs and accession numbers.

    Parse input file to create dataframe ('species_table') of scientific name and taxonomy ID.

    Pass the 'species_table' dataframe to collate_accession_numbers(), to retrieve all associated
    NCBI accession numbers. Store accession numbers in new column.

    If enabled, write out dataframe to .csv file.

    Return 'species_table' dataframe.
    """
    # Programme preparation:
    # Parse arguments
    # Check if namepsace isn't passed, if not parse command-line
    if argv is None:
        # Parse command-line
        parser = build_parser()
        args = parser.parse_args()
    else:
        args = build_parser(argv).parse_args()

    # Initiate logger
    # Note: log file only created if specified at cmdline
    if logger is None:
        logger = build_logger("Extract_genomes_NCBI", args)
    # logger = logging.getLogger("Extract_genomes_NCBI")
    logger.info("Run initated")

    # Add users email address from parser
    if args.user is None:
        logger.error(
            "No user email provided. Email MUST be providd. Terminating programme"
        )
        sys.exit(1)
    else:
        Entrez.email = args.user

    # If specified output directory for genomic files, create output directory
    if args.output is not sys.stdout:
        make_output_directory(args, logger)

    # Invoke main usage of programme
    # Create dataframe storing genus, species and NCBI Taxonomy ID, called 'species_table'
    species_table = parse_input_file(args.input_file, logger, args.retries)

    # Pull down accession numbers and GenBank files (if not disabled)
    species_table["NCBI Accession Numbers"] = species_table.apply(
        get_accession_numbers, args=(logger, args), axis=1
    )
    logger.info("Generated species table")

    # Write out dataframe
    if args.dataframe is not sys.stdout:
        write_out_dataframe(
            species_table, logger, args.dataframe, args.force, args.nodelete
        )
    else:
        species_table.to_csv(args.dataframe)

    # Program finished
    logger.info("Program finished and exiting")


def parse_input_file(input_filename, logger, retries):
    """Parse input file, returning dataframe of species names and NCBI Taxonomy IDs.

    Read input file. Calling functions as appropriate to retrieve scientific name
    or taxonomy ID if the other is given.

    Input file formating:
    Comments indicated with the first character of '#', perform no action.
    NCBI Taxonomy ID: Indicated with the first nine characters of 'NCBI:txid', pass
    to get_genus_species_name() to retrieve scientific name.
    Genus/species name: Indicated by lack of '#' and 'NCBI:txid', pass to get_tax_id()
    to return NCBI Taxonomy ID.

    Split genus and species name into separate variables, to be stored in separate
    dataframe columns.

    Generate a dataframe with three columns: 'Genus', 'Species' and 'NCBI Taxonomy ID'.

    :param input_filename: args, if specific name of input file, otherwise input taken from STDIN
    :param logger: logger object
    :param retries: parser argument, maximum number of retries excepted if network error encountered

    Return dataframe.
    """
    # open working_species_list.txt and extract lines, without newline character
    # then create genus name, species name and taxonomy ID tuplet
    all_species_data = []

    # test path to input file exists, if not exit programme
    if not input_filename.is_file():
        # report to user and exit programme
        logger.critical(
            (
                "Input file not found. Check filename, extension and directory is correct.\n"
                "Terminating program."
            ),
            exc_info=1,
        )
        sys.exit(1)

    # if path to input file exists proceed
    # parse input file
    with open(input_filename) as file:
        input_list = file.read().splitlines()

    # Parse input, retrieving tax ID or scientific name as appropriate
    line_count = 0
    for line in tqdm(input_list, desc="Reading lines"):
        line_count += 1

        if line.startswith("#"):
            continue

        line_data = parse_line(line, logger, line_count, retries)
        all_species_data.append(line_data)

    logger.info("Finished reading and closed input file")

    # create dataframe containing three columns: 'Genus', 'Species', 'NCBI Taxonomy ID'
    logger.info("Generating genus, species and taxonomy ID dataframe")
    species_table = pd.DataFrame(
        all_species_data, columns=["Genus", "Species", "NCBI Taxonomy ID"]
    )
    logger.info("Dataframe completed")
    return species_table


def parse_line(line, logger, line_count, retries):
    """Coordinate retrieval of scientific name or taxonomy ID.

    Read line from input file, calling functions as appropriate to retrieve
    scientific name or taxonomy ID.

    :param line: str, line from input file
    :param logger: logger object
    :line_count: number of line in input file - enable tracking if error occurs
    :param retries: parser argument, maximum number of retries excepted if network error encountered

    Return list of genus, species and taxonomy ID"""
    # For taxonomy ID retrieve scientific name
    if line.startswith("NCBI:txid"):
        gs_name = get_genus_species_name(line[9:], logger, line_count, retries)
        line_data = gs_name.split()
        line_data.append(line)
    # For scientific name retrieve taxonomy ID
    else:
        tax_id = get_tax_id(line, logger, line_count, retries)
        line_data = line.split()
        line_data.append(tax_id)

    return line_data


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
    # Retrieve scientific name
    with entrez_retry(
        logger, retries, Entrez.efetch, db="Taxonomy", id=taxonomy_id, retmode="xml"
    ) as handle:
        try:
            record = Entrez.read(handle)
        # if no record is returned from call to Entrez
        except TypeError:
            logger.error(
                (
                    f"Entrez failed to retrieve scientific name, for NCBI:txid{taxonomy_id}.\n"
                    "Potential typo in taxonomy ID, check input. Returned null value 'NA'."
                ),
                exc_info=1,
            )
            return "NA"

    # extract scientific name from record
    try:
        return record[0]["ScientificName"]

    except IndexError:
        logger.error(
            (
                f"Entrez failed to retrieve scientific name, for NCBI:txid{taxonomy_id}.\n"
                "Potential typo in taxonomy ID, check input. Returned null value 'NA'."
            ),
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
    if re.search(r"\d", genus_species):
        logger.warning(
            (
                f"Number with no 'NCBI:txid' prefix found line {line_number} ('{genus_species}').\n"
                "Maybe yypo in scientific name or ID missing 'NCBI:txid' prefix.\n"
                "Inturpretted as scientific name. Returning taxonomy ID as 'NA'.\n"
            ),
            exc_info=1,
        )
        return "NA"

    else:
        with entrez_retry(
            logger, retries, Entrez.esearch, db="Taxonomy", term=genus_species
        ) as handle:
            try:
                record = Entrez.read(handle)
            # if no record is returned from call to Entrez
            except TypeError:
                logger.error(
                    (
                        f"Entrez failed to retrieve scientific name, for NCBI:txid{taxonomy_id}.\n"
                        "Potential typo in taxonomy ID, check input. Returned null value 'NA'."
                    ),
                    exc_info=1,
                )
                return "NA"

    # extract taxonomy ID from record
    try:
        return "NCBI:txid" + record["IdList"][0]

    except IndexError:
        logger.error(
            (
                f"Entrez failed to retrieve taxonomy ID, for {genus_species}"
                "Potential typo in species name. Returning taxonomy ID as null value 'NA'"
            ),
            exc_info=1,
        )
        return "NA"


def get_accession_numbers(df_row, logger, args):
    """Return all NCBI accession numbers associated with NCBI Taxonomy ID.

    Use Entrez elink function to pull down the assembly IDs of all genomic
    assemblies in the NCBI Assembly database associated with the passed
    taxonomy ID.

    Use Entrez epost to post all assembly IDs to NCBI as a single query for
    subsequent Entrez efetch of all associated accession numbers.

    Format accession numbers into a human-readable list, 'all_accession_numbers'.

    Reminder of Pandas series structure:
    df_row[0]: Genus
    df_row[1]: Species
    df_row[2]: Taxonomy ID

    :param df_row: pd series, row from dataframe
    :param logger: logger object
    :param args: parser arguments

    Return list of NCBI accession numbers.
    """
    # If previously failed to retrieve the taxonomy ID cancel retrieval of accession numbers
    if df_row[2] == "NA":
        logger.warning(
            (
                f"Previously failed to retrieve taxonomy for {df_row[0][0]}.{df_row[1]}.\n"
                "Returning 'NA' for accession numbers."
            )
        )
        return "NA"

    # Retrieve all IDs of genomic assemblies for taxonomy ID

    logger.info(f"Retrieving assembly IDs for {df_row[2]}")
    # df_row[2][9:] removes 'NCBI:txid' prefix
    with entrez_retry(
        logger,
        args.retries,
        Entrez.elink,
        dbfrom="Taxonomy",
        id=df_row[2][9:],
        db="Assembly",
        linkname="taxonomy_assembly",
    ) as assembly_number_handle:
        assembly_number_record = Entrez.read(assembly_number_handle)

    # test record was returned, if failed to return exit retrieval of assembly IDs
    if assembly_number_record == "NA":
        logger.error(
            (
                f"Entrez failed to retrieve assembly IDs, for {df_row[2]}.\n"
                "Returning null value 'NA' for accession numbers."
            ),
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
                f"Entrez failed to retrieve assembly IDs, for {df_row[2]}."
                "Exiting retrieval of accession numbers, and returning null value 'NA'"
            ),
            exc_info=1,
        )
        return "NA"

    logger.info(f"Posting assembly IDs for {df_row[2]} to retrieve accession numbers")

    # compile list of ids in suitable format for epost
    id_post_list = str(",".join(assembly_id_list))
    # Post all assembly IDs to Entrez-NCBI for downstream pulldown of accession numbers
    epost_search_results = Entrez.read(
        entrez_retry(logger, args.retries, Entrez.epost, "Assembly", id=id_post_list)
    )

    # test record was returned, if failed to return exit retrieval of assembly IDs
    if epost_search_results == "NA":
        logger.error(
            (
                f"Entrez failed to post assembly IDs, for {df_row[2]}.\n"
                "Exiting retrieval of accession numbers, and returning null value 'NA'"
            ),
            exc_info=1,
        )
        return epost_search_results

    # Retrieve web environment and query key from Entrez epost
    epost_webenv = epost_search_results["WebEnv"]
    epost_query_key = epost_search_results["QueryKey"]

    logger.info(f"Retrieving accession numbers for {df_row[2]}")

    # Pull down all accession numbers
    ncbi_accession_numbers_list = []

    with entrez_retry(
        logger,
        args.retries,
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
                f"Entrez failed to retireve accession numbers, for {df_row[2]}."
                "Exiting retrieval of accession numbers, and returning null value 'NA'"
            ),
            exc_info=1,
        )
        return "NA"

    # Extract accession numbers from document summary
    for index_number in tqdm(
        range(len(accession_record["DocumentSummarySet"]["DocumentSummary"])),
        desc=f"Retrieving accessions ({df_row[2]})",
    ):
        try:
            new_accession_number = accession_record["DocumentSummarySet"][
                "DocumentSummary"
            ][index_number]["AssemblyAccession"]
            ncbi_accession_numbers_list.append(new_accession_number)

        except IndexError:
            total_assemblies = len(
                accession_record["DocumentSummarySet"]["DocumentSummary"]
            )
            logger.error(
                (
                    f"No accession number retrieved from NCBI for assembly ID {index_number}"
                    f"of {total_assemblies}.\n"
                    "Returning null value of 'NA'"
                ),
                exc_info=1,
            )
            return "NA"

        # If downloading of GenBank files is enabled, download Genbank files
        if args.genbank is True:
            get_genbank_files(
                new_accession_number,
                accession_record["DocumentSummarySet"]["DocumentSummary"][index_number][
                    "AssemblyName"
                ],
                logger,
                args,
            )

        index_number += 1

    logger.info(f"Finished processing retrieval accession numbers for {df_row[2]}")

    # Process accession numbers into human readable list for dataframe
    ncbi_accession_numbers = ", ".join(ncbi_accession_numbers_list)

    return ncbi_accession_numbers


def get_genbank_files(
    accession_number, assembly_name, logger, args, suffix="genomic.gbff.gz",
):
    """Coordiante download of GenBank from NCBI.

    :param accession_number: str, accession number
    :param assembly_name: str, name of assembly from NCBI record
    :param logger: logger object
    :param args: parser arguments
    :param suffix: str, suffix of file

    Return nothing.
    """
    # compile url for download
    genbank_url, filestem = compile_url(accession_number, assembly_name, logger, suffix)

    # if downloaded file is not to be written to STDOUT, compile output path
    if args.output is not sys.stdout:
        out_file_path = args.output / "_".join([filestem.replace(".", "_"), suffix])
    else:
        out_file_path = args.output

    # download GenBank file
    download_file(
        genbank_url, args, out_file_path, logger, accession_number, "GenBank file",
    )

    return


def compile_url(
    accession_number,
    assembly_name,
    logger,
    suffix,
    ftpstem="ftp://ftp.ncbi.nlm.nih.gov/genomes/all",
):
    """Compile url for file download.

    Reformate assembly name, replacing escape characters with underscors.
    This is because NCBI records can include a variety of escape characters
    in its records, but requires inclusion of underscores within its urls.

    Use ccession number and assembly name to generate file stem.

    Use the file stem to arquire the GCstem (region of the url which
    dictates the type of genome record, i.e. reference (GCF) or assembly (GCA)),
    and accession number block (in the format of nnn/nnn/nnn in the url).

    Regions of the url are compiled together with the ftpstem (prefix):
    "ftp://ftp.ncni.nlm.nih.gov/genomes/all/".

    :param accession_number: str, asseccion number of genomic assembly
    :param assembly_name: str, name of genomic assembly from NCBI record
    :param logger: logger object

    Return str, url required for download.
    """
    # Extract assembly name, removing alterantive escape characters
    escape_characters = re.compile(r"[\s/,#\(\)]")
    escape_name = re.sub(escape_characters, "_", assembly_name)

    # compile filstem
    filestem = "_".join([accession_number, escape_name])

    # separate out filesteam into GCstem, accession number intergers and discarded
    url_parts = tuple(filestem.split("_", 2))

    # separate identifying numbers from version number
    sub_directories = "/".join(
        [url_parts[1][i : i + 3] for i in range(0, len(url_parts[1].split(".")[0]), 3)]
    )

    # return url for downloading file
    return (
        "{0}/{1}/{2}/{3}/{3}_{4}".format(
            ftpstem, url_parts[0], sub_directories, filestem, suffix
        ),
        filestem,
    )


def download_file(
    genbank_url, args, out_file_path, logger, accession_number, file_type
):
    """Download file.

    :param genbank_url: str, url of file to be downloaded
    :param args: parser arguments
    :param out_file_path: path, output directory for file to be written to
    :param logger: logger object
    :param accession_number: str, accession number of genome
    :param file_type: str, denotes in logger file type downloaded

    Return nothing.
    """
    # Try URL connection
    try:
        response = urlopen(genbank_url, timeout=args.timeout)
    except HTTPError:
        logger.error(
            f"Failed to download {file_type} for {accession_number}", exc_info=1,
        )
        return
    except URLError:
        logger.error(
            f"Failed to download {file_type} for {accession_number}", exc_info=1,
        )
        return
    except timeout:
        logger.error(
            f"Failed to download {file_type} for {accession_number}", exc_info=1,
        )
        return

    if out_file_path.exists():
        logger.warning(f"Output file {out_file_path} exists, not downloading")
    else:
        # Download file
        logger.info("Opened URL and parsed metadata")
        file_size = int(response.info().get("Content-length"))
        bsize = 1_048_576
        try:
            with open(out_file_path, "wb") as out_handle:
                # Using leave=False as this will be an internally-nested progress bar
                with tqdm(
                    total=file_size,
                    leave=False,
                    desc=f"Downloading {accession_number} {file_type}",
                ) as pbar:
                    while True:
                        buffer = response.read(bsize)
                        if not buffer:
                            break
                        pbar.update(len(buffer))
                        out_handle.write(buffer)
        except IOError:
            logger.error(f"Download failed for {accession_number}", exc_info=1)
            return

        logger.info(
            f"Finished downloading GenBank file for {accession_number}", exc_info=1
        )

        return


def entrez_retry(logger, retries, entrez_func, *func_args, **func_kwargs):
    """Call to NCBI using Entrez.

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
                    f"Network error encountered during try no.{tries}.\nRetrying in 10s",
                    exc_info=1,
                )
                time.sleep(10)
            tries += 1

    if record is None:
        logger.error(
            "Network error encountered too many times. Exiting attempt to call to NCBI"
        )
        return

    return record


if __name__ == "__main__":
    main()
