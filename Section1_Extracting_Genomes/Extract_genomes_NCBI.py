#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Pull down genomic assemblies from NCBI database.

:cmd_args user_email: email address of user
:cmd_args --dataframe: output directory for dataframe
:cmd_args --force: force writing in output directory
:cmd_args --genbank: enable/disable GenBank file download
:cmd_args --input_file: path of input file
:cmd_args --log: enable log file writing
:cmd_args --nodelete: enable/disabling removal of exiting files in output
:cmd_args --output: output directory for downloaded files
:cmd_args --retries: maximum number of retries if network error occurs
:cmd_args --timeout: timeout limit of URL connection

:func main: generate a dataframe of scientific names, taxonomy IDs and accession numbers
:func build_logger: creates logger object
:func make_output_directory: create directory for genomic files to be written to
:func parse_input_file: parse input file
:func get_genus_species_name: retrieve scientific name from taxonomy ID
:func get_tax_id: retrieve NCBI taxonomy ID from scientific name
:func collate_accession_numbers: parse taxonomy ID column in dataframe
:func get_accession_numbers: retrieves all accessions associated to given taxonomy ID
:func get_genbank_files: organise download of genbank files
:func compile_URL: create URL for downloading file
:func compile_output_path: create file name and path to output directory for downloaded files
:func download_file: download file using provided URL
:func write_out_dataframe: write out species table as .csv file

Generates dataframe containing scientific names, taxonomy IDs and accession numbers.
Pulls down and stores genomic assemblies and GenBank filesfrom NCBI Assembly database.

Author:
Emma Hobbs

Contact
eemh1@st-andrews.ac.uk

Emma E. M. Hobbs,
Biomolecular Sciences Building,
University of St Andrews,
North Haugh Campus,
St Andrews,
KY16 9ST
Scotland,
UK

The MIT License
"""

import argparse
import datetime
import logging
import re
import shutil
import sys
import time

from pathlib import Path
from socket import timeout
from urllib.error import HTTPError, URLError
from urllib.request import urlopen

import pandas as pd

from Bio import Entrez
from tqdm import tqdm


def build_parser():
    """Return ArgumentParser parser for script."""
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
    # Add dataframe write out options
    # if not given dataframe written to STDOUT
    parser.add_argument(
        "-d",
        "--dataframe",
        type=Path,
        default=sys.stdout,
        help="Location of file for species table to be written to",
    )
    # Add option to force file over writting
    parser.add_argument(
        "-f",
        "--force",
        dest="force",
        action="store_true",
        default=False,
        help="Force file over writting",
    )
    # Add option to disable pull down of GenBank files
    parser.add_argument(
        "-g",
        "--genbank",
        type=bool,
        metavar="genbank file (.gbff)",
        default=True,
        help="Disable pulldown of GenBank (.gbff) files",
    )
    # Add input file name option
    # If not given input will be taken from STDIN
    parser.add_argument(
        "-i",
        "--input_file",
        type=Path,
        metavar="input file name",
        default=sys.stdin,
        help="Input filename",
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
    # Add option to prevent over writing of existing files
    # and cause addition of files to output directory
    parser.add_argument(
        "-n",
        "--nodelete",
        type=bool,
        metavar="not to delete exisiting files",
        default=False,
        help="enable/disable deletion of exisiting files",
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
        help="Output filename",
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
    # Add option to alter timeout allowance before cancelling downloading of files
    parser.add_argument(
        "-t",
        "--timeout",
        dest="timeout",
        action="store",
        default=10,
        help="Timeout for URL connections, in seconds",
    )

    return parser


def main():
    """Generate dataframe containing scientific names, taxonomy IDs and accession numbers.

    Parse input file to create dataframe ('species_table') of scientific name and taxonomy ID.

    Pass the 'species_table' dataframe to collate_accession_numbers(), to retrieve all associated
    NCBI accession numbers. Store accession numbers in new column.

    If enabled, write out dataframe to .csv file.

    Return 'species_table' dataframe.
    """
    # Capture date and time script is executed
    date = datetime.datetime.now()
    date_of_pulldown = date.strftime("%Y-%m-%d")
    time_of_pulldown = date.strftime("%H:%M")

    # Parse arguments
    parser = build_parser()
    args = parser.parse_args()

    # Add users email address from parser
    Entrez.email = args.user_email

    # Initiate logger
    # Note: log file only created if specified at cmdline
    logger = build_logger(
        "Extract_genomes_NCBI", args.log
    )
    # logger = logging.getLogger("Extract_genomes_NCBI")
    logger.info("Run initated")

    # If specified output directory for genomic files, create output directory
    if args.output is not sys.stdout:
        make_output_directory(args.output, logger, args.force, args.nodelete)

    # Invoke main usage of script
    # Create dataframe storing genus, species and NCBI Taxonomy ID, called 'species_table'
    species_table = parse_input_file(args.input_file, logger, args.retries)

    # Pull down accession numbers and GenBank files (if not disabled)
    species_table = collate_accession_numbers(species_table, logger, args)
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


def build_logger(script_name, log_file) -> logging.Logger:
    """Return a logger for this script.

    Enables logger for script, sets parameters and creates new file to store log.

    :param script_name: str, name of script
    :param log_file: parser argument, enable writing out of log file

    Return logger object.
    """
    logger = logging.getLogger(script_name)
    logger.setLevel(logging.WARNING)
    # Set format of loglines
    log_formatter = logging.Formatter(
        script_name + ": {} - {}".format("%(asctime)s", "%(message)s")
    )

    # Setup console handler to log to terminal
    console_log_handler = logging.StreamHandler()
    console_log_handler.setLevel(logging.WARNING)
    console_log_handler.setFormatter(log_formatter)
    logger.addHandler(console_log_handler)

    # Setup file handler to log to a file
    if log_file is not None:
        file_log_handler = logging.FileHandler(log_file)
        file_log_handler.setLevel(logging.INFO)
        file_log_handler.setFormatter(log_formatter)
        logger.addHandler(file_log_handler)

    return logger


def make_output_directory(output, logger, force, nodelete):
    """Create output directory for genomic files.

    Check if directory indicated for output existed already.
    If so check if force overwrite enabled. If not terminate programme.
    If so, check if deletion of exiting files was enabled.
    If so, exiting files in output directory are deleted.
    Create output directory, expecting error if already exists.

    :param output: Path, path to output directory
    :param logger: logger object
    :param force: bool, cmd-line args to enable/disable over writing existing directory
    :param nodelete: bool, cmd-line args to enable/disable deleting of existing files

    Return nothing.
    """
    logger.info("Checking if specified output directory for genomic files exists.")
    # If output directory specificed at cmd-line, check output directory does not already exist
    if output.exists():
        if force is False:
            logger.info(
                (
                    "Output directory already exists and forced overwrite not enabled.\n"
                    "Terminating program."
                )
            )
            sys.exit(1)
        else:
            if nodelete is False:
                logger.info(
                    (
                        "Output directory already exists and forced complete overwrite enabled.\n"
                        "Deleting existing content in outdir."
                    )
                )
                # delete existing content in outdir
                shutil.rmtree(output)
            else:
                logger.info(
                    (
                        "Output directory already exists and forced addition of files "
                        "to outdir enables."
                    )
                )
    # Recursively make output directory (directory may exist if
    # --force==True and --nodelete==True, so exist_ok is required)
    output.mkdir(exist_ok=True)
    return


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
        logger.info(
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
    # logger.info("(Retrieving scientific name for NCBI:txid{})".format(taxonomy_id))

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
                "Entrez failed to retrieve scientific name, for species in line {} of input file.\n"
                "Potential typo in taxonomy ID, check input. Returned 'NA'."
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
    if re.search(r"\d", genus_species):
        logger.warning(
            (
                "Warning: Number with no 'NCBI:txid' prefix found line {}.\n"
                "Typo in scientific name of ID missing 'NCBI:txid' prefix.\n"
                "Inturpretted as scientific name. Returning taxonomy ID as 'NA'.\n"
            ).format(line_number),
            exc_info=1,
        )
        return "NA"

    else:
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


def collate_accession_numbers(species_table, logger, args):
    """Return dataframe with column containing all associated NCBI accession numbers.

    Pass each Taxonomy ID in the 'NCBI Taxonomy ID' column of the 'species_table'
    dataframe to find all associated NCBI accession numbers.

    :param species_table: dataframe, dataframe containing scientific names and taxonomy IDs
    :param logger: logger object
    :param args: parser arguments

    Return modified dataframe, with four columns.
    """
    logger.info(
        "Aquiring accession numbers from NCBI Assembly database for NCBI Taxonomy IDs"
    )
    species_table["NCBI Accession Numbers"] = species_table.apply(
        lambda column: get_accession_numbers(
            column["NCBI Taxonomy ID"][9:], column, logger, args.retries, args
        ),
        axis=1,
    )

    logger.info("Finished adding accession numbers to dataframe")
    return species_table


def get_accession_numbers(taxonomy_id, df_row, logger, retries, args):
    """Return all NCBI accession numbers associated with NCBI Taxonomy Id.

    Use Entrez elink function to pull down the assembly IDs of all genomic
    assemblies in the NCBI Assembly database associated with the passed
    taxonomy ID.

    Use Entrez epost to post all assembly IDs to NCBI as a single query for
    subsequent Entrez efetch of all associated accession numbers.

    Format accession numbers into a human-readable list, 'all_accession_numbers'.

    :param taxonomy_id: str, NCBI taxonomy ID
    :param df_row: pd series, row from dataframe
    :param logger: logger object
    :param retries: parser argument, maximum number of retries excepted if network error encountered
    :param args: parser arguments

    Return list of NCBI accession numbers.
    """
    # If previously failed to retrieve the taxonomy ID cancel retrieval of accession numbers
    if taxonomy_id == "NA":
        logger.info(
            (
                "Previously failed to retrieve taxonomy for {} {}.\n"
                "Returning 'NA' for accession numbers."
            ).format(df_row[0], df_row[1])
        )
        return "NA"

    # Retrieve all IDs of genomic assemblies for taxonomy ID

    logger.info("Retrieving assembly IDs for NCBI:txid{}".format(taxonomy_id))

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
                "Entrez failed to retrieve assembly IDs, for NCBI:txid{}.\n"
                "Returning 'NA' for accession numbers for NCBI:txid{}"
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
                "Exiting retrieval of accession numbers for NCBI:txid{}"
            ).format(taxonomy_id, taxonomy_id),
            exc_info=1,
        )
        return "NA"

    logger.info(
        "Finished processing retrieval of assembly IDs for NCBI:txid{}".format(
            taxonomy_id
        )
    )

    # compile list of ids in suitable format for epost
    id_post_list = str(",".join(assembly_id_list))
    # Post all assembly IDs to Entrez-NCBI for downstream pulldown of accession numbers
    epost_search_results = Entrez.read(
        entrez_retry(logger, retries, Entrez.epost, "Assembly", id=id_post_list)
    )

    # test record was returned, if failed to return exit retrieval of assembly IDs
    if epost_search_results == "NA":
        logger.error(
            (
                "Entrez failed to post assembly IDs, for NCBI:txid{}."
                "Exiting retrieval of accession numbers for NCBI:txid{}"
            ).format(taxonomy_id, taxonomy_id),
            exc_info=1,
        )
        return epost_search_results

    # Retrieve web environment and query key from Entrez epost
    epost_webenv = epost_search_results["WebEnv"]
    epost_query_key = epost_search_results["QueryKey"]

    logger.info(
        "Finished processing positing of assembly IDs for fetching accession numbers"
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
                "Entrez failed to retireve accession numbers, for NCBI:txid{}."
                "Exiting retrieval of accession numbers for NCBI:txid{}"
            ).format(taxonomy_id, taxonomy_id),
            exc_info=1,
        )
        return accession_record

    # Extract accession numbers from document summary

    # Find total number of accession numbers
    number_of_accession_numbers = len(
        accession_record["DocumentSummarySet"]["DocumentSummary"]
    )

    # Retrieve accession number for assembly ID
    for index_number in tqdm(
        range(number_of_accession_numbers), desc=f"Retrieving accessions (txid:{taxonomy_id})"
    ):
        try:
            new_accession_number = accession_record["DocumentSummarySet"][
                "DocumentSummary"
            ][index_number]["AssemblyAccession"]
            ncbi_accession_numbers_list.append(new_accession_number)

        except IndexError:
            logger.error(
                "No accession number retrieved from NCBI for assembly ID {} of {}".format(
                    index_number, number_of_accession_numbers
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
                taxonomy_id,
                index_number,
                number_of_accession_numbers,
                logger,
                retries,
                args,
            )

        index_number += 1

    logger.info(
        "Finished processing retrieval accession numbers for NCBI:txid:{}".format(
            taxonomy_id
        )
    )

    # Process accession numbers into human readable list for dataframe
    ncbi_accession_numbers = ", ".join(ncbi_accession_numbers_list)

    return ncbi_accession_numbers


def get_genbank_files(
    accession_number,
    assembly_name,
    taxonomy_id,
    id_count,
    total_id_count,
    logger,
    retries,
    args,
):
    """Coordiante download of GenBank from NCBI.

    :param accession_number: str, accession number
    :param assembly_name: str, name of assembly from NCBI record
    :param taxonomy_id: str, taxonomy ID of host species
    :param id_count: int, index number of associated assembly ID in list of all IDs
    :param total_id_count: total number of retrieved assembly IDs for species
    :param logger: logger object
    :param retries: parser args, defines maximum number of retries if network error encountered
    :param args: parse arguments

    Return nothing.
    """
    logger.info(
        "Preparing to download GenBank file for assembly {} of {}".format(
            id_count, total_id_count
        )
    )

    # compile url for download
    genbank_url, filestem = compile_url(accession_number, assembly_name, logger)

    # if downloaded file is not to be written to STDOUT, compile output path
    if args.output is not sys.stdout:
        out_file_path = compile_output_path(
            args.output, filestem, "genomic.gbff.fz", logger
        )
    else:
        out_file_path = args.output

    # download GenBank file
    download_file(
        genbank_url, args.timeout, out_file_path, logger, accession_number,
    )

    logger.info(
        "Finished downloading GenBank file for assembly {} of {}".format(
            id_count, total_id_count
        )
    )

    return


def compile_url(accession_number, assembly_name, logger):
    """Compile url for file download.

    Reformate assembly name, replacing escape characters with underscors.
    This is because NCBI records can include a variety of escape characters
    in its records, but requires inclusion of underscores within its urls.

    Use ccession number and assembly name to generate file stem.

    Use the file stem to arquire the GCstem (region of the url which
    dictates the type of genome record, i.e. reference or assembly), and
    accession number block (in the format of nnn/nnn/nnn in the url).

    Regions of the url are compiled together with the ftpstem (prefix):
    "ftp://ftp.ncni.nlm.nih.gov/genomes/all/".

    :param accession_number: str, asseccion number of genomic assembly
    :param assembly_name: str, name of genomic assembly from NCBI record
    :param logger: logger object

    Return str, url required for download.
    """
    logger.info("Compiling URL for download")
    # Extract assembly name, removing alterantive escape characters
    escape_characters = re.compile(r"[\s/,#\(\)]")
    escape_name = re.sub(escape_characters, "_", assembly_name)
    # compile filstem
    filestem = "_".join([accession_number, escape_name])
    # separate out filesteam into GCstem, intergers and version number
    gcstem, accession_block, _ = tuple(filestem.split("_", 2))
    # separate identifying numbers from version number
    accession_intergers = accession_block.split(".")[0]
    # generate accession number block for url in format (nnn\nnn\nnn)
    url_accession_block = "/".join(
        [accession_intergers[i : i + 3] for i in range(0, len(accession_intergers), 3)]
    )
    # return url for downloading file
    return (
        f"ftp://ftp.ncbi.nlm.nih.gov/genomes/all/{gcstem}/{url_accession_block}/{filestem}/{filestem}_genomic.gbff.gz",
        filestem,
    )


def compile_output_path(output, filestem, suffix, logger):
    """Compile path for output file.

    Standardised file name of accession_number.suffix.
    Suffix for GenBank file is _genomic.gbff.gz

    :param output: cmd-args, path to directory for download files to be written
    :param filestem: str, file stem of GenBank file
    :param suffix: str, file extension
    :param logger: logger object

    Return path.
    """
    logger.info("Compiling output path for file")
    return output / "_".join([filestem.replace(".", "_"), suffix])


def download_file(genbank_url, timeout_limit, out_file_path, logger, accession_number):
    """Download file.

    :param genbank_url: str, url of file to be downloaded
    :param timeout_limit: int, cmd-line args, timout out of URL connection (s)
    :param out_file_path: path, output directory for file to be written to
    :param accession: str, accession number of genome

    Return nothing.
    """
    # Try URL connection
    try:
        response = urlopen(genbank_url, timeout=timeout_limit)
    except HTTPError:
        logger.error(
            "Failed to download GenBank file for {}".format(accession_number),
            exc_info=1,
        )
        return
    except URLError:
        logger.error(
            "Failed to download GenBank file for {}".format(accession_number),
            exc_info=1,
        )
        return
    except timeout:
        logger.error(
            "Download timed out, thus failed to download GenBank file for {}".format(
                accession_number
            ),
            exc_info=1,
        )
        return

    # Download file
    logger.info("Opened URL and parsed metadata")
    file_size = int(response.info().get("Content-length"))
    bytes_downloaded = 0
    bsize = 1_048_576
    try:
        logger.info("Downloading file. Accession: {}".format(accession_number))
        with open(out_file_path, "wb") as out_handle:
            # Using leave=False as this will be an internally-nested progress bar
            with tqdm(total=file_size, leave=False, desc=accession_number) as pbar:
                while True:
                    buffer = response.read(bsize)
                    if not buffer:
                        break
                    pbar.update(len(buffer))
    except IOError:
        logger.error("Download failed for {}".format(accession_number), exc_info=1)
        return

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


def write_out_dataframe(species_table, logger, outdir, force, nodelete):
    """Write out dataframe to output directory.

    :param species_table: pandas dataframe
    :param logger: logger object
    :param outdir: cmd-args, Path, output directory
    :param force: booleon, cmd-line argument to enable/disable over writing of existing files
    :param nodelete: boolean, cmd-line args to enable/disable deleting of existing files in outdir

    return Nothing.
    """
    # Check if overwrite of existing directory will occur
    logger.info("Checking if output directory for dataframe already exists")
    if outdir.exists():
        if force is False:
            logger.warning(
                "Specified directory for dataframe already exists.\nExiting writing out dataframe."
            )
            return ()
        else:
            logger.warning(
                "Specified directory for dataframe already exists.\nForced overwritting enabled."
            )
    logger.info("Writing out species dataframe to directory")
    species_table.to_csv(outdir)
    return ()


if __name__ == "__main__":
    main()
