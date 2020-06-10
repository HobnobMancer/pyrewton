#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Create summary of annotated CAZy classes in GenBank files.

:cmd_args df_input: path, path to input dataframe
:cmd_args force: bool, force overwriting files in output directory
:cmd_args genbank: path, path to directory containing GenBank files
:cmd_args log: path, path to direct writing out log file
:cmd_args nodelete: not delete existing files in output directory
:cmd_args output: path, path to output directory

:func build_parser: Build paser to allow cmd-line operation
:func main: Coordinate calling of other functions
:func build_logger: Build logger
:func check_input: Check paths to input dataframe and GenBank files are valid
:func make_output_directory: Establish output directory
:func get_input_df: parse input dataframe
:func create_dataframe: build dataframe summarising CAZy annotation in GenBank files
:func create_df_foundation: Parse input dataframe row
:func build_df_foundation: Compile row data for dataframe
:func get_protein_data: Retrieve protein name and IDs

Generate summary dataframe and of annotated CAZy classes in all GenBank
files associated with a given species.

Author:
Emma E. M. Hobbs

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
import gzip
import logging
import shutil
import sys

from pathlib import Path

import pandas as pd
import seaborn as sns

from Bio import SeqIO
from bioservices import UniProt
from tqdm import tqdm


def build_parser():
    """Return ArgumentParser parser for script."""
    # Create parser object
    parser = argparse.ArgumentParser(
        prog="cazy_genbank_summary.py",
        description="Generate summary of CAZy annotation in GenBank files",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # Add arguments to parser
    # Add option to specific input directory for dataframe
    parser.add_argument(
        "-d",
        "--df_input",
        type=Path,
        metavar="input datafram name",
        default=sys.stdin,
        help="input dataframe path",
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
    # Add option to specific input directory for GenBank files
    parser.add_argument(
        "-g",
        "--genbank",
        type=Path,
        metavar="GenBank file directory",
        default=sys.stdin,
        help="GenBank file path directory",
    )
    # Add option to specific directory for log to be written out to
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
        dest="nodelete",
        action="store_true",
        default=False,
        help="enable/disable deletion of exisiting files",
    )
    # Add option to specific directory for output to be written to
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        metavar="output file name",
        default=sys.stdout,
        help="output filename",
    )

    return parser


def main():
    """docstring summary.

    Detail.

    Return.
    """
    # Create parser object for cmd-line ctrl
    parser = build_parser()
    args = parser.parse_args()

    # Initiate logger
    # Note: log file only created if specified at cmdline
    build_logger("cazy_genbank_summary", args.log)
    logger = logging.getLogger("cazy_genbank_summary")
    logger.info("Run initated")

    # Check inputs are valid
    check_input(args, logger)
    logger.info("Inputs accepted")

    # If specified output directory for genomic files, create output directory
    if args.output is not sys.stdout:
        make_output_directory(args, logger)

    # Open input dataframe
    logger.info("Opening input dataframe")
    input_df = get_input_df(args.df_input, logger)

    # Build dataframe
    cazy_summary_df = create_dataframe(input_df, args, logger)

    # Create summary charts of CAZy annotation distribution

    # Write out dataframe

    return


def build_logger(script_name, log_file,) -> logging.Logger:
    """Return a logger for this script.

    Enables logger for script, sets parameters and creates new file to store log.

    :param script_name: str, name of script
    :param log_file: parser argument, enable writing out of log file

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


def check_input(args, logger):
    """Check paths to input dataframe and GenBank files is valid.

    :param args: parser arguments
    :param logger: logger object

    Return nothing if paths are valid.
    """
    logger.info("Checking path to input dataframe is valid")
    if (args.df_input).is_file() is False:
        logger.info(
            (
                "Input dataframe not found. Check filename, extension and directory are correct."
                "\nTerminating program."
            ),
            exc_info=1,
        )
        sys.exit(1)

    logger.info("Checking path to GenBank file containing directory is valid")
    if (args.genbank).exists is False:
        logger.info(
            (
                "GenBank file directory not found. Check correct directory was provided."
                "\nTerminating program."
            ),
            exc_info=1,
        )
        sys.exit(1)

    return


def make_output_directory(args, logger):
    """Create output directory for genomic files.

    Check if directory indicated for output existed already.
    If so check if force overwrite enabled. If not terminate programme.
    If so, check if deletion of exiting files was enabled.
    If so, exiting files in output directory are deleted.
    Create output directory, expecting error if already exists.

    :param args: parser arguments
    :param logger: logger object

    Return nothing.
    """
    logger.info("Checking if specified output directory for genomic files exists.")
    # If output directory specificed at cmd-line already exists, and 'force' not enabled
    if (args.output).exists():
        if (args.force) is False:
            logger.info(
                "Output directory already exists and forced overwrite not enabled.\n"
                "Terminating programme."
            )
            sys.exit()
        # If output directory exists and 'force' overwrite enabled
        else:
            # 'Nodelete' not enabled so delete output directory contents
            if (args.nodelete) is False:
                logger.info(
                    "Output directory already exists and forced complete overwrite enabled.\n"
                    "Deleting existing content in outdir."
                )
                # delete existing content in outdir
                shutil.rmtree(args.output)
            # 'Nodelete' enabled, don't empty directory
            else:
                logger.info(
                    "Output directory already exists and forced addition of files"
                    "to outdir enabled."
                )
    # Recursively make output directory
    try:
        (args.output).mkdir(exist_ok=True)
    except OSError:
        # this will occur if directory already exists
        # ignored if forced over write enabled
        if args.force is True:
            return ()
        else:
            logger.error(
                "OSError occured while creating output directory for genomic files.\n"
                "Terminating programme."
            )
            sys.exit()
    return


def get_input_df(input_df, logger):
    """Open input dataframe (df).

    Input dataframe must contain at least columns titled:
    'Genus', 'Species', 'NCBI Taxonomy ID', and 'NCBI Accession Numbers'.

    Return dataframe.
    """
    input_df = pd.read_csv(
        input_df,
        header=0,
        names=["Genus", "Species", "NCBI Taxonomy ID", "NCBI Accession Numbers"],
    )
    return input_df


def create_dataframe(input_df, args, logger):
    """Build datafame.

    Iterate over input dataframe row wise. This allows for converting human
    readable accession number list from a string to a Python list.

    Per species, retrieve all protein names and IDs for every accession number,
    therefore, will return multiple rows with the same accession number,
    but unique protein ID.

    Append data for all rows for a single species to the tuple
    'all_foundation_data', so that all data it compiled together and can be
    used simultaneously to populate a pandas dataframe without risk of
    data missalignment or over writing.

    Then iterate over protein ID column to retrieve CAZy data if available,
    allowing for faster retrieval of data from UniProt.

    Any failed to retrieve data will be returned as pandas null value 'NA'.

    :param input_df: pandas dataframe
    :param args: parser arguments
    :param logger: logger object

    Return dataframe.
    """
    # Create empty dataframe to add data to
    cazy_summary_df = pd.DataFrame(
        columns=[
            "Genus",
            "Species",
            "NCBI Taxonomy ID",
            "NCBI Accession Number",
            "NCBI Protein ID",
            "Locus Tag",
            "Gene Locus",
            "NCBI Recorded Function",
        ]
    )

    # Retrieve data for dataframe foundation and add to empty dataframe
    cazy_summary_df = pd.concat(
        (input_df.apply(get_df_foundation_data, args=(args, logger), axis=1)),
        ignore_index=True,
    )

    print(cazy_summary_df)

    # Add CAZy data to dataframe
    # if cazy class returned full section will be titled 'cazy class',
    # if familied returned use 'cazy family' instead
    # CAZy_summary_df["Cazy Class", "Function"] = CAZy_summary_df.apply(
    #     lambda column: get_cazy_data(column["Protein ID"], logger), axis=1
    # )


def get_df_foundation_data(df_row, args, logger):
    """Prepare row data to create dataframe.

    Coordinate retrieval of protein data from GenBank files for every accession
    number for the species passed to the function.

    Store data in a dataframe: Genus, Species, Tax ID, Accession number,
    protein ID, locus tag, gene location, product.

    Each row in dataframe contains a unique protein, and thus multiple rows
    will have the same accession number.

    Any failed to retrieve data will be returned as pandas null value 'NA'.

    Reminder of panda series (referred to as df_row) structure:
    df_row[0] = "Genus"
    df_row[1] = "Species"
    df_row[2] = "Taxonomy ID"
    df_row[3] = "Accession list" - this is human readable list, stored as a str.

    :param df_row: row from input_df (dataframe)
    :param args: parser arguments
    :param logger: logger object

    Return dataframe.
    """
    # Create empty dataframe to store data in
    protein_data_df = pd.DataFrame(
        columns=[
            "Genus",
            "Species",
            "NCBI Taxonomy ID",
            "NCBI Accession Number",
            "NCBI Protein ID",
            "Locus Tag",
            "Gene Locus",
            "NCBI Recorded Function",
        ]
    )

    logger.info(
        (
            "Adding scientific name, accession numbers, protein names and IDs\n"
            f"to dataframe for {df_row[0][0]}.{df_row[1]}, {df_row[2]}"
        )
    )

    # convert human readable list of accession numbers into Python list
    accession_list = df_row[3].split(", ")

    # open GenBank file of each accession number in the list and retrieve
    # all protein data in that GenBank file, stored as a tuple with each
    # list in the tuple containing data for a unique protein
    for accession in accession_list:
        logger.info(f"Retrieving data for {accession}")
        protein_data = get_protein_data(accession, args.genbank, logger)

        # check if any data retrieved
        if len(protein_data) == 0:
            logger.warning(
                (
                    f"No protein data retrieved for {accession} from GenBank file.\n"
                    "Most likely cause is GenBank file contained no CDS type features."
                )
            )
            # Add null values to dataframe for the accession number
            new_df_row = {
                "Genus": df_row[0],
                "Species": df_row[1],
                "NCBI Taxonomy ID": df_row[2],
                "NCBI Accession Number": accession,
                "NCBI Protein ID": "NA",
                "Locus Tag": "NA",
                "Gene Locus": "NA",
                "NCBI Recorded Function": "NA",
            }
            protein_data_df.append(new_df_row, ignore_index=True)

        # if data was returned add to list of all row data
        else:
            # For each unique protein in the GenBank file create a new row in
            # dataframe. The data for each unique protein is stored as a single
            # list in the tuple protein_data
            protein_index = 0  # index number in protein_data tuple
            for protein_index in range(len(protein_data)):
                logger.info("adding protein data to row data")

                # Compile data for new row to be added to dataframe
                new_df_row = {
                    "Genus": df_row[0],
                    "Species": df_row[1],
                    "NCBI Taxonomy ID": df_row[2],
                    "NCBI Accession Number": accession,
                    "NCBI Protein ID": protein_data[protein_index][0],
                    "Locus Tag": protein_data[protein_index][1],
                    "Gene Locus": protein_data[protein_index][2],
                    "NCBI Recorded Function": protein_data[protein_index][3],
                }

                # Add new row to dataframe
                protein_data_df.append(new_df_row, ignore_index=True)

                protein_index += 1
                logger.info(f"finished collecting data for {accession}")

    return protein_data_df


def get_protein_data(accession_number, genbank_input, logger):
    """Retrieve protein ID, locus tag and function from GenBank file.

    From each record the protein ID, locus tag, location and annotated
    function is retrieved, and stored as a list.

    Lists wil be added to a single tuple containing all protein data.

    Any failed to retrieve data will be returned as pandas null value 'NA'.

    :param accession_number: str
    :param genbank_input: path, path to directory containing GenBank files
    :param logger: logger object

    Return tuple.
    """
    # check if accession number was provided
    if accession_number == "NA":
        logger.warning(
            (
                f"Null value ('NA') was contained in cell for {accession_number},"
                "exiting retrieval of protein data.\nReturning null ('NA') value"
                "for all protein data"
            )
        )
        return ["NA", "NA", "NA", "NA"]

    # replace '.' with '_' to match format in GenBank file name
    file_stem = accession_number.replace(".", "_")

    # create empty list to store file entries, to allow checking if multiple files were retrieved
    gb_file = []

    # retrieve all files from directory
    files_in_entries = (
        entry for entry in Path(genbank_input).iterdir() if entry.is_file()
    )
    for item in files_in_entries:
        # search for accession number's GenBank file
        if item.name.startswith(f"{file_stem}") and item.name.endswith(".gbff.gz"):
            gb_file.append(item)

    # check file was retrieved, not multiple or none
    if len(gb_file) == 0:
        logger.warning(
            (
                f"Retrieved 0 files for {accession_number}.\n"
                "Returning null ('NA') value for all protein data"
            )
        )

    elif len(gb_file) > 1:
        logger.warning(
            (
                f"Retrieved multiple files for {accession_number}.\n"
                "Returning null ('NA') value for all protein data"
            )
        )
        return ["NA", "NA", "NA", "NA"]

    # check if files is empty
    if gb_file[0].stat().st_size == 0:
        logger.warning(
            (
                f"GenBank file retrieved for {accession_number} is empty.\n"
                "Returning null ('NA' value for all protein data"
            )
        )
        return ["NA", "NA", "NA", "NA"]

    all_protein_data = []
    # Retrieve protein data
    with gzip.open(gb_file[0], "rt") as handle:
        logger.info(f"opening gb file for {accession_number}")
        # create list to store all protein data retrieved from GenBank file, making it a tuple
        for gb_record in SeqIO.parse(handle, "genbank"):
            for (index, feature) in enumerate(gb_record.features):
                # empty protein data list so as not to contaminate data of next protein
                protein_data = []
                # Parse over only protein encoding features (type = 'CDS')
                if feature.type == "CDS":
                    # extract protein ID
                    logger.info("retrieving protein data for CDS feature")
                    protein_data.append(
                        get_record_feature(feature, "protein_id", logger)
                    )
                    # extract locus tag
                    protein_data.append(
                        get_record_feature(feature, "locus_tag", logger)
                    )
                    # extract location
                    protein_data.append(get_record_feature(feature, "location", logger))
                    # extract annotated function of product
                    protein_data.append(get_record_feature(feature, "product", logger))

                    # add protein data to total protein data list, only if data was retrieved
                    if len(protein_data) == 4:
                        if protein_data != ["NA", "NA", "NA", "NA"]:
                            all_protein_data.append(protein_data)
                        else:
                            logger.warning(
                                f"No data retrieved from CDS type feature, index: {index}",
                                exc_info=1,
                            )
                    else:
                        logger.warning(
                            (
                                f"Error occured during retrieval of data from feature, {index}\n"
                                f"for {accession_number}. Returning no protein data"
                            )
                        )

    return all_protein_data


def get_record_feature(feature, qualifier, logger):
    """Retrieve data from GenBank record feature.

    :param feature: feature object, GenBank file record feature
    :param qualifier: str, key of feature attribute
    :param logger: logger object

    Return data from GenBank record feature, or "NA" if failed to retrieve.
    """
    # if called to extract location, extract location as human readable list
    if qualifier == "location":
        try:
            location_list = []
            for item in feature.location.parts:
                location_list.append(str(item))
            compiled_location = str(",".join(location_list))
            return compiled_location
        except AttributeError:
            logger.warning(
                "Failed to retrieve feature location, returning 'NA'", exc_info=1
            )
            return "NA"
    else:
        try:
            data = feature.qualifiers[qualifier][0]
            return data
        except KeyError:
            logger.warning(
                f"Failed to retrieve feature {qualifier}, returning 'NA'", exc_info=1
            )
            return "NA"


def get_cazy_data(protein_id, logger):
    # use protein_id to call to UniProt
    # if CAZy link return class - may return full family rather than class - no worries
    # if no CAZy link return 'NA'
    return  # list, first CAZy class and then protein function


def create_summary_chart(cazy_fam_column, logger):
    # Use seaborn to create summary chart
    return


def write_out_df(df, output, logger):
    # write out df to specified directory
    return


def write_out_chart(chart, output, logger):
    # write out chart to specified directory
    return


if __name__ == "__main__":
    main()
