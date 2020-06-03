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
        prog="CAZy_genbank_summary.py",
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
    build_logger("CAZy_genbank_summary", args.log)
    logger = logging.getLogger("CAZy_genbank_summary")
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
    CAZy_summary_df = create_dataframe(input_df, args, logger)

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
    if (args.input_df).is_file() is False:
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
    # Retrieve data for dataframe foundation: genus, species, accession number, protein name,
    # protein ID
    all_foundation_data = []  # empty list to store all data for foundation dataframe
    all_foundation_data.append(
        input_df.apply(get_df_foundation_data, args=(args, logger)), axis=1
    )

    # create foundation dataframe
    CAZy_summary_df = pd.Dataframe(
        all_foundation_data,
        columns=[
            "Genus",
            "Species",
            "NCBI Accession Number",
            "NCBI Protein ID",
            "Locus Tag",
            "NCBI Recorded Function",
        ],
    )

    # Add CAZy data to dataframe
    # if cazy class returned full section will be titled 'cazy class',
    # if familied returned use 'cazy family' instead
    CAZy_summary_df["Cazy Class", "Function"] = CAZy_summary_df.apply(
        lambda column: get_cazy_data(column["Protein ID"], logger), axis=1
    )


def get_df_foundation_data(df_row, args, logger):
    """Prepare row data to create dataframe.

    Retrieve all row data for single species, as a tuple. Each row is represented
    as a unique list in tuple. Each row/list containing genus, species, accession
    number, protein name and protein ID, with a unique protein name and ID per
    row/list.

    Will return multiple rows with the same accession number, but unique
    protein ID.

    Any failed to retrieve data will be returned as pandas null value 'NA'.

    Reminder of panda series (referred to as df_row) structure:
    df_row[0] = "Genus"
    df_row[1] = "Species"
    df_row[2] = "Taxonomy ID"
    df_row[3] = "Accession list" - this is human readable list, stored as a str.

    :param df_row: row from input_df (dataframe)
    :param args: parser arguments
    :param logger: logger object

    Return tuple.
    """
    logger.info(
        (
            "Adding scientific name, accession numbers, protein names and IDs\n"
            f"to dataframe for {df_row[0][0]}.{df_row[1]}, {df_row[2]}"
        )
    )

    # Gather row data
    row_foundation = []  # genus, species and accession
    # complete row data includes genus, species, accession, protein name, protein ID
    all_rows_data = []  # data for all rows

    accession_list = df_row[3].split(", ")

    count = 1
    for accession in tqdm(accession_list, desc="Compiling data"):
        row_foundation.append([df_row[0], df_row[1], accession])
        protein_data = get_protein_data(accession, args.genbank, logger)
        # construct data for completed row, then compile all dataframe data
        index = 0
        for index in range(len(protein_data)):
            # add protein ID, locus tag and function as individual items so as not to create a tuple
            all_rows_data.append(
                row_foundation.append(
                    protein_data[index][0],
                    protein_data[index][1],
                    protein_data[index][2],
                )
            )
            index += 1
            count += 1

    return all_rows_data


def get_protein_data(accession_number, genbank_input, logger):
    """Retrieve protein ID, locus tag and function from GenBank file.

    From each record the protein ID, locus tag and function will be retrieved,
    and stored as a list.

    Lists wil be added to a single tuple containing all protein data.

    Any failed to retrieve data will be returned as pandas null value 'NA'.

    :param accession_number: str
    :param genbank_input: path, path to directory containing GenBank files
    :param logger: logger object

    Return tuple.
    """
    # use accession number to open associated genbank file
    # use SeqIO to extract protein ID and protein name
    # return as list for dataframe
    # will need to check if protein ID already present in dataframe, and is don't add

    # check if accession number was provided
    if accession_number == "NA":
        logger.warning(
            (
                f"Null value ('NA') was contained in cell for {accession_number},"
                "exiting retrieval of protein data.\nReturning null ('NA') value."
            )
        )
        return ["NA", "NA", "NA"]

    # Retrieve GenBank (gb) file
    gb_file = list(Path(genbank_input).glob(rf"{accession_number}*.gbff.fz"))

    # check file was retrieved, not multiple or none
    if len(gb_file) > 1 or len(gb_file) == 0:
        # log error and return 'NA' for protein name and protein ID
        logger.warning(
            (
                f"Failed to retrieve GenBank file for {accession_number}.\n"
                "Returning null ('NA') value for protein name and ID"
            )
        )
        return ["NA", "NA", "NA"]

    else:
        # Retrieve protein data
        with gzip.open("GCA_001515345_1_ASM151534v1_genomic.gbff.gz", "rt") as handle:
            # create list to store all protein data retrieved from GenBank file, making it a tuple
            all_protein_data = []
            for gb_record in SeqIO.parse(handle, "genbank"):
                # create a list to store data for a single protein
                # make sure list is empty before parsing first feature
                protein_data = []
                for (index, feature) in enumerate(gb_record.features):
                    # empty protein data list so as not to contaminate data of next protein
                    protein_data = []
                    # Parse over only protein encoding features (type = 'CDS')
                    if feature.type == "CDS":
                        protein_data.append(feature.qualifiers["protein_id"][0])
                        protein_data.append(feature.qualifiers["locus_tag"][0])
                        protein_data.append("; ".join(feature.qualifiers["product"]))
                        # add protein data to total protein data list
                        all_protein_data.append(protein_data)

    return all_protein_data


def get_cazy_data(protein_id, logger):
    # use protein_id to call to UniProt
    # if CAZy link return class - may return full family rather than class - no worries
    # if no CAZy link return 'NA'
    return  # list, first CAZy class and then protein function


def create_summary_chart(CAZy_fam_column, logger):
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
