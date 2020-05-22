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
:func get_dataframe: build dataframe summarising CAZy annotation in GenBank files
:func create_df_foundation: Parse input dataframe row
:func build_df_foundation: Compile row data for dataframe
:func get_protein_data: Retrieve protein name and IDs

Generate summary dataframe and of annotated CAZy classes in all GenBank
files associated with a given species.

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
import logging
import shutil
import sys
from pathlib import Path

import pandas as pd
import seaborn as sns
from Bio import SeqIO
from bioservices import UniProt


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
        type=bool,
        metavar="force file overwritting",
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
        type=bool,
        metavar="not to delete exisiting files",
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
    CAZy_summary_df = get_dataframe(input_df, args, logger)

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
                shutil.rmtree(output)
            # 'Nodelete' enabled, don't empty directory
            else:
                logger.info(
                    "Output directory already exists and forced addition of files"
                    "to outdir enabled."
                )
    # Recursively make output directory
    try:
        output.mkdir(exist_ok=True)
    except OSError:
        # this will occur if directory already exists
        # ignored if forced over write enabled
        if force is True:
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
        names=["Genus", "Species", "NCBI Taxonomy ID", "NCBI Accession Numbers"]
    )
    return input_df


def get_dataframe(input_df, args, logger):
    """Build datafame.

    :param input_df:
    :param args:
    :param logger:

    Return dataframe.
    """
    # Build dataframe foundation: genus, species, accession number, protein name, protein ID
    all_foundation_data = []  # empty list to store all data for foundation dataframe
    all_foundation_data.append(input_df.apply(
        lambda df_row: create_df_foundation(df_row, args, logger),
        axis=1,
        )
    )
    CAZy_summary_df = pd.Dataframe(all_foundation_data, columns=[
        "Genus",
        "Species",
        "NCBI Accession Number",
        "Protein Name",
        "Protein ID"
        ])

    # Add CAZy data to dataframe:
    # Vectorise over protein ID column to call to NCBI to see if CAZy linked, and return class,
    # stored in a new column 'CAZy class' and protein funciton 'Function'
    # if cazy class returned full section will be titled 'cazy class',
    # if familied returned use 'cazy family' instead
    CAZy_summary_df["Cazy Class", "Function"] = CAZy_summary_df.apply(
        lambda column: get_cazy_family(column["Protein ID"], logger), axis=1
    )


def create_df_foundation(df_row, args, logger):
    """Prepare row data to create dataframe.

    Parse row from input dataframe as pandas series.
    row[0] = Genus
    row[1] = Species
    row[2] = NCBI Taxonomy ID
    row[3] = NCBI Accession numbers

    Retrieve all row data for single species, as a tuple. Each row is represented
    as a unique list in tuple. Each row/list containing genus, species, accession
    number, protein name and protein ID, with a unique protein name and ID per
    row/list.

    :param df_row: row from input_df (dataframe)
    :param args: parser arguments
    :param logger: logger object

    Return tuple.
    """
    logger.info(
        (
            "Adding scientific name, accession numbers, protein names and IDs\n"
            f"to dataframe for {df_row[0][0]}.{df_row[1]}"
        )
    )
    # Separate accession numbers in string, to form a list of accession numbers
    # with a unique accession number per index number
    single_species_data = get_species_data(
        df_row[0],
        df_row[1],
        df_row[3].split(", "),
        args,
        logger,
    )


def get_species_data(genus, species, accession_list, args, logger):
    """Build dataframe of genus, species and unique accession number per row.

    :param genus: str, genus name
    :param species: str, species name
    :param accession_list: list, list of accession numbers
    :param args: parser arguments
    :param logger: logger object

    Return dataframe.
    """
    # Row data
    row_foundation = []  # genus, species and accession
    complete_row_data = []  # genus, speccies, accession, protein name, protein ID
    all_rows_data = []  # data for all rows

    for accession in tqdm(accession_list, desc="Compiling data"):
        row_foundation.append([genus, species, accession])
        protein_data = get_protein_data(accession, len(accession_list), args.genbank, logger)
        # construct data for completed row, then compile all dataframe data
        for name_id_pair in protein_data:
            # add protein name and ID as individual items so as not create list within a list
            all_rows_data.append(row_foundation.append(name_id_pair[0], name_id_pair[1]))

    return all_rows_data


def get_protein_data(accession_number, total_accession, genbank_input, logger):
    """Retrieve protein names and ID from GenBank file.

    From each record the protein name and ID will be retrieved, and stored as
    a list.

    Lists wil be added to a single tuple containing all protein data.

    :param accession_number: str
    :param total_accession: int, total number of accessions
    :param genbank_input: path, path to directory containing GenBank files
    :param logger: logger object

    Return tuple.
    """
    # use accession number to open associated genbank file
    # use SeqIO to extract protein ID and protein name
    # return as list for dataframe
    # will need to check if protein ID already present in dataframe, and is don't add

    # check if accession number was provided
    if accession_number == 'NA':
        logger.info(
            (
                "Null values ('NA') was contained in cell, exiting retrieval of protein data.\n"
                "Returning null ('NA') value."
            )
        )
        return ('NA')

    count = 1
    for count in range(total_accession):
        if genbank_input.glob(accession_number + '*genomic.gbff.gz'):
            logger.info(
                (
                    f"Opening GenBank file for {accession_number},"
                    f"accession {count} of {total_accession}"
                )
            )
            # Open and extract protein ID and protein name

            protein_name = a
            protein_id = a

        else:
            logger.info(
                (
                    f"Could not find GenBank file for {accession_number},"
                    "exiting retrieval of protein data.\nReturning null ('NA') value."
                )
            )
            return ('NA')

    return(protein_name, protein_id)

    # create search term for file # look up re again
    file_term = accession_number.replace('.', '_').re...

    # create path to genbank file
    genbank_path = genbank_input + '/' + file_term

    # check if file exists
    if file_term.exists():  # alter so if does not exist, look up exists
        # log warning
        return('NA')

    return(protein_name, protein_id)


def get_cazy_family(protein_id, logger):
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
