#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Create summary of annotated CAZy classes in GenBank files.

<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
=======
>>>>>>> update module docstring
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
<<<<<<< HEAD
<<<<<<< HEAD
:func create_dataframe: build dataframe summarising CAZy annotation in GenBank files
:func create_df_foundation: Parse input dataframe row
:func build_df_foundation: Compile row data for dataframe
:func get_genbank_protein_data: Retrieve protein name and IDs

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
import io
import logging
import shutil
import sys

=======
:func main:
=======
:cmd_args:...
>>>>>>> create parser function and update build logger

:func main:...
=======
:func get_dataframe: build dataframe summarising CAZy annotation in GenBank files
=======
:func create_dataframe: build dataframe summarising CAZy annotation in GenBank files
>>>>>>> add gb_file search and error logging
:func create_df_foundation: Parse input dataframe row
:func build_df_foundation: Compile row data for dataframe
:func get_protein_data: Retrieve protein name and IDs
>>>>>>> update module docstring

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
<<<<<<< HEAD
>>>>>>> First draft of script
=======

>>>>>>> remove lambda
from pathlib import Path

import pandas as pd
import seaborn as sns
<<<<<<< HEAD
<<<<<<< HEAD
=======
>>>>>>> remove lambda

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
=======
from Bio import SeqIO
from bioservices import UniProt


def build_parser():
    """Return ArgumentParser parser for script."""
    # Create parser object
    parser = argparse.ArgumentParser(
        prog="CAZy_genbank_summary.py",
        description="Generate summary dataframes and bar charts of CAZy annotation in GenBank files",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

<<<<<<< HEAD
    # Define cmd-line args
>>>>>>> First draft of script
=======
    # Add arguments to parser
>>>>>>> create parser function and update build logger
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
<<<<<<< HEAD
<<<<<<< HEAD
        dest="force",
        action="store_true",
=======
        type=bool,
        metavar="force file overwritting",
>>>>>>> First draft of script
=======
        dest="force",
        action="store_true",
>>>>>>> remove lambda
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
<<<<<<< HEAD
<<<<<<< HEAD
        dest="nodelete",
        action="store_true",
=======
        type=bool,
        metavar="not to delete exisiting files",
>>>>>>> First draft of script
=======
        dest="nodelete",
        action="store_true",
>>>>>>> remove lambda
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

<<<<<<< HEAD
<<<<<<< HEAD
=======
>>>>>>> create parser function and update build logger
    return parser


def main():
    """docstring summary.

    Detail.

    Return.
    """
    # Create parser object for cmd-line ctrl
    parser = build_parser()
<<<<<<< HEAD
=======
    # Parse arguments into args variable
>>>>>>> First draft of script
=======
>>>>>>> create parser function and update build logger
    args = parser.parse_args()

    # Initiate logger
    # Note: log file only created if specified at cmdline
<<<<<<< HEAD
    build_logger("cazy_genbank_summary", args.log)
    logger = logging.getLogger("cazy_genbank_summary")
    logger.info("Run initated")

    # Check inputs are valid
    check_input(args, logger)
=======
    build_logger("CAZy_genbank_summary", args.log)
    logger = logging.getLogger("CAZy_genbank_summary")
    logger.info("Run initated")

    # Check inputs are valid
    check_input(args.df_input, args.genbank, logger)
>>>>>>> First draft of script
    logger.info("Inputs accepted")

    # If specified output directory for genomic files, create output directory
    if args.output is not sys.stdout:
<<<<<<< HEAD
        make_output_directory(args, logger)

    # Open input dataframe
    logger.info("Opening input dataframe")
    input_df = get_input_df(args.df_input, logger)

    # Build dataframe
<<<<<<< HEAD
    cazy_summary_df = create_dataframe(input_df, args, logger)
=======
    CAZy_summary_df = create_dataframe(input_df, args, logger)
>>>>>>> add gb_file search and error logging

    # Create summary charts of CAZy annotation distribution

    # Write out dataframe

    return


def build_logger(script_name, log_file,) -> logging.Logger:
    """Return a logger for this script.

    Enables logger for script, sets parameters and creates new file to store log.

    :param script_name: str, name of script
    :param log_file: parser argument, enable writing out of log file
=======
        make_output_directory(args.output, logger, args.force, args.nodelete)

    # Open input dataframe
    logger.info("Opening input dataframe")
    input_df = pd.read_csv(
        args.df_input,
        header=0,
        names=["Genus", "Species", "NCBI Taxonomy ID", "NCBI Accession Numbers"],
    )  # is it necessary to save to a variable?

    # create CAZy summary per species in input dataframe
    # dataframe written to output in function not main().
    input_df.apply(
        lambda df_row: create_cazy_df(df_row, args.genbank, args.output, logger),
        axis=1,
    )
    return


def build_logger(script_name, log_file,) -> logging.Logger:
    """Return a logger for this script.

    Enables logger for script, sets parameters and creates new file to store log.
<<<<<<< HEAD
    
    :param script_name: Name of script
    :param custom_string: Additional string parsed from cmdline by user - required for log
                    file to be written out
>>>>>>> First draft of script
=======

    :param script_name: str, name of script
    :param log_file: parser argument, enable writing out of log file
>>>>>>> create parser function and update build logger

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


<<<<<<< HEAD
def check_input(args, logger):
    """Check paths to input dataframe and GenBank files is valid.

    :param args: parser arguments
    :param logger: logger object

    Return nothing if paths are valid.
    """
    logger.info("Checking path to input dataframe is valid")
    if (args.df_input).is_file() is False:
=======
def check_input(dataframe, genbank, logger):
    """Summary...

    Process...
    
    :param...

    Return...
    """
    logger.info("Checking path to input dataframe is valid")
    if dataframe.is_file() is False:
        # report to user and exit programme
>>>>>>> First draft of script
        logger.info(
            (
                "Input dataframe not found. Check filename, extension and directory are correct."
                "\nTerminating program."
            ),
            exc_info=1,
        )
        sys.exit(1)
<<<<<<< HEAD

    logger.info("Checking path to GenBank file containing directory is valid")
    if (args.genbank).exists is False:
=======
    logger.info("Checking path to GenBank file containing directory is valid")
    if genbank.exists is False:
>>>>>>> First draft of script
        logger.info(
            (
                "GenBank file directory not found. Check correct directory was provided."
                "\nTerminating program."
            ),
            exc_info=1,
        )
        sys.exit(1)
<<<<<<< HEAD

    return


def make_output_directory(args, logger):
=======
    
    return


def make_output_directory(output, logger, force, nodelete):
>>>>>>> First draft of script
    """Create output directory for genomic files.

    Check if directory indicated for output existed already.
    If so check if force overwrite enabled. If not terminate programme.
    If so, check if deletion of exiting files was enabled.
    If so, exiting files in output directory are deleted.
    Create output directory, expecting error if already exists.

<<<<<<< HEAD
    :param args: parser arguments
    :param logger: logger object
=======
    :param output: Path, path to output directory
    :param logger: logger object
    :param force: bool, cmd-line args to enable/disable over writing existing directory
    :param nodelete: bool, cmd-line args to enable/disable deleting of existing files
>>>>>>> First draft of script

    Return nothing.
    """
    logger.info("Checking if specified output directory for genomic files exists.")
<<<<<<< HEAD
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
=======
    # If output directory specificed at cmd-line, check output directory does not already exist
    if output.exists():
        if force is False:
            logger.info(
                "Output directory already exists and forced overwrite not enabled.\nTerminating program."
            )
            sys.exit()
        else:
            if nodelete is False:
                logger.info(
                    "Output directory already exists and forced complete overwrite enabled.\nDeleting existing content in outdir."
                )
                # delete existing content in outdir
                shutil.rmtree(output)
            else:
                logger.info(
                    "Output directory already exists and forced addition of files to outdir enables."
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
                "OSError occured while creating output directory for genomic files.\nTerminating programme."
>>>>>>> First draft of script
            )
            sys.exit()
    return


<<<<<<< HEAD
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

<<<<<<< HEAD
    -- Dev note: add explanation for for loops and not apply:
    httpS://ys-l.github.io/posts/2015/08/28/how-not-to-use-pandas-apply

    Iterate over input dataframe row wise. This allows for converting human
    readable accession number list from a string to a Python list.

    Per species, retrieve all protein names and IDs for every accession number,
=======
    Iterate over input dataframe row wise. This allows for converting human
    readable accession number list from a string to a Python list.

<<<<<<< HEAD
    Per species, retrieve all protein names and IDs for every accession number, 
>>>>>>> add gb_file search and error logging
=======
    Per species, retrieve all protein names and IDs for every accession number,
>>>>>>> remove lambda
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
<<<<<<< HEAD
<<<<<<< HEAD
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
    df_index = 0
    for df_index in range(len(input_df["Genus"])):
        cazy_summary_df = cazy_summary_df.append(
            get_df_foundation_data(input_df.iloc[df_index], args, logger),
            ignore_index=True,
        )
        df_index += 1

    print("=====\nFoundation dataframe:\n", cazy_summary_df)

    # Build empty dataframe to store data retrieved from UniProtKB
    # GO = Gene Ontology (GO) project
    all_uniprotkb_data_df = pd.DataFrame(
        columns=[
            "UniProt entry ID",
            "UniProt entry name",
            "UniProt assigned protein names",
            "Length (Aa)",
            "Mass (Da)",
            "Domains",
            "Domain count",
            "UniProt linked protien families",
            "GO IDs",
            "GO molecular function",
            "GO biologocial process",
        ]
=======
    # Retrieve data for dataframe foundation: genus, species, accession number, protein name, protein ID
=======
    # Retrieve data for dataframe foundation: genus, species, accession number, protein name,
    # protein ID
>>>>>>> remove lambda
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
>>>>>>> add gb_file search and error logging
    )

    # parse current cazy_summary_df, retreiving UniProtKB data for each protein
    df_index = 0
    for df_index in range(len(cazy_summary_df["Genus"])):
        all_uniprotkb_data_df = all_uniprotkb_data_df.append(
            get_uniprotkb_data(cazy_summary_df.iloc[df_index], logger),
            ignore_index=True,
        )
        df_index += 1

    print("=====\nUniProt dataframe:\n", all_uniprotkb_data_df)

<<<<<<< HEAD
    # join UniProt dataframe to foundation dataframe
    cazy_summary_df = pd.concat([cazy_summary_df, all_uniprotkb_data_df], axis=1)

    print("=====\nFoundatoin + UniProt dataframe:\n", all_uniprotkb_data_df)

    return cazy_summary_df


def get_df_foundation_data(df_row, args, logger):
    """Prepare row data to create dataframe.

    Coordinate retrieval of protein data from GenBank files for every accession
    number for the species passed to the function.

<<<<<<< HEAD
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
=======
def get_df_foundation_data(genus, species, accession_list, args, logger):
=======
def get_df_foundation_data(df_row, args, logger):
>>>>>>> remove lambda
    """Prepare row data to create dataframe.

    Retrieve all row data for single species, as a tuple. Each row is represented
    as a unique list in tuple. Each row/list containing genus, species, accession
    number, protein name and protein ID, with a unique protein name and ID per
    row/list.
>>>>>>> add gb_file search and error logging

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

    # convert human readable list of accession numbers into Python list
    accession_list = df_row[3].split(", ")

    # open GenBank file for each accession number in the list and retrieve
    # all protein data in that GenBank file, stored as a tuple with each
    # list in the tuple containing data for a unique protein
    for accession in accession_list:
        protein_data = get_genbank_protein_data(accession, args, logger)  # tuple

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
            protein_data_df = protein_data_df.append(new_df_row, ignore_index=True)

        # if data was returned add to dataframe, with unique protein per row
        else:
            # For each unique protein in the GenBank file create a new row in
            # dataframe. The data for each unique protein is stored as a single
            # list in the tuple protein_data
            protein_index = 0  # index number in protein_data tuple
            for protein_index in tqdm(
                range(len(protein_data)),
                desc=f"Getting proteins {df_row[2]}-{accession}",
            ):
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
                protein_data_df = protein_data_df.append(new_df_row, ignore_index=True)
                protein_index += 1

    return protein_data_df


def get_genbank_protein_data(accession_number, args, logger):
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
<<<<<<< HEAD
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

    # retrieve GenBank file for accession number
    gb_file = get_genbank_file(
        accession_number, args, logger
    )  # list with GenBank file with index [0]

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

    # create empty list to store protein data
    all_protein_data = []

    # Retrieve protein data from GenBank file
    with gzip.open(gb_file[0], "rt") as handle:
        # create list to store all protein data retrieved from GenBank file, making it a tuple
        for gb_record in SeqIO.parse(handle, "genbank"):
            for (index, feature) in enumerate(gb_record.features):
                # empty protein data list so as not to contaminate data of next protein
                protein_data = []
                # Parse over only protein encoding features (type = 'CDS')
                if feature.type == "CDS":
                    # extract protein ID
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
                        # if null value was returned for every feature attribute log error
                        # and don't add to all_protein_data list
                        if protein_data == ["NA", "NA", "NA", "NA"]:
                            logger.warning(
                                f"No data retrieved from CDS type feature, index: {index}",
                                exc_info=1,
                            )
                        # if some data retrieved, add to all_protein_list
                        else:
                            all_protein_data.append(protein_data)
                        # if some data was retrieved all to all_protein_data list
                        if protein_data != ["NA", "NA", "NA", "NA"]:
                            all_protein_data.append(protein_data)

                    else:
                        # error occured in that one of the appending actions failed to append
                        # and would lead to misalignment in the dataframe if added to the
                        # all_protein_data list
                        logger.warning(
                            (
                                f"Error occured during retrieval of data from feature, {index}\n"
                                f"for {accession_number}. Returning no protein data"
                            )
                        )

    return all_protein_data


def get_genbank_file(accession, args, logger):
    """Retrieve GenBank file for accession number in local dir.

    :param accession: str, accession number of GenBank file
    :param args: parser arguments
    :param logger: logger object

    Return list of length 1, containing path to GenBank file.
    """
    # replace '.' with '_' to match format in GenBank file name
    file_stem = accession.replace(".", "_")

    # create empty list to store file entries, to allow checking if multiple files were retrieved
    gb_file = []
=======
    logger.info(
        (
            "Adding scientific name, accession numbers, protein names and IDs\n"
            f"to dataframe for {df_row[0][0]}.{df_row[1]}, {df_row[2]}"
        )
    )

    # convert human readable list of accession numbers into Python list
    accession_list = df_row[3].split(", ")

    # create empty list to store all row_data lists, producing a tuple
    all_rows_data = []  # data for all rows

    # open GenBank file of each accession number in the list and retrieve
    # all protein data in that GenBank file, stored as a tuple with each
    # list in the tuple containing data for a unique protein
    for accession in tqdm(accession_list, desc="Compiling data"):
        protein_data = get_protein_data(accession, args.genbank, logger)

        # For each unique protein in the GenBank file, create a new row in
        # dataframe by compiling the protein's data into a single list
        # and adding the list to the all_rows_data tuple
        tuple_index = 0
        list_index = 0
        for tuple_index in range(len(protein_data)):
            # create empty list to store data for new dataframe row
            row_data = []

            # add genus, species, taxonomy ID and accession number to row_data
            row_data.append(df_row[0])
            row_data.append(df_row[1])
            row_data.append(df_row[2])
            row_data.append(df_row[3])

            # add protein data one item at a time, so they can populate different
            # columns in the data frame
            for list_index in range(len(protein_data[tuple_index])):
                row_data.append(protein_data[tuple_index][list_index])
                list_index += 1

            # add row_data to all_row_data tuple
            all_rows_data.append(row_data)
            tuple_index += 1

    #  adding to a list of row_data containing:
    # genus, species, accession number, protein ID, locus tag, location
    # and function

    # Which is then added to tuple of all row data, with each item being
    # a unique row containing a unique protein

    count = 1
    for accession in tqdm(accession_list, desc="Compiling data"):
        row_foundation.append([df_row[0], df_row[1], accession])
        protein_data = get_protein_data(accession, args.genbank, logger)
        # construct data for completed row, then compile all dataframe data
        index = 0
        for index in range(len(protein_data)):
            # add protein ID, locus tag, location and function as individual items so as not to create a tuple
            all_rows_data.append(
                row_foundation.append(
                    protein_data[index][0],
                    protein_data[index][1],
                    protein_data[index][2],
                    protein_data[index][3],
                )
            )
            index += 1
            count += 1
>>>>>>> add gb_file search and error logging

    # retrieve all files from directory
    files_in_entries = (
        entry for entry in Path(args.genbank).iterdir() if entry.is_file()
    )
    for item in files_in_entries:
        # search for accession number's GenBank file
        if item.name.startswith(f"{file_stem}") and item.name.endswith(".gbff.gz"):
            gb_file.append(item)

    return gb_file

<<<<<<< HEAD
=======
def get_protein_data(accession_number, genbank_input, logger):
<<<<<<< HEAD
    """Retrieve protein names and ID from GenBank file.
>>>>>>> add gb_file search and error logging

def get_record_feature(feature, qualifier, logger):
    """Retrieve data from GenBank record feature.
=======
    """Retrieve protein ID, locus tag and function from GenBank file.

<<<<<<< HEAD
    From each record the protein ID, locus tag and function will be retrieved,
    and stored as a list.
>>>>>>> add parsing of genbank files
=======
    From each record the protein ID, locus tag, location and annotated
    function is retrieved, and stored as a list.
>>>>>>> move feature data retrieval to separate function

    :param feature: feature object, GenBank file record feature
    :param qualifier: str, key of feature attribute
    :param logger: logger object

<<<<<<< HEAD
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


def get_uniprotkb_data(df_row, logger):
    """Coordinate retrieval of protein data from UniProtKB.

    Reminder of structure of df_row (row from cazy_summary_df) when
    retrieving specific data
    df_row[0] = genus
    df_row[1] = species
    df_row[5] = protein unique locus tag

    :param df_row: pandas series, row from cazy_summary_df dataframe
=======
    Any failed to retrieve data will be returned as pandas null value 'NA'.

    :param accession_number: str
    :param genbank_input: path, path to directory containing GenBank files
>>>>>>> add gb_file search and error logging
    :param logger: logger object

    Return dataframe.
    """
    # create query term
    query = f'{df_row[5]} AND organism:"{df_row[0]} {df_row[1]}"'

    # establish data to be retrieved from UniProt
    columnlist = (
        "id,entry name,protein names,length,mass,domains,domain,"
        "families,"
        "go-id,go(molecular function),go(biological process)"
    )

    # open connection to UniProt() and convert result into pandas df
    search_result_df = pd.read_table(
        io.StringIO(UniProt().search(query, columns=columnlist))
    )

    # rename columns to match to indicate UniProtKB source of data
    search_result_df.rename(
        columns={
            "Entry": "UniProtKB Entry ID",
            "Entry name": "UniProtKB Entry Name",
            "Protein names": "UniProtKB Protein Names",
            "Length": "Length (Aa)",
            "Mass": "Mass (Da)",
            "Protein families": "UniProtKB Linked Protein Families",
        }
    )

    return search_result_df


def get_cazy_data(protein_id, logger):
    # use protein_id to call to UniProt
    # if CAZy link return class - may return full family rather than class - no worries
    # if no CAZy link return 'NA'
    return  # list, first CAZy class and then protein function


def create_summary_chart(cazy_fam_column, logger):
=======
def create_cazy_df(df_row, genbank_input, output, logger):
    """Create dataframe of annotated proteins in GenBank files.
    
    Parse row from dataframe as pandas series.
    row[0] = Genus
    row[1] = Species
    row[2] = NCBI Taxonomy ID
    row[3] = NCBI Accession numbers

    First letter of genus, species and taxonomy ID form name of
    dataframe summarising current CAZy annotation for species.

    :param df_row: row from input_df (dataframe)
    :param genbank_input:
    :param output:
    :param logger:

    Return ...
    """
    # Create dataframe. Name formate: 'abbreviated_scientific_name_taxID_CAZy_summary_df'
    # with one column: 'Accession Number', with a unique accession number per row
    logger.info(f"Creating dataframe summarising CAZy annotations in GenBank,\n     for {df_row[0]} {df_row[1]}, {df_row[2]}")
    logger.info("Retrieving scientific name and taxonomy ID")
    df_name = df_name = (
        df_row[0][0] + "_" + df_row[1] + "_" + df_row[2][9:] + "_CAZy_summary_df"
    )

    df_name = pd.DataFrame(df_row[3].split(", "), columns=["Accession Number"])

    # vectorise over accession_number column, using get_protein_data to fill out two new columns:
    # 'Protein ID' and 'Protein Name'
    logger.info(
        f"Retrieving protein names and IDs from annotations for {df_row[0]} {df_row[1]}, {df_row[2]}"
    )
    df_name["Protein Name", "Protein ID"] = df_name.apply(
        lambda column: get_protein_data(
            column["Accession Number"], len(df_name["Accession Number"]), genbank_input, logger
        ),
        axis=1,
    )

    # vectorise over protein ID column to call to NCBI to see if CAZy linked, and return class,
    # stored in a new column 'CAZy class'
    # if cazy class returned full section will be titled 'cazy class',
    # if familied returned use 'cazy family' instead
    df_name["Cazy Class"] = df_name.apply(
        lambda column: get_cazy_family(column["Protein ID"], logger), axis=1
    )

    # create bar chart summarise cazy class distribution/annotation frequency
    # if families were returned in the function called us the startswith
    # function to determine the family
    # check what input is required for seaborn
    chart = create_summary_chart(df_name["CAZy Family"], logger)

    # if enabled write out chart to output directory
    write_out_chart(chart, output, logger)

    # if enabled write out df to output directory
    write_out_df(df_name, output, logger)
    return


def get_protein_data(accession_number, total_accession, genbank_input, logger):
    """Summary...

    Process...

    args...

    Return...
    """
    # use accession number to open associated genbank file
    # use SeqIO to extract protein ID and protein name
    # return as list for dataframe
    # will need to check if protein ID already present in dataframe, and is don't add

    # check if accession number was provided
<<<<<<< HEAD
    if accession_number == 'NA':
<<<<<<< HEAD
        logger.info("Null values ('NA') was contained in cell, exiting retrieval of protein data.\nReturning null ('NA') value.")
        return ('NA')

    count = 1
    for count in range(total_accession):
        if genbank_input.glob(accession_number + '*genomic.gbff.gz'):
            logger.info(f"Opening GenBank file for {accession_number}, accession {count} of {total_accession}")
            # Open and extract protein ID and protein name

            protein_name = a
            protein_id = a
            
        else:
            logger.info(f"Could not find GenBank file for {accession_number}, exiting retrieval of protein data.\nReturning null ('NA') value.")
            return ('NA')
            
    return(protein_name,protein_id)



    # create search term for file # look up re again
    file_term = accession_number.replace('.', '_').re...

    
    # create path to genbank file
    genbank_path = genbank_input + '/' + file_term
=======
=======
    if accession_number == "NA":
<<<<<<< HEAD
>>>>>>> remove lambda
        logger.info(
=======
        logger.warning(
>>>>>>> add parsing of genbank files
            (
                f"Null value ('NA') was contained in cell for {accession_number},"
                "exiting retrieval of protein data.\nReturning null ('NA') value"
                "for all protein data"
            )
        )
        return ["NA", "NA", "NA", "NA"]

    # Retrieve GenBank (gb) file
    gb_file = list(Path(genbank_input).glob(rf"{accession_number}*.gbff.fz"))

    # check file was retrieved, not multiple or none
    if len(gb_file) > 1 or len(gb_file) == 0:
        # log error and return 'NA' for protein name and protein ID
        logger.warning(
            (
                f"Failed to retrieve GenBank file for {accession_number}.\n"
                "Returning null ('NA') value for all protein data"
            )
        )
<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
        return (['NA', 'NA'])
>>>>>>> add gb_file search and error logging
=======
        return ["NA", "NA"]
>>>>>>> remove lambda
=======
        return ["NA", "NA", "NA"]
>>>>>>> add parsing of genbank files
=======
        return ["NA", "NA", "NA", "NA"]
>>>>>>> add data for unique protein to df row

    else:
        all_protein_data = []
        # Retrieve protein data
        with gzip.open(gb_file[0], "rt") as handle:
            # create list to store all protein data retrieved from GenBank file, making it a tuple
            for gb_record in SeqIO.parse(handle, "genbank"):
                for (index, feature) in enumerate(gb_record.features):
                    # empty protein data list so as not to contaminate data of next protein
                    protein_data = []
                    # Parse over only protein encoding features (type = 'CDS')
                    if feature.type == "CDS":
                        # extract protein ID
                        protein_data.append(
                            get_record_feature(feature, "protein_id", logger)
                        )
                        # extract locus tag
                        protein_data.append(
                            get_record_feature(feature, "locus_tag", logger)
                        )
                        # extract location
                        protein_data.append(
                            get_record_feature(feature, "location", logger)
                        )
                        # extract annotated function of product
                        protein_data.append(
                            get_record_feature(feature, "product", logger)
                        )

                        # add protein data to total protein data list, only if data was retrieved
                        if protein_data != ["NA", "NA", "NA", "NA"]:
                            all_protein_data.append(protein_data)
                        else:
                            logger.warning(
                                f"No data retrieved from CDS type feature, index: {index}",
                                exc_info=1,
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
    return  # return family as str


def create_summary_chart(CAZy_fam_column, logger):
>>>>>>> First draft of script
    # Use seaborn to create summary chart
    return


<<<<<<< HEAD
def write_out_df(df, output, logger):
    # write out df to specified directory
    return


def write_out_chart(chart, output, logger):
    # write out chart to specified directory
=======
def write_out_chart(chart, output, logger):
    # write out chart to specified directory
    return


def write_out_df(df, output, logger):
    # write out df to specified directory
>>>>>>> First draft of script
    return


if __name__ == "__main__":
    main()
