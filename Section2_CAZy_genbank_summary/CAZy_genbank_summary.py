#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Create summary of annotated CAZy classes in GenBank files.

:cmd_args:...

:func main:...

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
        description="Generate summary dataframes and bar charts of CAZy annotation in GenBank files",
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
    check_input(args.df_input, args.genbank, logger)
    logger.info("Inputs accepted")

    # If specified output directory for genomic files, create output directory
    if args.output is not sys.stdout:
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


def check_input(dataframe, genbank, logger):
    """Summary...

    Process...
    
    :param...

    Return...
    """
    logger.info("Checking path to input dataframe is valid")
    if dataframe.is_file() is False:
        # report to user and exit programme
        logger.info(
            (
                "Input dataframe not found. Check filename, extension and directory are correct."
                "\nTerminating program."
            ),
            exc_info=1,
        )
        sys.exit(1)
    logger.info("Checking path to GenBank file containing directory is valid")
    if genbank.exists is False:
        logger.info(
            (
                "GenBank file directory not found. Check correct directory was provided."
                "\nTerminating program."
            ),
            exc_info=1,
        )
        sys.exit(1)
    
    return


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
            )
            sys.exit()
    return


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
    if accession_number == 'NA':
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

    # check if file exists
    if file_term.exists():  # alter so if does not exist, look up exists
        # log warning
        return('NA')

    return(protein_name, protein_id)


def get_cazy_family(protein_id, logger):
    # use protein_id to call to UniProt
    # if CAZy link return class - may return full family rather than class - no worries
    return  # return family as str


def create_summary_chart(CAZy_fam_column, logger):
    # Use seaborn to create summary chart
    return


def write_out_chart(chart, output, logger):
    # write out chart to specified directory
    return


def write_out_df(df, output, logger):
    # write out df to specified directory
    return


if __name__ == "__main__":
    main()
