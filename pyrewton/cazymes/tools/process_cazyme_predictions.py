#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author:
# Emma E. M. Hobbs
#
# Contact
# eemh1@st-andrews.ac.uk
#
# Emma E. M. Hobbs,
# Biomolecular Sciences Building,
# University of St Andrews,
# North Haugh Campus,
# St Andrews,
# KY16 9ST
# Scotland,
# UK
#
# The MIT License
"""Process output from dbCAN, CUPP and eCAMI into dataframes.

:cmd_args input: (Required) path to directory containing prediction tool outputs
:cmd_args tool: (Required) 'dbcan', 'cupp', 'ecami', define which tools output is to be processed
:cmd_args force: Force enabling of writing existing output directory
:cmd_args log: Enable writing of log file
:cmd_args nodelete: Enable not deleting content in existing output directory
:cmd_args output: Define output directory to which output dataframes are written
:cmd_args verbose: Change logger level from WARNING to INFO - verbose logging

:func main: Build parser, build logger and coordinate script operation
:func parse_inputs: parse input dataframe and retrieve files/dirs in cwd
:func write_dbcan_dfs: write df for dbcan, HMMER, DIAMOND and Hotpep
:func parse_dbcan_overview_file: parse 'overview.txt' output from dbCAN
:func standardise_dbcan_results: remove added data from predicated CAZy family
:func get_dbcan_consensus: get the consensus results from dbCAN
:func write_cupp_df: write CUPP results in dataframe
:func parse_cupp_output: parse output from CUPP, to create dataframe data
:func write_ecami_df: write out dataframe of eCAMI results
:func parse_ecami_output: parse eCAMI output file, to create dataframe data

Write out dataframe containing gene IDs and result from the given CAZyme
prediction tool. The output from one CAZyme prediction tool can be processed
at a time.
"""

import logging
import re

import pandas as pd

from os import path
from pathlib import Path
from typing import List, Optional

from pyrewton.file_io import make_output_directory, write_out_pre_named_dataframe
from pyrewton.loggers import build_logger
from pyrewton.parsers.parser_process_cazyme_predictions import build_parser


def main(argv: Optional[List[str]] = None, logger: Optional[logging.Logger] = None):
    """Build parser and logger, and coordinate creation of dataframes."""

    # Programme preparation
    # Build parser
    # Check if namepsace isn't passed, if not parse command-line
    if argv is None:
        # Parse command-line
        args = build_parser().parse_args()
    else:
        args = build_parser(argv).parse_args()

    # Build logger
    # Note: log file only created if specified at cmdline
    if logger is None:
        logger = build_logger("process_cazyme_predictions", args)
    logger.info("Run initated")

    # Make output directory to store all create dataframes
    make_output_directory(args.output, logger, args.force, args.nodelete)

    # Retrieve accession numbers of genomic assemblies to identify the tools' output directories
    genomic_accessions, path_list = parse_input(args, logger)

    # Build dataframe from dbCAN output
    # Create consensus dataframe and dataframe per tool
    if args.tool == "dbcan":
        write_dbcan_dfs(genomic_accessions, path_list, args, logger)

    # Build dataframe from CUPP output
    if args.cupp == "cupp":
        write_cupp_df(genomic_accessions, path_list, args, logger)

    # Build dataframe from eCAMI output
    if args.ecami == "ecami":
        write_ecami_df(genomic_accessions, path_list, args, logger)


def parse_input(args, logger):
    """Parse input csv file to retrieve accession numbers and retrieve dir or file listing.

    :param args: parser object
    :param logger: logger object

    Return list of accession numbers and list of directories in cwd.
    """
    df = pd.read_cvs(args.accessions, header=0, index_col=0)
    accessions = df["NCBI Accession Numbers"]

    # Generate list of directories or files in specified input parent directory
    path_list = []  # list to store directory names
    parent_dir = Path(args.input)  # get path to parent directory of tools outputs

    if args.tools == "dbcan":
        for entry in parent_dir.iterdir():
            # add the path for each DIRECTORY in the parent directory to the directory list
            if path.isdir(entry):
                path_list.append(entry)

    elif args.tools == "cupp":
        for entry in parent_dir.iterdir():
            # add the path for each .fasta.log FILE in the parent directory to the path list
            if path.isfile(entry) and str(entry).endswith(".fasta.log"):
                path_list.append(entry)

    elif args.tools == "ecami":
        for entry in parent_dir.iterdir():
            # add the path for each .txt file in the parent directory to the path list
            if path.isfile(entry) and str(entry).endswith(".txt"):
                path_list.append(entry)

    return accessions, path_list


def write_dbcan_dfs(accession_numbers, directories, args, logger):
    """Build dataframe from dbCAN output file.

    Create summary dataframe of dbCAN consensus result and per tool.
    Write out dataframes to csv files.

    :param accession_numbers: list, list of genomic assembly accession numbers
    :param directories: list, list of directory paths in args.input
    :param args: parser object
    :param logger: logger object

    Return nothing.
    """
    for accession in accession_numbers:
        # format accession to match format in dir name
        accession.replace(".", "_")
        for directory in directories:
            if directory.find(accession):
                # Write path to the 'overview.txt' in the directory containing the working
                # accession number
                overview_file_path = directory / "overview.txt"

                # Create datframes containing results from the overview.txt file
                # in order of dbcan_df, diamond_df, hmmer_df, hotpep_df
                dataframes = parse_dbcan_overview_file(overview_file_path, logger)

                # Remove k-mer cluster labeling from Hotpep results
                dataframes[3].apply(
                    standardise_dbcan_results, args=("Hotpep", logger), axis=1
                )

                # Remove standardise HMMER predicated CAZy class/family formating
                dataframes[2].apply(
                    standardise_dbcan_results, args=("HMMER", logger), axis=1
                )
                dataframes[0].apply(
                    standardise_dbcan_results,
                    args=("dbCAN consensus CAZyme prediction", logger),
                    axis=1,
                )

                # Write out dataframes to csv files
                dataframe_names = ["dbCAN", "DIAMOND", "HMMER", "Hotpep"]
                index = 0
                for index in range(len(dataframes)):
                    write_out_pre_named_dataframe(
                        dataframes[index],
                        f"{dataframe_names[index]}_{accession}_output.csv",
                        logger,
                        args.output,
                        args.force,
                    )
                    index += 1

    return


def parse_dbcan_overview_file(overview_file_path, logger):
    """Read 'overview.txt' from dbCAN, and produce dataframe.

    :param args: parser object
    :param logger: logger object

    Return list of dataframes.
    """
    with open(overview_file_path) as fh:
        file_lines = fh.read().splitlines()

    df_rows = []
    for line in file_lines:
        df_rows.append(line.split())

    df_data = df_rows[1:]  # exclude header row
    df = pd.DataFrame(
        df_data, columns=["Gene ID", "HMMER", "Hotpep", "DIAMOND", "#ofTools"]
    )

    # Create separate dataframes for DIAMOND, HMMER and Hotpep
    diamond_df = df[["Gene ID", "DIAMOND CAZyme prediction"]]
    hmmer_df = df[["Gene ID", "HMMER CAZyme prediction"]]
    hotpep_df = df[["Gene ID, Hotpep CAZyme prediction"]]

    # Create consensus dbCAN dataframe (result when 2 or more tools match predications)
    # rename "#ofTools" column to facilitate use of .contains
    dbcan_df = get_dbcan_consensus(df, logger)

    return [dbcan_df, diamond_df, hmmer_df, hotpep_df]


def standardise_dbcan_results(df_row, column_name, logger):
    """Format results to only include predicated CAZy family.

    Remove numbers in brackets.
    Remove EC numbers.
    Remove non-standardised CAZy family labels.

    :param df: pandas dataframe
    :param column_name: str, name of the column containing results
    :param logger: logger object

    Return dataframe.
    """
    # Remove content within parentheses
    column_content = df_row[f"{column_name}"].split("+")

    index = 0
    for index in range(len(column_content)):
        # Remove bracketed numbers
        column_content[index] = re.sub(r"\(.*\)", "", column_content[index])
        # Remove non-standardised class names, e.g "_Chitin_synth_1"
        column_content[index] = re.sub(r"_\D*", "", column_content[index])
        index += 1
    new_column_content = "+".join(column_content)
    df_row[f"{column_name}"] = new_column_content

    return df_row


def get_dbcan_consensus(df, logger):
    """Get dataframe of consensus dbCAN results.

    :param df: pandas dataframe, dbCAN output formated into dataframe
    :param logger: logger object

    Return dataframe of consensus dbCAN output, where 2 or more tools match.
    """
    # Retrieve the rows from the dataframe that have 2 or 3 results
    # recorded in the '#ofTools' column
    df = df.rename(columns={"#ofTools": "Matches"})
    consensus_found = df.Matches.str.contains("2" or "3")
    df = df[consensus_found]

    # Create consensus result to summarise all
    consensus_result = []
    index = 0
    for index in range(len(df)):
        df_row = df.iloc[index]

        row_data = []  # empty list to store data for dbCAN df row

        if df_row["HMMER"] != "-":
            row_data.append(df_row["HMMER"])
        elif df_row["Hotpep"] != "-":
            row_data.append(df_row["Hotpep"])
        else:
            row_data.append(df_row["DIAMOND"])
        consensus_result.append(row_data)
        index += 1

    dbcan_df = df["Gene ID"]
    dbcan_df["dbCAN consensus CAZyme prediction"] = consensus_result

    return dbcan_df


def write_cupp_df(accession_numbers, files, args, logger):
    """Formate output from CUPP into Pandas dataframe.

    :param accession_numbers: list, list of genomic assembly accession numbers
    :param directories: list, list of directory paths in args.input
    :param args: parser object
    :param logger: logger object

    Return nothing.
    """
    for accession in accession_numbers:
        # format accession to match format in dir name
        accession.replace(".", "_")
        for entry in files:
            if entry.find(accession):
                # Collect predicated CAZyme families from output
                # Store in tuple with each list containing the predicated CAZyme
                # family of a unique protein
                dataframe_data = parse_cupp_output(entry, logger)

                df = pd.DataFrame(
                    dataframe_data, columns=["Gene ID", "CUPP CAZyme prediction"]
                )

                write_out_pre_named_dataframe(
                    df, "0HERER0", logger, args.output, args.force
                )

    return


def parse_cupp_output(output_file, logger):
    """Parse CUPP output file, creating tuple of predicated CAZymes.

    Each list in the tuple contains two items, the first item being
    the gene ID, and the second the predicated CAZy family.

    :param output_file: path, path to CUPP output file
    :param logger: logger object

    Return tuple.
    """
    dataframe_data = []  # empty tuple to store dataframe data
    with open(output_file) as fh:
        file_lines = fh.read().splitlines
    for line in file_lines:
        line_data = []  # list to store CAZyme family prediction from working line
        line = line.split()
        line_data.append(line[0])  # add gene ID
        line_data.append(line[2])  # add CAZy family prediction
        dataframe_data.append(line_data)

    return dataframe_data


def write_ecami_df(accession_numbers, files, args, logger):
    """Retrieve data from CUPP output and store in dataframe.

    :param accession_numbers: list, list of genomic assembly accession numbers
    :param directories: list, list of directory paths in args.input
    :param args: parser object
    :param logger: logger object

    Return nothing.
    """
    for accession in accession_numbers:
        # format accession to match format in dir name
        accession.replace(".", "_")
        for entry in files:
            if entry.find(accession):
                # Collect predicated CAZyme families from output
                # Store in tuple with each list containing the predicated CAZyme
                # family of a unique protein
                dataframe_data = parse_ecami_output(entry, logger)

                # Construct dataframe
                ecami_df = pd.DataFrame(dataframe_data, columns=["Gene ID", "eCAMI"])

                # Write out dataframe containing all predicated CAZymes for genomic accession
                write_out_pre_named_dataframe(
                    ecami_df, "0HERE0", logger, args.output, args.force
                )

    return


def parse_ecami_output(output_file, logger):
    """Parse eCAMI output file, creating tuple of predicated CAZymes.

    Each list in the tuple contains two items, the first item being
    the gene ID, and the second the predicated CAZy family.

    :param output_file: path, path to eCAMI output file
    :param logger: logger object

    Return tuple.
    """
    dataframe_data = []  # empty tuple to store dataframe data
    with open(output_file) as fh:
        file_lines = fh.read().splitlines
    for line in file_lines:
        line_data = []  # list to store CAZyme family prediction from working line
        if line.startswith(">"):
            # Retrieve the gene ID
            search_result = re.search(r">\D+\d*?\.\d\|", line)
            gene_id = search_result.group()
            gene_id = re.sub(r">|\|", "", gene_id)

            # Remove gene ID from current working line
            line = line.replace(search_result.group(), "")

            # Retrieve predicated CAZy families from working line
            search_result = re.search(r"\D+.+?\|(\d+| )", line)
            line = search_result.group()

            # Add formated gene ID and predicated CAZyme family to line data
            line_data.append(gene_id)
            line_data.append(line[:-2].replace("|", "+"))
            dataframe_data.append(line_data)

    return dataframe_data


if __name__ == "__main__":
    main()
