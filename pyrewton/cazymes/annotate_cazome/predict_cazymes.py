#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) University of St Andrews 2022
# (c) University of Strathclyde 2022
# (c) James Hutton Institute 2022
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
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
"""Automate running CAZyme classifiers to predict CAZymes"""


import logging
import os
import re
import sys

from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from tqdm import tqdm
from typing import List, Optional

from pyrewton.cazymes.annotate_cazome.tools.parse_tools import (
    parse_cupp_output,
    parse_dbcan_output,
    parse_ecami_output
)
from pyrewton.cazymes.annotate_cazome.tools import invoke_prediction_tools
from pyrewton.utilities import build_logger
from pyrewton.utilities.file_io import make_output_directory, write_out_pre_named_dataframe


@dataclass
class Query:
    """Data of prediction tool query.

    :param fasta: path, path to query fasta file
    :param tax_id: str, NCBI taxonomy ID of host species
    :param source: str, source of proteins within fasta file
    :param prediction_dir: path, path to which prediction tool output is written
    """

    fasta: Path  # path to FASTA file containing proteins for prediction
    tax_id: str  # NCBI taxonomy id, prefix NCBI:txid
    source: str  # source of protein sequences, genomic assembly or database
    prediction_dir: Path  # path to dir containing outputs from all prediciton tools

    def __str__(self):
        """Representation of object"""
        return(
            (
                f"<tax-id={self.tax_id} source={self.source} "
                "fasta={self.fasta} prediction_dir={self.prediction_dir}"
            )
        )


def main(args: Optional[List[str]] = None, logger: Optional[logging.Logger] = None):
    """Set up parser and logger, and coordinate prediction of CAZymes."""
    # Initiate logger
    # Note: log file only created if specified at cmdline
    if logger is None:
        logger = build_logger("predict_cazymes", args)

    # check current working directory, to make sure can access the CAZyme prediction tools
    check_cwd(logger)

    # create paths to tools
    tool_paths = {}  # {tool : path}
    for toolname in args.tools:

        tool_dir = args.tool_dir / toolname

        if (tool_dir.exists()) and (tool_dir.is_dir()):
            # get absolute path
            tool_paths[toolname.lower()] = tool_dir.resolve()

        else:
            # user may have included capitals
            # try lower case
            
            tool_dir = args.tool_dir / toolname.lower()
            if (tool_dir.exists()) and (tool_dir.is_dir()):
                tool_paths[toolname.lower()] = tool_dir.resolve()

            else:
                logger.error(
                    f"Could not find directory for {toolname} in {args.tool_dir}\n"
                    f"Tried: {args.tool_dir / toolname}\n"
                    f"And: {args.tool_dir / toolname.lower()}\n"
                    "Check the correct directory has been provided.\n"
                    "Terminating program as can't find tool directories"
                )
                sys.exit(1)

    # If specified output directory for genomic files, create output directory
    if args.output is not sys.stdout:
        make_output_directory(args.output, logger, args.force, args.nodelete)

    # invoke prediction tools and build prediciton Query instances
    predictions = get_predictions(tool_paths, args, logger)

    # standardist output, statistical evaluation performance (if CAZy data provided), and write
    # summary report
    index = 0
    for index in tqdm(range(len(predictions)), desc="Standardising tools outputs"):
        write_prediction_report(predictions[index], args, logger)

    logger.info("Program finished, and no terminating.")


def check_cwd(logger):
    """Check user invoked script in correct directory so can access the prediciton tools.

    If user is not within any of pyrewton/cazymes/prediction module directories then
    raise error and terminate program. Excepted directories for starting at:
    pyrewton/cazymes/prediction
    pyrewton/cazymes/tools

    :param logger: logger object

    Return nothing.
    """
    current_path = os.getcwd()

    if current_path.endswith("pyrewton/cazymes/prediction"):
        return

    elif current_path.endswith("pyrewton/cazymes/prediction/tools"):
        logger.warning(
            (
                "Script invoked in 'tools' directory.\n"
                "Changed current working directory to pyrewton/cazymes/prediction."
            )
        )
        os.chdir('..')
        return

    else:
        logger.error(
            (
                "Script invokved in wrong directory.\n"
                "To be able to access dbCAN, CUPP and eCAMI please invokve script when cwd is\n"
                "pyrewton/cazymes/prediction within the pyrewton programm. Terminating program."
            )
        )
        sys.exit(1)


def get_predictions(tool_paths, args, logger):
    """Build prediction queries and invoke prediction tools for each query.

    :param tool_paths: dict {tool name : path to dir containing tool's db files}
    :param args: parser object
    :param logger: logger object

    Return list of Query class objects, queries to prediction tools.
    """
    # create list of paths to all fasta files in input directory
    all_fasta_paths = get_fasta_paths(args, logger)

    # create empty list to store all instances of Prediction class objects
    predictions = []

    # for each FASTA file invoke dbCAN, CUPP and eCAMI
    for file_path in tqdm(all_fasta_paths, desc="Invoking tools for FASTA file"):  # make tqdm
        # retrieve data on source of protein sequences and species taxonomy ID
        protein_source = get_protein_source(file_path, logger)
        tax_id = get_tax_id(file_path, logger)

        time_stamp = datetime.now().strftime("%Y-%m-%d--%H-%M-%S")

        # name output dir to store prediction output and statistical evaluation
        if tax_id is None:
            outdir_name = f"cazyme_predictions_{protein_source}_{time_stamp}"
        else:
            outdir_name = f"cazyme_predictions_{tax_id}_{protein_source}_{time_stamp}"

        # create output_dir for given input FASTA file within user specified parent directory
        output_path = args.output / outdir_name
        make_output_directory(output_path, logger, args.force, args.nodelete)

        # create Query class object to store data on the query made to the prediction tools
        prediction_tool_query = Query(file_path, tax_id, protein_source, output_path)

        # invoke prediction tools and retrieve paths to the prediction tools outputs
        full_outdir_path = invoke_prediction_tools(prediction_tool_query, logger, tool_paths)

        # update outdir path with full path
        prediction_tool_query.prediction_dir = full_outdir_path

        predictions.append(prediction_tool_query)

    return predictions


def get_fasta_paths(args, logger):
    """Retrieve paths to call FASTA files in input dir.

    :param args: parser object
    :param logger: logger object

    Returns list of paths to fasta files.
    """
    # create empty list to store the file entries, to allow checking if no files returned
    fasta_file_paths = []

    # retrieve all files from input directory
    files_in_entries = (
        entry for entry in Path(args.input).iterdir() if entry.is_file()
    )
    # retrieve paths to fasta files
    for item in files_in_entries:
        # search for fasta file extension
        if item.name.endswith(".fasta") or item.name.endswith(".fa"):
            fasta_file_paths.append(item)

    # check files were retrieved from input directory
    if len(fasta_file_paths) == 0:
        logger.warning(
            (
                "No FASTA files retrieved from input directory.\n"
                "Check the path to the input directory is correct.\n"
                "Terminanting program."
            )
        )
        sys.exit(1)

    return fasta_file_paths


def get_protein_source(file_path, logger):
    """Retrieve source of protein sequences from FASTA file path.

    :param file_path: path, path to fasta file
    :param logger: logger object

    Return string, source of protein sequences.
    """
    # Attempt to retrieve genomic assembly from FASTA file path
    search_result = re.search(r"GC(A|F)_\d+?_\d+?", str(file_path), re.IGNORECASE)

    try:
        protein_source = search_result.group()
        return protein_source
    except AttributeError:
        # search for uniprot query used to retrieve proteins
        search_result = re.search(r"uniprot(.+?\.)", str(file_path), re.IGNORECASE)
        try:
            protein_source = search_result.group()[:-1]  # remove terminal '.'
            return protein_source
        except AttributeError:
            protein_source = "unknown_source"
            return protein_source


def get_tax_id(file_path, logger):
    """Retrieve taxonomy ID from fasta file path.

    :param file_path: path, path to fasta file path
    :param logger: logger object

    Return string, taxonomy ID of proteins' host species.
    """
    # search for first taxonomy ID format
    search_result = re.search(r"txid\d+?\D", str(file_path), re.IGNORECASE)

    try:
        tax_id = search_result.group()[:-1]
        return tax_id
    except AttributeError:      
        # search for other taxonomy ID format
        search_result = re.search(r"taxonomy__\d+?__", str(file_path), re.IGNORECASE)
        try:
            tax_id = search_result.group()[:-2]
            return tax_id
        except AttributeError:
            return


def write_prediction_report(prediction, args, logger):
    """Coordinate standardising prediction tools outputs and writing summary report.

    If the data from CAZy is provided a statistical evaluation of the prediction tools, by using
    the data from CAZy as the 'known truth', will be performed.

    :param prediction: Query class object
    :param args: cmd args parser
    :param logger: logger object


    Return nothing
    """
    # Standardise the output from each prediction tool
    dbcan_stnd_df, hmmer_stnd_df, hotpep_stnd_df, diamond_stnd_df, cupp_stnd_df, ecami_stnd_df = standardise_prediction_outputs(
        prediction.prediction_dir, args, logger
    )
    
    if (
        (dbcan_stnd_df is None) and
        (hmmer_stnd_df is None) and
        (hotpep_stnd_df is None) and
        (diamond_stnd_df is None) and
        (cupp_stnd_df is None) and
        (ecami_stnd_df is None)
    ):
        return  # error in retrieving output was logged in get_output_files()
    
    # Perform statistical evaluation of performance if CAZy data provided
    # if args.cazy is not None:
        # pass

    # Produce summary report of predictions (number of CAZymes, numbers per class etc.)

    return


def standardise_prediction_outputs(out_dir, args, logger):
    """Coordinate standardising prediction tools outputs.

    Produces a dataframe for eCAMI, CUPP, Hotpep, HMMER, DIAMOND and a consensus result for dbCAN
    where consensus means any result where at least 2 tools (out of Hotpep, HMMER and DIAMOND)
    match.

    :param out_dir: path to directory containing respective prediction query's output files
    :param args: cmd args parser
    :param logger: logger object

    Return 6 Pandas dataframes.
    """
    # retrieve the paths to prediction tool output files in the output directory for 'prediction'
    raw_output_files = get_output_files(out_dir, logger)  # returns dict

    if raw_output_files is None:
        return  None, None, None, None, None, None  # error was logged in get_output_files()

    # Standardise the output from each prediction tool
    dbcan_stnd_df, hmmer_stnd_df, hotpep_stnd_df, diamond_stnd_df = (
        parse_dbcan_output.parse_dbcan_output(raw_output_files["dbcan"], logger)
    )

    # retrieve predicated EC number for Hotpep
    if hotpep_stnd_df is not None:
        hotpep_stnd_df = parse_dbcan_output.add_hotpep_ec_predictions(
            raw_output_files["hotpep"],
            hotpep_stnd_df,
            logger,
        )

    # standardise output from CUPP and eCAMI
    cupp_stnd_df = parse_cupp_output.parse_cupp_output(raw_output_files["cupp_raw"], logger)

    ecami_stnd_df = parse_ecami_output.parse_ecami_output(raw_output_files["ecami_raw"], logger)

    # Write out dataframes to disk
    write_out_dataframes(dbcan_stnd_df, "dbcan_stnd_df", out_dir, args, logger)
    write_out_dataframes(hmmer_stnd_df, "hmmer_stnd_df", out_dir, args, logger)
    write_out_dataframes(hotpep_stnd_df, "hotpep_stnd_df", out_dir, args, logger)
    write_out_dataframes(diamond_stnd_df, "diamond_stnd_df", out_dir, args, logger)
    write_out_dataframes(cupp_stnd_df, "cupp_stnd_df", out_dir, args, logger)
    write_out_dataframes(ecami_stnd_df, "ecami_stnd_df", out_dir, args, logger)

    return dbcan_stnd_df, hmmer_stnd_df, hotpep_stnd_df, diamond_stnd_df, cupp_stnd_df, ecami_stnd_df


def write_out_dataframes(df, df_name, out_dir, args, logger):
    """Coordinate writing out standardise dataframes to disk.

    Performs check and logs if no dataframe was produced for a given CAZyme prediction tool.

    :param df: pandas dataframe
    :param df_name: path to location where dataframe is to be written (excludes .csv extension)
    :param args: args parser object
    :param logger: logger object

    Return nothing.
    """
    if df is None:
        logger.warning(
            "No standardised dataframe was produced for\n"
            f"{out_dir} / {df_name}"
        )
        return

    if len(df["cazy_family"]) == 0:
        logger.warning(
            "Empty standardised dataframe created for\n"
            f"{out_dir} / {df_name}"
        )

    write_out_pre_named_dataframe(df, df_name, logger, out_dir, args.force)
    return


def get_output_files(output_dir, logger):
    """Retrieve files to CAZyme prediction tool output files.

    :param output_dir: path to the directory containing prediction tool outputs
    :param logger: logger object

    Return dictionary, keyed by tool name and valued by path to respective output file
    """
    # create empty dictionary to store file to output files, keyed by prediciton tool name
    output_dict = {}

    prediction_tools = ["dbcan", "hotpep", "cupp", "ecami"]

    # retrieve all files in the output directory for the current working FASTA file
    files_in_outdir = (
        entry for entry in output_dir.iterdir() if entry.is_file()
    )

    try:
        if len(list(files_in_outdir)) == 0:
            logger.error(
                    "Did not retrieve any prediction tool output files in\n"
                    f"{output_dir}\n"
                    "Not producing and output for this dir"
            )
            return
    except AttributeError:  # raised if files_in_outdir is None
        logger.error(
                "Did not retrieve any prediction tool output files in\n"
                f"{output_dir}\n"
                "Not producing and output for this dir"
        )
        return

    # Retrieve paths to specific prediction tool's output files
    for entry in files_in_outdir:
        # retrieve dbCAN overview.txt file
        if entry.name.endswith("overview.txt"):
            output_dict["dbcan"] = entry

        # retrieve hotpep.out file
        elif entry.name.endswith("hotpep.out"):
            output_dict["hotpep"] = entry

        # retrieve cupp output file
        elif entry.name.endswith("cupp_output.fasta"):
            output_dict["cupp"] = entry

        # retrieve ecami output file
        elif entry.name.endswith("ecami_output.txt"):
            output_dict["ecami"] = entry

    # check output file for each prediction tool was retrieved
    for tool in prediction_tools:
        try:
            output_dict[tool]
        except KeyError:
            output_dict[tool] = "NA"

    return output_dict


if __name__ == "__main__":
    main()
