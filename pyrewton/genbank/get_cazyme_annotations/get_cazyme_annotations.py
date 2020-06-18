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
"""Create summary of annotated CAZy classes in GenBank files.

:cmd_args df_input: path, path to input dataframe
:cmd_args force: bool, force overwriting files in output directory
:cmd_args genbank: path, path to directory containing GenBank files
:cmd_args log: path, path to direct writing out log file
:cmd_args nodelete: not delete existing files in output directory
:cmd_args output: path, path to output directory

:func main: Coordinate calling of other functions
:func create_dataframe: build dataframe summarising CAZy annotation in GenBank files
:func create_df_foundation: Parse input dataframe row
:func build_df_foundation: Compile row data for dataframe
:func get_genbank_protein_data: Retrieve protein name and IDs

Generate summary dataframe and of annotated cazymess in all GenBank
files directly linked to a given species.
"""

import gzip
import io
import logging
import re
import sys

from typing import List, Optional

import pandas as pd

from Bio import SeqIO
from bioservices import UniProt
from pandas.errors import EmptyDataError
from tqdm import tqdm
from urllib.error import HTTPError

from pyrewton.file_io import make_output_directory, input_dir_get_cazyme_annotations
from pyrewton.loggers import build_logger
from pyrewton.parsers.parser_get_cazyme_annotations import build_parser


def main(argv: Optional[List[str]] = None, logger: Optional[logging.Logger] = None):
    """docstring summary.

    Detail.

    Return.
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
        logger = build_logger("get_cazyme_annotations", args)
    logger.info("Run initated")

    # If specified output directory, create output directory
    if args.output is not sys.stdout:
        try:
            make_output_directory(args, logger)
        except FileExistsError:
            logger.error("Output directory %s already exists (exiting)" % args.output)
            sys.exit(1)

    # Open input dataframe
    logger.info("Opening input dataframe %s", args.df_input)
    input_df = pd.read_csv(args.df_input, header=0, index_col=0)

    # Build dataframe
    cazy_summary_df = create_dataframe(input_df, args, logger)

    # Create summary charts of CAZy annotation distribution

    # Write out dataframe


def create_dataframe(input_df, args, logger):
    """Build datafame.

    -- Dev note: add explanation for for loops and not apply:
    httpS://ys-l.github.io/posts/2015/08/28/how-not-to-use-pandas-apply

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
            "Protein Sequence",
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

    # these are debugging purposes and will not be included in final version
    print("=====Foundation dataframe======\n", cazy_summary_df, "\n")
    cazy_summary_df.to_csv("foundation_dataframe.csv")

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
    )

    # parse current cazy_summary_df, retreiving UniProtKB data for each protein
    df_index = 0
    for df_index in tqdm(
        range(len(cazy_summary_df["Genus"])), desc="Retrieving UniProt data"
    ):
        all_uniprotkb_data_df = all_uniprotkb_data_df.append(
            get_uniprotkb_data(cazy_summary_df.iloc[df_index], logger),
            ignore_index=True,
        )
        df_index += 1

    print("=====\nUniProt dataframe:\n", all_uniprotkb_data_df)

    # join UniProt dataframe to foundation dataframe
    cazy_summary_df = pd.concat([cazy_summary_df, all_uniprotkb_data_df], axis=1)

    print("=====\nFoundatoin + UniProt dataframe:\n", all_uniprotkb_data_df)

    return cazy_summary_df


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
            "Protein Sequence",
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
                "Protein Sequence": "NA",
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
                    "Protein Sequence": protein_data[protein_index][4],
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
    # check if accession number was provided
    if accession_number == "NA":
        logger.warning(
            (
                f"Null value ('NA') was contained in cell for {accession_number},"
                "exiting retrieval of protein data.\nReturning null ('NA') value"
                "for all protein data"
            )
        )
        return ["NA", "NA", "NA", "NA", "NA"]

    # retrieve GenBank file for accession number
    gb_file = input_dir_get_cazyme_annotations.get_genbank_file(
        accession_number, args, logger
    )  # list with GenBank file with index [0]

    # If retrieving of GenBank file failed, return 'NA' for all protein data
    # for accession number
    if gb_file is None:
        # error logging performd in get_genbank_file()
        return ["NA", "NA", "NA", "NA", "NA"]

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
                    # extract protein sequence
                    protein_data.append(
                        get_record_feature(feature, "translation", logger)
                    )

                    # add protein data to total protein data list, only if data was retrieved
                    if len(protein_data) == 5:
                        # if null value was returned for every feature attribute log error
                        # and don't add to all_protein_data list
                        if protein_data == ["NA", "NA", "NA", "NA", "NA"]:
                            logger.warning(
                                f"No data retrieved from CDS type feature, index: {index}",
                                exc_info=1,
                            )
                        # if some data retrieved, add to all_protein_list
                        else:
                            all_protein_data.append(protein_data)
                        # if some data was retrieved all to all_protein_data list
                        if protein_data != ["NA", "NA", "NA", "NA", "NA"]:
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


def get_uniprotkb_data(df_row, logger):
    """Coordinate retrieval of protein data from UniProtKB.

    Reminder of structure of df_row (row from cazy_summary_df) when
    retrieving specific data
    df_row[0] = genus
    df_row[1] = species
    df_row[5] = protein unique locus tag

    :param df_row: pandas series, row from cazy_summary_df dataframe
    :param logger: logger object

    Return dataframe.
    """
    # establish data to be retrieved from UniProt
    columnlist = (
        "id,entry name,protein names,length,mass,domains,domain,"
        "families,"
        "go-id,go(molecular function),go(biological process)"
    )

    # This dictionary will be used to populate "blank"/"empty" databases when
    # an error is thrown. Iterables are used as values to avoid problems with
    # "ValueError: If using all scalar values, you must pass an index"
    blank_data = {
        "UniProtKB Entry ID": ["NA"],
        "UniProtKB Entry Name": ["NA"],
        "UniProtKB Protein Names": ["NA"],
        "EC number": ["NA"],
        "Length (Aa)": ["NA"],
        "Mass (Da)": ["NA"],
        "Domains": ["NA"],
        "Domain count": ["NA"],
        "UniProtKB Linked Protein Families": ["NA"],
        "Gene ontology IDs": ["NA"],
        "Gene ontology (molecular function)": ["NA"],
        "Gene ontology (biological process)": ["NA"],
    }

    try:
        # open connection to UniProt(), search and convert result into pandas df
        logger.info(df_row)
        query = f'{df_row[5]} AND organism:"{df_row[0]} {df_row[1]}"'
        search_result = UniProt().search(
            query, columns=columnlist,
        )  # returns empty string for no result
        logger.info(search_result)
        search_result_df = pd.read_table(io.StringIO(search_result))

    except HTTPError:
        logger.warning(
            (
                f"Network error occured when searching UniProt for locus tag:{df_row[5]}.\n"
                "Returning null value 'NA' for all UniProt data"
            )
        )
        return pd.DataFrame(blank_data)

    except EmptyDataError:
        # No UniProt entries found for locus tag, return null data for
        logger.warning(
            (
                f"No data returned from UniProt for locus tag:{df_row[5]}.\n"
                "Returning null value 'NA' for all UniProt data"
            )
        )
        return pd.DataFrame(blank_data)

    # check if multiple entries were returned
    if len(search_result_df) > 1:
        logger.warning(
            (
                f"Multiple hits returned from UniProt for locus tag:{df_row[5]}.\n"
                "Returning null value 'NA' for all UniProt data"
            )
        )
        return pd.DataFrame(blank_data)

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

    # Retrieve EC number from 'UniProtKB Protein Names'
    # or return 'NA' if not included
    EC_search = re.findall(
        r"\(EC [\d-]\d*\.[\d-]\d*\.[\d-]\d*\.[\d-]\d*\)", search_result_df[2]
    )
    if EC_search is None:
        EC_number = "NA"
    else:
        # compiall EC together incase multiple are given
        # and remove EC numbers from protein name
        EC_number = ""
        for EC in EC_search:
            search_result_df[2] = search_result_df[2].replace(EC, "")
            EC = EC.replace("(", "")
            EC = EC.replace(")", "")
            EC_number += EC
            # Do I want to remove the EC number becuase if multiple are given,
            # maybe helpful to see in the protein name section which EC number
            # corresponds to which function
    # Add EC number to dataframe0
    search_result_df.insert(3, "EC number", EC_number)

    return search_result_df


if __name__ == "__main__":
    main()
