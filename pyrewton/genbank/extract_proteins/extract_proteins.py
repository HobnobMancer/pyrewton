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
"""Retrieve all protein annotations from GenBank files.

:cmd_args df_input: path, path to input dataframe - required
:cmd_args genbank: path, path to directory containing GenBank files - required

:cmd_args --force: bool, force overwriting files in output directory
:cmd_args --log: path, path to direct writing out log file
:cmd_args --nodelete: not delete existing files in output directory
:cmd_args --output: path, path to output directory
:cmd_args --output_df: path, path to output dataframe, including file name
:cmd_args --verbose: boolean, enable verbose logging

:func main: Coordinate calling of other functions
:func create_dataframe: build dataframe summarising CAZy annotation in GenBank files
:func create_df_foundation: Parse input dataframe row
:func build_df_foundation: Compile row data for dataframe
:func get_genbank_protein_data: Retrieve protein name and IDs

Generate summary dataframe and of all annotated proteins in all GenBank
files directly linked to a given species.
"""

import gzip
import logging
import re
import sys

from pathlib import Path
from typing import List, Optional

import pandas as pd

from Bio import SeqIO
from saintBioutils.utilities.file_io import get_paths
from tqdm import tqdm

from pyrewton.utilities import config_logger
from pyrewton.utilities.parsers.cmd_parser_compile_annotations import build_parser
from pyrewton.utilities.file_io import make_output_directory, write_out_dataframe


def main(argv: Optional[List[str]] = None, logger: Optional[logging.Logger] = None):
    """Coordinate the retrieval of protein annotations from GenBank (.gbff) files.

    Including building parser, logger and output directory.

    Return dataframe of protein data.
    """
    if argv is None:
        parser = build_parser()
        args = parser.parse_args()
    else:
        parser = build_parser(argv)
        args = parser.parse_args()

    if logger is None:
        config_logger(args)
    logger = logging.getLogger(__name__)

    # If specified output directory, create output directory to write FASTA files too
    if args.output is not sys.stdout:
        make_output_directory(args.output, args.force, args.nodelete)

    # Open input dataframe
    logger.info("Opening input dataframe %s", args.input_df)
    input_df = pd.read_csv(args.input_df, header=0, index_col=0)

    # build dict of genomic assembly paths {acc: path}
    genome_dict = get_genomic_assembly_paths(args)

    # Build dataframe containing extracted proteome data
    protein_annotation_df = create_dataframe(input_df, genome_dict, args)

    # Write out dataframe
    if args.output_df is not None:
        write_out_dataframe(protein_annotation_df, args.output_df, args.force)

    # Write out FASTA files
    index = 0
    for index in tqdm(
        range(len(protein_annotation_df["Genus"])), desc="Writing protein to FASTA"
    ):
        df_row = protein_annotation_df.iloc[index]
        write_fasta(df_row, logger, args)
        index += 1

    logger.info("Programme finsihed. Terminating.")


def get_genomic_assembly_paths(args):
    """Compile dict of paths to genomic assemblies
    
    :param args: cmd-line args parser
    
    Return dict {accession: Path}
    """
    logger = logging.getLogger(__name__)

    # retrieve path to all genomic assemblies in the input dir
    genome_paths = get_paths.get_file_paths(args.genbank, suffixes=[".gbff.gz"])

    if len(genome_paths) == 0:
        logger.error(
            f"No genomic assemblies retrieved from {args.genbank}\n"
            "Make sure the input path is correct\n"
            "Terminating program"
        )
        sys.exit(1)

    genome_dict = {}  # {accession: Path}
    
    for _path in genome_paths:
        # extract the accession
        accession = f"{_path.name.split('_')[0]}_{_path.name.split('_')[1]}.{_path.name.split('_')[2]}"

        try:
            genome_dict[accession]
            logger.warning(
                f"Retrieved multiple genomic assemblies with the accession {accession}\n"
                f"Will extract proteins from the file {_path}\n"
                f"NOT {_path}\n"
                "As this was the first genome retrieved with this accession"
            )
        
        except KeyError:
            genome_dict[accession] = _path

    return genome_dict


def create_dataframe(input_df, genome_dict, args):
    """Build datafame containing all protein annotations in GenBank files.

    Iterate over input dataframe row wise. This allows for converting human
    readable accession number list from a string to a Python list.

    Per species, retrieve all protein names and IDs for every accession number,
    therefore, will return multiple rows with the same accession number, but unique protein ID.

    Append data for all rows for a single species to the tuple 'all_foundation_data', so that all
    data it compiled together and can be used simultaneously to populate a pandas dataframe without
    risk of data missalignment or over writing.

    :param input_df: pandas dataframe
    :param genome_dict: dict keyed by GenBank accession and valued by Path() instance to a genomic
        assembly
    :param args: parser arguments

    Return dataframe.
    """
    # Create empty dataframe to add data to
    protein_annotation_df = pd.DataFrame(
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
        protein_annotation_df = protein_annotation_df.append(
            compile_annotations(input_df.iloc[df_index], genome_dict, args),
            ignore_index=True,
        )
        df_index += 1

    return protein_annotation_df


def compile_annotations(df_row, genome_dict, args):
    """Prepare row data to create dataframe.

    Coordinate retrieval of protein data from GenBank files for every accession number for the
    species passed to the function. Store data in a dataframe: Genus, Species, Tax ID, Accession
    number, protein ID, locus tag, gene location, product.

    Each row in dataframe contains a unique protein, and thus multiple rows will have the same
    accession number. Any failed to retrieve data will be returned as pandas null value 'NA'.

    Reminder of panda series (referred to as df_row) structure:
    df_row[0] = "Genus"
    df_row[1] = "Species"
    df_row[2] = "Taxonomy ID"
    df_row[3] = "Accession list" - this is human readable list, stored as a str.

    :param df_row: row from input_df (dataframe)
    :param genome_dict: dict keyed by GenBank accession and valued by Path() instance to a genomic
        assembly
    :param args: parser arguments

    Return dataframe.
    """
    logger = logging.getLogger(__name__)

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
        protein_data = get_annotations(accession, genome_dict, args)  # list

        # check if any data retrieved
        if len(protein_data) == 0:
            logger.warning(
                f"No protein data retrieved for {accession} from GenBank file.\n"
                "Most likely cause is GenBank file contained no CDS type features."
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


def get_annotations(accession_number, genome_dict, args):
    """Retrieve protein ID, locus tag and function from GenBank file.

    From each record the protein ID, locus tag, location and annotated function is retrieved, and
    stored as a list. Lists wil be added to a single tuple containing all protein data. Any failed
    to retrieve data will be returned as pandas null value 'NA'.

    :param accession_number: str
    :param genome_dict: dict keyed by GenBank accession and valued by Path() instance to a genomic
        assembly
    :param genbank_input: path, path to directory containing GenBank files

    Return tuple.
    """
    logger = logging.getLogger(__name__)

    # check if accession number was provided
    if accession_number == "NA":
        logger.warning(
            f"Null value ('NA') was contained in cell for {accession_number},"
            "exiting retrieval of protein data.\nReturning null ('NA') value"
            "for all protein data"
        )
        return ["NA", "NA", "NA", "NA", "NA"]

    try:
        gb_file = genome_dict[accession_number]
    except KeyError:
        # no genomic assembly retrieved for this accession
        return ["NA", "NA", "NA", "NA", "NA"]

    # create empty list to store protein data
    all_protein_data = []

    # Retrieve protein data from GenBank file

    with gzip.open(gb_file, "rt") as handle:
        for gb_record in SeqIO.parse(handle, "genbank"):
            for (index, feature) in enumerate(gb_record.features):
                # empty protein data list so as not to contaminate data of next protein
                protein_data = []
                
                # Parse over only protein encoding features (type = 'CDS')
                if feature.type == "CDS":
                    protein_data.append(
                        get_record_feature(feature, "protein_id", accession_number)
                    )
                    protein_data.append(
                        get_record_feature(feature, "locus_tag", accession_number)
                    )
                    protein_data.append(
                        get_record_feature(feature, "location", accession_number)
                    )
                    # extract annotated function of product
                    protein_data.append(
                        get_record_feature(feature, "product", accession_number)
                    )
                    # extract protein sequence
                    protein_data.append(
                        get_record_feature(feature, "translation", accession_number)
                    )

                    # add protein data to total protein data list, only if data was retrieved
                    if len(protein_data) == 5:
                        # if null value was returned for every feature attribute log error
                        # and don't add to all_protein_data list
                        if protein_data == ["NA", "NA", "NA", "NA", "NA"]:
                            logger.warning(
                                    f"No data retrieved from CDS type feature, index: {index}, "
                                    f"accession: {accession_number}"
                            )
                        # if some data retrieved, add to all_protein_list
                        else:
                            all_protein_data.append(protein_data)

                    else:
                        # error occured in that one of the appending actions failed to append
                        # and would lead to misalignment in the dataframe if added to the
                        # all_protein_data list
                        logger.warning(
                            f"Error occured during retrieval of data from feature, {index}\n"
                            f"for {accession_number}. Returning no protein data"
                        )

    return all_protein_data


def get_record_feature(feature, qualifier, accession):
    """Retrieve data from GenBank record feature.

    :param feature: feature object, GenBank file record feature
    :param qualifier: str, key of feature attribute

    Return data from GenBank record feature, or "NA" if failed to retrieve.
    """
    logger = logging.getLogger(__name__)
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
                f"Failed to retrieve feature location, returning 'NA', accession: {accession}",
            )
            return "NA"
    else:
        try:
            data = feature.qualifiers[qualifier][0]
            return data
        except KeyError:
            logger.warning(
                f"Failed to retrieve feature {qualifier}, returning 'NA', accession: {accession}",
            )
            return "NA"


def write_fasta(df_row, args, filestem="genbank_proteins"):
    """Write out FASTA file.

    :param df_row: row from pandas df of UniProt search results
    :param filestem: str, FASTA file name
    :param args: parser arguments

    Returns nothing.
    """
    # Create file content

    # FASTA sequences have 60 characters per line, add line breakers into protein sequence
    # to match FASTA format
    sequence = df_row["Protein Sequence"]
    sequence = "\n".join([sequence[i : i + 60] for i in range(0, len(sequence), 60)])
    # Retrieve protein ID
    protein_id = df_row["NCBI Protein ID"] + " " + df_row["Locus Tag"]

    file_content = f">{protein_id} \n{sequence}\n"

    # Create file name

    # Retrieve Taxonomy ID and add txid prefix is present
    tax_id = str(df_row["NCBI Taxonomy ID"])
    if tax_id.startswith("NCBI:txid") is False:
        tax_id = "txid" + tax_id
    else:
        tax_id = tax_id.replace("NCBI:txid", "txid")
    tax_id.replace(" ", "_")
    # Retrieve accession number
    accession = df_row["NCBI Accession Number"]
    # remove characters that could make file names invalid
    invalid_file_name_characters = re.compile(r'[,;./ "*#<>?|\\:]')
    tax_id = re.sub(invalid_file_name_characters, "_", tax_id)
    accession = re.sub(invalid_file_name_characters, "_", accession)

    # Create output path
    if args.output is not sys.stdout:
        output_path = args.output / f"{filestem}_{tax_id}_{accession}.fasta"
        # Write out data to Fasta file
        with open(output_path, "a") as fh:
            fh.write(file_content)
    else:
        binary_file_content = bytearray(file_content, "utf8")
        sys.stdout.buffer.write(binary_file_content)

    return


if __name__ == "__main__":
    main()
