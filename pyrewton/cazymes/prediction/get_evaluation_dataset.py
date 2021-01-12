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
"""Retrieves datasets for evaluating the CAZyme prediction tools: dbCAN, CUPP and eCAMI.

:cmd_args   :

:func    :

Writes out a FASTA per candidate species, containing all the protein sequences to analysed by the
CAZymes prediction tools. The datasets contain an equal number of CAZymes to non-CAZymes.
"""

import gzip
import logging
import os
import random
import re
import sys

from pathlib import Path
from typing import List, Optional
from tqdm import tqdm

import numpy as np
import pandas as pd

from Bio import SeqIO
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from sqlalchemy.ext.declarative import declarative_base

from pyrewton.utilities import build_logger
from pyrewton.utilities.cmd_get_evaluation_dataset import build_parser
from pyrewton.utilities.file_io import make_output_directory


# Use the declarative system
Base = declarative_base()
Session = sessionmaker()


def main(argv: Optional[List[str]] = None, logger: Optional[logging.Logger] = None):
    """Coordinate preparation for script, and terminating the programme when finished."""
    # programme preparation

    # build cmd-line arguments parser
    # Check if namepsace isn't passed, if not parse command-line
    if argv is None:
        # Parse command-line
        parser = build_parser()
        args = parser.parse_args()
    else:
        parser = build_parser(argv)
        args = parser.parse_args()

    # build logger
    # Note: log file is only created if specified at cmd-line
    if logger is None:
        logger = build_logger("get_ncbi_genomes", args)
    logger.info("Run initated")

    # open session to local CAZy database (this can be curated by using the
    # cazy_webscraper, available at https://github.com/HobnobMancer/cazy_webscraper)
    session = get_cazy_db_session(args, logger)

    # retrieve paths to FASTA files
    fasta_files = get_fasta_file_paths(args, logger)

    # create dataset per FASTA file
    for fasta in fasta_files:
        build_protein_dataframe(fasta, session, args, logger)

    logger.info(
        "Finished creating a dataset for each input FASTA file.\n"
        "Terminating Progamme."
    )


def get_cazy_db_session(args, logger):
    """Retrieve an open database session to a local CAZy database.

    :param args: cmd-line args parser
    :param logger: logger object

    Return open database session.
    """
    logger.info("Opening session to the local copy of the CAZy database.")

    # check if database file exists
    if os.path.exists(args.database) is False:
        logger.error(
            "Path provided for local CAZy database does not exist.\n"
            "Cannot proceed without a valid local CAZy database.\n"
            "Terminating program."
        )
        sys.exit(1)

    # build database engine
    try:
        engine = create_engine(f"sqlite+pysqlite:///{args.database}", echo=False)
        Base.metadata.create_all(engine)
        Session.configure(bind=engine)
    except Exception:
        logger.error(
            (
                "Was unable to open local CAZy database. The causing error is presented below\n"
                "Cannot proceed without a valid local CAZy database.\n"
                "Terminating program."
            ),
            exc_info=1,
        )

    return Session()


def get_fasta_file_paths(args, logger):
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


def build_protein_dataframe(fasta_path, session, args, logger):
    """Build a dataframe containing the protein data for the current working input FASTA file."""
    logger.info(f"Create dataset for {fasta_path}")

    # create dictionary to store data that will go in the database
    protein_dict = {"protein_data": [], "sequence": []}

    # open the current FASTA file and populate dictionary with its data
    for seq_record in SeqIO.parse(fasta_path, "fasta"):
        protein_dict["protein_data"].append(seq_record.id)
        protein_dict["sequence"].append(seq_record.seq)

    # build dataframe of protein sequences
    protein_df = pd.DataFrame(protein_dict)

    # add new column populate with CAZyme ('1') and non-CAZyme ('0') classification
    protein_df = get_cazy_classification(protein_df, session, logger)

    # write out a random dataset
    get_dataset(protein_df, args, logger)

    return


def get_cazy_classification(protein_df, session, logger):
    """For each protein check if classified as a CAZyme by CAZy, and its annotated CAZy families.

    Adds two new columnd to the protein dataframe:
    cazyme_classification: 1=CAZymes, 0=non-CAZymes
    cazy_families: list of CAZy families of CAZymes, each family is a separate str in the list.

    Return protein dataframe containing CAZyme classification.
    """
    logger.info("Adding CAZy data to protein dataframe")

    index = 0
    for index in tqdm(
        range(len(protein_df["protein_data"])),
        desc="Adding CAZy data to protein df"
    ):
        df_row = protein_df.iloc[index]
        genbank_accession = df_row[0].split(" ")[0]

        # query local CAZy database to see if GenBank accession is contained
        query = session.query(Genbank).filter_by(genbank_accession=genbank_accession).all()

        # if protein is a CAZyme
        if len(query) != 0:
            cazyme_classification = 1
            # Retrieve CAZy families

        # if protein is not a CAZyme
        else:
            cazyme_classification = 0
            cazy_families = np.nan

        # add CAZy classification and families to protein df
        try:
            protein_df.insert(3, "cazyme_classification", cazyme_classification)
            protein_df.insert(4, "cazy_families", cazy_families)
        except ValueError:
            logger.warning(
                "Failed to insert CAZyme classification into the protein dataframe,\n"
                "column already present"
            )

    return protein_df


def get_dataset(protein_df, fasta_file_path, args, logger):
    """Retrieve dataset of equal number of CAZymes and non-CAZymes."""
    # create dataframe only of non-CAZymes
    non_cazyme_rows = protein_df["cazyme_classification"] == 0
    non_cazyme_df = protein_df[non_cazyme_rows]

    # create dataframe of only CAZymes
    cazyme_rows = protein_df["cazyme_classification"] == 1
    cazyme_df = protein_df[cazyme_rows]

    if args.sample_size == None:
        sample_size = (len(cazyme_df["protein_data"]) / 2)
    else:
        sample_size = (args.sample_size / 2)

    # select random rows from non-CAZymes df
    non_cazyme_df = non_cazyme_df.sample(n=sample_size)

    # select random rows from CAZymes df
    cazyme_df = cazyme_df.sample(n=sample_size)

    # combine the non-CAZyme and CAZyme dfs
    dfs = [non_cazyme_df, cazyme_df]
    protein_df = pd.concat(dfs)

    # add the proteins in the new protein dataframe to a FASTA file in a random order
    protein_order = random.sample(
        range(0, len(protein_df["protein_data"])),
        (len(protein_df["protein_data"]))
    )

    for index in protein_order:
        add_protein_to_output_fasta(protein_df.iloc[index], fasta_file_path, args, logger)

    return


def add_protein_to_output_fasta(df_row, fasta_file_path, args, logger):
    """Write out the dataset of CAZymes and non-CAZymes to a single FASTA file."""
    return
