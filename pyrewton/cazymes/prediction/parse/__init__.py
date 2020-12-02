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
"""Module to parse the output from dbCAN, CUPP and eCAMI. Create standarised output."""


import json
import re

import pandas as pd


def get_cazy_accessions(args, logger):
    """Retrieve dataframe of CAZymes catalogued in CAZy and dictionary of GenBank synonyms.

    Returns the CAZy family, GenBank and UniProt columns from the CAZy dataframe. GenBank synonyms
    are the GenBank accessions listed under a CAZyme in addition to the first GenBank accession
    listed for the CAZyme in CAZy and which is always hyperlinked to GenBank. The GenBank synonyms
    are frequently annotated as 'identical proteins' to the CAZyme first GenBank accession in
    GenBank.

    :param args: cmd args parser
    :param logger: logger object

    Return Pandas dataframe.
    """
    # open the cazy df
    try:
        cazy_df = pd.read_csv(args.stat, header=0, index_col=0)
    except FileNotFoundError:
        logger.error(
            "Could not find file containing CAZy dataframe.\n"
            "Check the given path is correct\n."
            "Statistical evaluation will not be performed."
            "Only reports summarising the predictions will be produced.\n"
        )
        return None, None

    # Parse the CAZy dataframe to retrieve only the CAZy family, GenBank and UniProt columns
    cazy_df = cazy_df.drop(["Protein_name", "EC#", "Source_organism", "PDB/3D"], axis=1)

    genbank_synonyms = get_genbank_synonyms(args, logger)

    return cazy_df, genbank_synonyms


def get_genbank_synonyms(args, logger):
    """Retrieve dictionary of GenBank synonyms from JSON file.

    :param args: cmd args parser
    :param logger: logger object

    Return dictionary of GenBank synonyms.
    """
    # Compile path to GenBank synonym dictionary
    cazy_path = str(args.cazy)
    try:
        cazy_filename = re.match(r"cazy_.+?\.csv", cazy_path, re.IGNORECASE).group()
        dict_filename = cazy_filename.replace("cazy_", "cazy_genbank_synonyms_")
        dict_filename = dict_filename.replace(".csv", ".json")
        dict_path = cazy_path.replace(cazy_filename, dict_filename)
    except AttributeError:
        logger.warning(
            "Did not find a JSON file of GenBank synonyms from CAZy.\n"
            "Continuing with evaluating CAZyme prediction tools without GenBank synonyms."
        )
        dict_path = None
        genbank_synonyms = {}

    # Open JSON file
    if dict_path is not None:
        try:
            with open(dict_path, "r") as fh:
                genbank_synonyms = json.load(fh)
        except FileNotFoundError:
            logger.error(
                "Could not open the JSON file of GenBank synonyms from CAZy.\n"
                "Continuing with evaluating CAZyme prediction tools without GenBank synonyms."
            )

    return genbank_synonyms
