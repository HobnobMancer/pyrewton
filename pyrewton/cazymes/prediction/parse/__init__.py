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


import pandas as pd


def get_cazy_accessions(args, logger):
    """Retrieve GenBank and UniProt accession numbers of proteins in CAZy dataframe.

    :param args: cmd args parser
    :param logger: logger object

    Return two lists of accession numbers.
    """
    # open the cazy df
    cazy_df = 

    # create empty list to store accessions in
    genbank_accessions = []
    uniprot_accessions = []

    # get genbank accession numbers
    index = 0
    for index in range(len(cazy_df["Protein_name"])):
        cell_content = cazy_df[index, 4]
        cell_content = cell_content.split(",\n")
        i = 0
        for i in range(len(cell_content)):
                accession = cell_content[i][:(cell_content[1].find(" >
                genbank_accessions.append(accession)

    # get uniprot accession numbers
    index = 0
    for index in range(len(cazy_df["Protein_name"])):
        cell_content = cazy_df[index, 5]
        cell_content = cell_content.split(",\n")
        i = 0
        for i in range(len(cell_content)):
                accession = cell_content[i][:(cell_content[1].find(" >
                uniprot_accessions.append(accession)

    return genbank_accessions, uniprot_accessions
