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
"""Contains function to parse the output from run-dbCAN. Create standarised output.

:func parse_cupp_output: Coordinate parsing the CUPP output log file
:func get_cupp_domain_range: Retrieve the domain ranges of predicted domains
:func get_cupp_ec_number: Retrieve the predicted EC numbers of predicted domains
"""

import re

import pandas as pd
import numpy as np


def parse_cupp_output(log_file_path, logger):
    """Parse the output from the output log file from CUPP and write out data to a dataframe.

    Retrieves the protein accession/name/identifier, predicated CAZy family, predicated CAZy
    subfamily, predicated EC number and predicated range of domain within the protein sequence
    (the index of the first and last residues of the domain).

    :param log_file_path: path, path to the output log file
    :param logger: logger object

    Return Pandas dataframe containing CUPP output
    """
    try:
        with open(log_file_path, "r") as lfh:
            log_file = lfh.read().splitlines()
    except FileNotFoundError:
        logger.warning(
            "Could not find CUPP output file\n"
            f"{log_file_path}"
            "Not producing standardised result for this file"
        )
        return None

    # build an empty dataframe to add predication outputs to
    cupp_df = pd.DataFrame(columns=[
        "protein_accession",
        "cazy_family",
        "cazy_subfamily",
        "ec_number",
        "domain_range",
    ])

    # add in predictions from log file to dataframe
    for line in log_file:
        prediction_output = line.split("\t")

        # retrieve domain range if given
        domain_range = get_cupp_domain_range(prediction_output, log_file_path, logger)

        # retrieve EC number if given
        ec_number = get_cupp_ec_number(prediction_output)

        # retrieve predicated CAZy subfamily
        subfam = prediction_output[-1]
        # change empty string to null value
        if len(subfam) == 0:
            subfam = np.nan
        else:
            # When writing the log file sometimes the ':' raw data is not converted to '_'
            # ensure this change is made to standardise CAZy subfamily naming
            subfam = subfam.replace(":", "_")
            subfam = subfam.replace("+", ", ")  # separate subfams ', ' for human readability

        # build dict to enable easy building of df
        prediction = {
            "protein_accession": [prediction_output[0]],
            "cazy_family": [prediction_output[1]],
            "cazy_subfamily": [subfam],
            "ec_number": [ec_number],
            "domain_range": [domain_range],
        }

        prediction_df = pd.DataFrame(prediction)
        cupp_df = cupp_df.append(prediction_df)

    return cupp_df


def get_cupp_domain_range(prediction_output, log_file_path, logger):
    """Retrieve the amino acid domain_range from CUPP output log file.

    :param prediciton_output: list of items from log file line.
    :param log_file_path: path to CUPP output log file
    :param logger: logger object

    Return string if domain range given, or null value if not.
    """
    for item in prediction_output:
        if item.find("..") != -1:
            domain_range = item
            break
        else:
            domain_range = np.nan

    if type(domain_range) != float:  # if domain_range != np.nan
        try:
            re.match(r"\d+?\.\.\d+", domain_range).group()
        except AttributeError:
            logger.warning(
                "Incorrect parsing to retrieve domain ranges for\n"
                f"{prediction_output[0]} in\n"
                f"{log_file_path}"
            )
            domain_range = np.nan

    return domain_range


def get_cupp_ec_number(prediction_output):
    """Retrieve the predicted EC numbers from the CUPP output log file.

    EC numbers are represented as "CAZy_fam:EC_number" in the log file.
    If multiple EC numbers are predicated and the CAZy families of each are
    given then these are separated by '-'. If multiple EC numbers are predicated
    for the same CAZy family, these are separated by '&'.

    :param prediciton_output: list of items from log file line.

    Return string if EC numbers are given, or null value if not.
    """
    ec_data = prediction_output[-2].split(":")

    if len(ec_data) == 1:
        ec_number = np.nan

    elif len(ec_data) == 2:
        ec_numbers = []
        ec_split = ec_data[1].split("&")  # multiple predicated EC numbers are separated by "&"

        for item in ec_split:
            # check formating, e.g. somtimes 'Unknown' is written in stead of an EC number
            try:
                re.match(r"\d+?\.(\d+?|\*)\.(\d+?|\*)\.(\d+?|\*)", item)
            except AttributeError:
                ec_numbers.append(np.nan)
                continue
            ec_numbers.append(item)

        ec_number = ", ".join(ec_numbers)

    else:
        ec_numbers = []
        for item in ec_data:
            # check first character is a digit and thus item contains an EC number
            try:
                re.match(r"\d.+", item).group()
            except AttributeError:
                continue

            # check if item contains aditional info
            if item.find("-") != -1:
                cut_off = item.find("-")
                ec = item[:cut_off]
                ec_split = ec.split("&")

                for i in ec_split:
                    try:
                        re.match(r"\d+?\.(\d+?|\*)\.(\d+?|\*)\.(\d+?|\*)", i)
                    except AttributeError:
                        continue
                    ec_numbers.append(i)

            else:
                ec = item.split("&")
                for i in ec:
                    try:
                        re.match(r"\d+?\.(\d+?|\*)\.(\d+?|\*)\.(\d+?|\*)", i)
                    except AttributeError:
                        continue
                    ec_numbers.append(i)

        ec_number = ", ".join(ec_numbers)

    # standardise missing digits in ec number to -
    if type(ec_number) != float:  # if ec_number != np.nan
        ec_number = ec_number.replace("*", "-")

    return(ec_number)
