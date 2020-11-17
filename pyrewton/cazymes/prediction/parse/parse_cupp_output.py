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
        ec_number = get_cupp_ec_number(prediction_output, log_file_path, logger)

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


def get_cupp_ec_number(prediction_output, log_file_path, logger):
    """Retrieve the predicted EC numbers from the CUPP output log file.

    EC numbers are represented as "CAZy_fam:EC_number" in the log file.
    If multiple EC numbers are predicated and the CAZy families of each are
    given then these are separated by '-'. If multiple EC numbers are predicated
    for the same CAZy family, these are separated by '&'.

    :param prediciton_output: list of items from log file line.
    :param log_file_path: path to CUPP output file
    :param logger: logger object

    Return string if EC numbers are given, or null value if not.
    """
    # Separate CAZy families from their respective EC numbers
    ec_data = prediction_output[-2].split(":")

    # create empty list to store all predicted EC numbers in
    ec_numbers = []

    # If multiple CAZy families are listed these are separed by '-'
    # and will appear on the end of an item in ec_data

    for item in ec_data:
        item = item.strip()

        # check that it isn't a CAZy family:
        try:
            re.match(r"\d", item).group()
        except AttributeError:  # raised if first character is not a digit, thus not an EC number
            continue

        # check if CAZy family is appendaged to EC number
        if item.find("-") != -1:
            item = item[:(item.find("-"))]  # remove appendaged CAZy family

        # separate out EC numbers if multiple were predicted
        item = item.split("&")

        # parse each predicted EC number
        for ec in item:
            # check formating, e.g. somtimes 'Unknown' is written in stead of an EC number
            if (ec == "Unknown") or (ec == "unknown"):
                print("unknown=", ec)
                continue
            else:
                try:
                    re.match(r"\d+?\.(\d+?|\*)\.(\d+?|\*)\.(\d+?|\*)", ec).group()
                    #  standardise missing digits in ec number to '-'
                    ec = ec.replace("*", "-")
                    ec_numbers.append(ec)
                except AttributeError:
                    logger.warning(
                        f"Non-standard EC# for {prediction_output[0]}, {ec} in\n"
                        f"{log_file_path}"
                    )

    # Convert lists into more human readable strings
    if len(ec_numbers) == 0:
        ec_number = np.nan
    else:
        ec_number = ", ".join(ec_numbers)

    return(ec_number)
