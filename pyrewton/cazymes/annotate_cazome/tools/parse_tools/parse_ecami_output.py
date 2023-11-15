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
"""Contains function to parse the output from run-dbCAN. Create standarised output.

:func parse_ecami_output: Coordinate parsing the eCAMI output text file
:func get_ecami_additional_info: Retrieve predicted EC numbers and additional domains
"""

import re

import pandas as pd
import numpy as np


def parse_ecami_output(txt_file_path, logger):
    """Parse the output from the output text file from eCAMI and write out data to a dataframe.

    Retrieves the protein accession/name/identifier, predicated CAZy family, predicated CAZy
    subfamily, predicated EC number and additional/other domains predicated to also be within the
    protein sequence, indicating prediciton of a multiple module enzymes.

    :param text_file_path: path, path to the output text file
    :param logger: logger object

    Return Pandas dataframe containing eCAMI output
    """
    try:
        with open(txt_file_path, "r") as fh:
            ecami_file = fh.read().splitlines()
    except FileNotFoundError:
        logger.warning(
            "Could not find eCAMI output file\n"
            f"{txt_file_path}"
            "Not producing standardised result for this file"
        )
        return None

    # build an empty dataframe to add predication outputs to
    ecami_df = pd.DataFrame(columns=[
        "protein_accession",
        "cazy_family",
        "cazy_subfamily",
        "ec_number",
        "additional_domains",
    ])

    # parse the outputs so in format suitable for final dataframe
    for line in ecami_file:
        if line.startswith(">"):  # identifies new protein
            prediction_output = line.split("\t")

            # Retrieve the CAZy family
            cazy_family = prediction_output[1].split(":")[0]

            # Check CAZy family is formated correctly
            try:
                re.match(r"\D{2,3}\d+", cazy_family.group())
            except AttributeError:
                logger.warning(
                    f"Non-standardised CAZy family name for {line[0][1:]} {cazy_family} in\n"
                    f"{txt_file_path}\n"
                    f"Returning no predicted CAZy family for {line[0][1:]}"
                )
                cazy_family = np.nan

            # Retrieve subfamily, additional domains and EC numbers from the additional info
            additional_domains, ec_number, cazy_subfam = get_ecami_additional_info(
                prediction_output,
                cazy_family,
                txt_file_path,
                logger,
            )

            # build dict to enable easy building of df
            prediction = {
                "protein_accession": [prediction_output[0][1:]],
                "cazy_family": [cazy_family],
                "cazy_subfamily": [cazy_subfam],
                "ec_number": [ec_number],
                "additional_domains": [additional_domains],
            }

            prediction_df = pd.DataFrame(prediction)
            ecami_df = ecami_df.append(prediction_df)

    return ecami_df


def get_ecami_additional_info(prediction_output, cazy_family, txt_file_path, logger):
    """Retrieve additional predicated domains and EC numbers from eCAMI output file.

    :param prediction_output: list of predicted items for a given protein
    :param cazy_family: str, predicted CAZy family

    Return list of additional domains and list of EC numbers as strings, items separated by ', '.
    """
    additional_domains = []
    ec_numbers = []
    subfams = []

    # separate the data stored in the 'additional' data section of the line in the eCAMI output file
    additional_info = prediction_output[2].split("|")

    for item in additional_info:
        item = item.split(":")[0]  # drop the eCAMI group number
        item = item.strip()

        # check if it is an EC number
        if item.find(".") != -1:
            ec_numbers.append(item)
            continue

        # check if it is a subfamily or family prediction
        try:
            re.match(r"(\D{2,3}\d+)|(\D{2,3}\d+_\d+)", item).group()
        except AttributeError:
            logger.warning(
                f"Non-standardised output for {prediction_output[0][1:]} {item} in\n"
                f"{txt_file_path}"
                "Returning no CAZy family/subfamily"
            )
            continue

        # Check if subfamily
        if item.find("_") != -1:  # predicated CAZy subfamily
            # check format
            try:
                re.match(r"\D{2,3}\d+_\d+", item).group()
                # check if it is the subfamily of the predicated CAZy family
                if cazy_family == item[:item.find("_")]:
                    subfams.append(item)
                else:
                    additional_domains.append(item)
            except AttributeError:
                logger.warning(
                    "Non-standardised CAZy subfamily name for "
                    f"{prediction_output[0][1:]} {item} in\n"
                    f"{txt_file_path}"
                    "Returning subfamily"
                )

        else:  # predicated CAZy family
            # check format
            try:
                re.match(r"\D{2,3}\d+", item).group()
                additional_domains.append(item)
            except AttributeError:
                logger.warning(
                    f"Non-standardised CAZy family name for {prediction_output[0][1:]} {item} in\n"
                    f"{txt_file_path}"
                    "Returning no CAZy family"
                )

    # convert empty lists to null values if necessary, or make list more human readable
    if len(additional_domains) == 0:
        additional_domains = np.nan
    else:
        additional_domains = ", ".join(additional_domains)

    if len(ec_numbers) == 0:
        ec_numbers = np.nan
    else:
        ec_numbers = ", ".join(ec_numbers)

    if len(subfams) == 0:
        subfams = np.nan
    else:
        subfams = ", ".join(subfams)

    return(additional_domains, ec_numbers, subfams)