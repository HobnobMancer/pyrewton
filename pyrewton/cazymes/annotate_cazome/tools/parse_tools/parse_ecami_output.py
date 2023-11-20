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

    Return dict
    {toolname: {protein acc: {family: [domain range]}}}
    or None if fails
    """
    tool_predictions = {'eCAMI': {}}

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

    # parse the outputs so in format suitable for final dataframe
    for line in ecami_file:
        if line.startswith(">"):  # identifies new protein
            prediction_output = line.split("\t")
            prot_acc = prediction_output[0]

            cazy_families = set()
            for predictions in prediction_output[1:]:
                domains = predictions.split("|")
                for domain in domains:
                    if domain.split(":")[0].startswith(("G", "P", "C", "A")):  # do not add EC numbers
                        cazy_families.add(domain.split(":")[0])

            for cazy_family in cazy_families:
                try:
                    tool_predictions['eCAMI'][prot_acc][cazy_family] = [None]
                except KeyError:
                    tool_predictions['eCAMI'][prot_acc] = {cazy_family: None}

    return tool_predictions
