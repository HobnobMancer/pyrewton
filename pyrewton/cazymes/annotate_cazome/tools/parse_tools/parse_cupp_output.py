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

    Return dict
    {toolname: {protein acc: {family: [domain range]}}}
    or None if fails
    """
    tool_predictions = {'CUPP': {}}

    try:
        cupp_df = pd.read_table(log_file_path, header=None)
    except FileNotFoundError:
        logger.warning(
            "Could not find CUPP output file\n"
            f"{log_file_path}"
            "Not add annotations from CUPP to db"
        )
        return None

    # add in predictions from log file to dataframe
    for ri in range(len(cupp_df)):
        prediction_output = cupp_df.iloc[ri]
        prot_acc = prediction_output[0].split(" ")[0]
        cazy_family = prediction_output[1]

        # retrieve domain range if given
        domain_range = get_cupp_domain_range(prediction_output[6], log_file_path, logger)

        try:
            tool_predictions['CUPP'][prot_acc]

            try:
                tool_predictions['CUPP'][prot_acc][cazy_family].append(domain_range)
            except KeyError:
                tool_predictions['CUPP'][prot_acc][cazy_family] = [domain_range]

        except KeyError:
            tool_predictions['CUPP'][prot_acc] = {cazy_family: None}
        
    return tool_predictions


def get_cupp_domain_range(domain_range, log_file_path, logger):
    """Retrieve the amino acid domain_range from CUPP output log file.

    :param domain_range: str from domain range column in cupp df
    :param log_file_path: path to CUPP output log file
    :param logger: logger object

    Return string if domain range given, or null value if not.
    """
    if type(domain_range) != float:  # if domain_range != np.nan
        try:
            re.match(r"\d+?\.\.\d+", domain_range).group()
        except AttributeError:
            logger.warning(
                "Incorrect parsing to retrieve domain ranges for\n"
                f"{domain_range} in\n"
                f"{log_file_path}"
                "Retunring no domain range"
            )
            domain_range = None

    return domain_range
