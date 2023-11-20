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

:func parse_dbcan_output: Coordinate parsing dbCAN overview.txt, build df for each prediction tool
:func parse_hmmer_output: Parse the HMMER output in overview.txt
:func parse_hotpep_output: Parse the Hotpep output in overview.txt
:func parse_diamond_output: Parse the DIAMOND output in overview.txt
:func get_dbcan_consensus: Get the consensus predictions across dbCAN
:func add_hotpep_ec_predictions: Add the EC# predictions from Hotpep.out to the Hotpep dataframe
"""

import re

import pandas as pd
import numpy as np

from tqdm import tqdm


def parse_dbcan_output(overview_file_path, logger, dbcan_version):
    """Parse the output from the dbCAN overview.txt file, write a dataframe for each prediction tool

    :param overview_file_path: path, path to the output 'overview.txt' file
    :param logger: logger object
    :param dbcan_version: int, firt number of version number

    Return dict
    {toolname: {protein acc: {family: domain range}}}
    or None if fails
    """
    tool_predictions = {}
    if dbcan_version == 2:
        kmer_tool = 'Hotpep'
    elif dbcan_version == 3:
        kmer_tool = 'eCAMI'
    else:
        kmer_tool = 'dbCAN-sub'

    try:
        df = pd.read_table(overview_file_path/"overview.txt")
    except FileNotFoundError:
        logger.error(f"Could not find overview.txt file in {overview_file_path.name}\nSkipping output dir")
        return

    for ri in tqdm(range(len(df)), desc=f"Parsing {overview_file_path.name}"): 
        row = df.iloc[ri]
        protein_acc = row[0]

        # get one dict per tool {fam: [domain ranges]}
        if list(df.columns)[1].startswith('EC#'):  # versions >= 3 include additional EC# col
            hmmer_fams = parse_hmmer_output(row[2], logger)
            hotpep_fams = parse_hotpep_output(row[3], logger)  # works for hotpep and ecami
            diamond_fams = parse_diamond_output(row[4], logger)
        else:
            hmmer_fams = parse_hmmer_output(row[1], logger)
            hotpep_fams = parse_hotpep_output(row[2], logger)
            diamond_fams = parse_diamond_output(row[3], logger)


        # retrieve consensus results for dbCAN
        if row[-1] != "1":   # check ''#ofTools', don't include if no consensus
            # consensus is any prediction at at least two tools agree upon
            consensus_cazy_fams = get_dbcan_consensus(
                list(hmmer_fams.keys()),
                list(hotpep_fams.keys()),
                list(diamond_fams.keys()),
            )  # set of consensus CAZy family annotations
        else:
            consensus_cazy_fams = []

        try:
            tool_predictions[protein_acc]
        except KeyError:
            tool_predictions[protein_acc] = {
                'dbCAN': {},
                'HMMER': {},
                'DIAMOND': {},
                kmer_tool: {},
            }

        for family in consensus_cazy_fams:
            try:
                tool_predictions[protein_acc]['dbCAN'][family].append(None)
            except KeyError:
                tool_predictions[protein_acc]['dbCAN'][family] = [None]

        for family in hmmer_fams:
            for domain_range in hmmer_fams[family]:
                try:
                    tool_predictions[protein_acc]['HMMER'][family].append(domain_range)
                except KeyError:
                    tool_predictions[protein_acc]['HMMER'][family] = [domain_range]

        for family in diamond_fams:
            try:
                tool_predictions[protein_acc]['DIAMOND'][family].append(None)
            except KeyError:
                tool_predictions[protein_acc]['DIAMOND'][family] = [None]

        for family in hotpep_fams:
            try:
                tool_predictions[protein_acc][kmer_tool][family].append(None)
            except KeyError:
                tool_predictions[protein_acc][kmer_tool][family] = [None]


    return tool_predictions


def parse_hmmer_output(tool_data, logger):
    """Parse the output from HMMER from the overview.txt file.

    Retrieve the protein accession, predicated CAZy families, accomanying subfamiles and domain
    ranges. Multiple domains can be predicated, and the families, subfamiles and ranges will be
    collected as as lists, with the index for each item in each list being associated with the
    matching items (by index) in the other lists. For example, the first item in every list will
    contain the prediciton for the first predicated domain, the second item the second prediction
    etc.

    Two dataframes are produced in parallel. The first is to build a human readble dataframe for
    writing out to a .csv file. The second stores the data in lists to build the dbCAN consensus
    dataframe later on.

    :param tool_data: str, contents from HMMER cell in dbCAN df
    :param logger: logger object

    Return dict {fam: [domain ranges]}
    """
    # Retrieve predictions for HMMER and separate out each the predicated domains
    hmmer_predictions = tool_data.split("+")

    fam_annotations = {}  # {fam: [domain ranges]}

    # check if no CAZy prediction was made for protein represented in 'line'
    # if no family/subfamily (and thus domain ranges) predicted result will be '-' in 'line'
    if hmmer_predictions[0] != "-":  # CAZyme was predicted

        for domain in hmmer_predictions:
            # separate the predicted CAZy family/subfamily from the domain range
            domain = domain.split("(")

            # Get predicted CAZy family and subfamily if given, called domain_name
            cazy_family = domain[0].split("_")[0]  # drop CAZy subfamily

            # Check for none standard CAZy family naming, and standardise
            # known culprits are GT2_Chitin.. and GT2_Glyco.. though there may be others
            try:
                non_standard_name = re.match(r"\D{2,3}\d+_\D", cazy_family).group()
                # log the non-standard name incase the user wants it later
                cazy_family = non_standard_name[:non_standard_name.find("_")]
            except AttributeError:
                pass

            if cazy_family.startswith(("G", "P", "C", "A")) is False:  # filter out EC numbers
                continue

            # Get domain range
            # Get predicted AA range for domain, called domain_range
            if len(domain) == 1:  # if no domain predicted then no domain ranges predicted
                continue

            else:
                domain_range = domain[1]
                try:
                    re.match(r"\d+-\d+\)", domain_range).group()
                    # remove terminal ')' and standardise range separation
                    domain_range = domain_range[:-1].replace("-", "..")
                except AttributeError:
                    domain_range = None
                    logger.warning(
                        "Non-standardised domain range in HMMER for "
                        f"{tool_data}, {domain_range}"
                        f"- Returning no domain range for {domain_range}"
                    )

            # check family format
            try:
                re.match(r"\D{2,3}\d+", cazy_family).group()
                
                try:
                    fam_annotations[cazy_family].append(domain_range)
                except KeyError:
                    fam_annotations[cazy_family] = [domain_range]

            except AttributeError:
                logger.warning(
                    "Non-standardised CAZy family name in HMMER for "
                    f"{tool_data}, {cazy_family}"
                    f" - Returning no CAZy family."
                )

    return fam_annotations


def parse_hotpep_output(tool_data, logger):
    """Parse the output from Hotpep from the overview.txt file.

    Retrieve the protein accession, predicated CAZy families, accomanying subfamiles.
    Multiple domains can be predicated, and the families and subfamiles will be collected as
    as lists, with the index for each item in each list being associated with the matching items
    (by index) in the other lists. For example, the first item in every list will contain the
    prediciton for the first predicated domain, the second item the second prediction etc.

    Two dataframes are produced in parallel. The first is to build a human readble dataframe for
    writing out to a .csv file. The second stores the data in lists to build the dbCAN consensus
    dataframe later on.

    :param tool_data: str, data from Hotpep cell in dbCAN df
    :param logger: logger object

    Return dict {fam: [domain ranges]}
    """
    # Separate out each of the predicted CAZy families/subfamiles
    hotpep_predictions = tool_data[2].split("+")

    fam_annotations = {}  # {fam: [domain ranges]}

    # remove the k-mer cluster number
    for domain in hotpep_predictions:

        domain = domain.split("(")[0]  # drop the k-mer group number

        if domain == "-":  # no CAZy family predicted
            continue

        if domain.startswith(("G", "P", "C", "A")):  # filter out EC numbers
            cutoff = domain.find("_")
            cazy_family = domain[:cutoff].split("_")[0]  # drop CAZy subfamily

            try:
                re.match(r"\D{2,3}\d+", cazy_family).group()
                fam_annotations[cazy_family] = [None]
            except AttributeError:
                logger.warning(
                        "Non-standardised CAZy family name in overview.txt for Hotpep for "
                        f"{tool_data}, {cazy_family}"
                        f" - Returning no CAZy family for {cazy_family}"
                    )

    return fam_annotations


def parse_diamond_output(tool_data, logger):
    """Parse the output from DIAMOND from the overview.txt file.

    Retrieve the protein accession, predicated CAZy families, accomanying subfamiles.
    Multiple domains can be predicated, and the families and subfamiles will be collected as
    as lists, with the index for each item in each list being associated with the matching items
    (by index) in the other lists. For example, the first item in every list will contain the
    prediciton for the first predicated domain, the second item the second prediction etc.

    Two dataframes are produced in parallel. The first is to build a human readble dataframe for
    writing out to a .csv file. The second stores the data in lists to build the dbCAN consensus
    dataframe later on.

    :param tool_data: str, contents from HMMER cell in dbCAN df
    :param logger: logger object

    Return dict {fam: [domain ranges]}
    """
    # Separate out each of the CAZy family/subfamily predictions from DIAMOND
    diamond_prediction = tool_data.split('+')

    fam_annotations = {}  # {fam: [domain ranges]}

    for domain in diamond_prediction:
        if domain == '-':  # no CAZy family of subfamily predicted
            continue

        # separate predicted family and subfamily
        cutoff = domain.find("_")
        cazy_family = domain[:cutoff]  # drop subfamily
        
        # check family and subfamily formats are correct
        try:
            re.match(r"\D{2,3}\d+", cazy_family).group()
            if cazy_family.startswith(("G", "P", "C", "A")):  # skip EC numbers
                fam_annotations[cazy_family] = [None]
        except AttributeError:
            logger.warning(
                    "Non-standardised CAZy family name in overview.txt for DIAMOND for "
                    f"{tool_data}, {cazy_family}"
                    f" - Returning no CAZy family for {cazy_family}"
                )

    return fam_annotations


def get_dbcan_consensus(hmmer, hotpep, diamond):
    """Get consensus results across HMMER, Hotpep and DIAMOND.

    A consensus result is defined by a result at at least two of the tools in dbCAN have predicated.
    Stores the consensus results in a list, with no duplicates.

    Retrieves list of items common to all three tools, and another list of items common to two of
    the tools. Then builds a dataframe storing the data.

    :param hmmer: list of predictions from HMMER
    :param hotpep: list of predictions from Hotpep
    :param diamond: list of predictions from DIAMOND

    Returns a list.
    """
    # Retrieve list of items predicated by all three tools
    consensus = list(set(hmmer) & set(hotpep) & set(diamond))

    # add domains predicated by two of the tools
    consensus += list(set(hmmer) & set(hotpep))
    consensus += list(set(hmmer) & set(diamond))
    consensus += list(set(hotpep) & set(diamond))

    # remove duplicates
    return set(consensus)
