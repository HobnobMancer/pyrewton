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
import numpy as np


# Functions for parsing the output from dbCAN, producing a dataframe each for the consensus dbCAN
# result, HMMER, Hotpep and DIAMOND

def parse_dbcan_output(overview_file_path):
    """Parse the output from the run_dbCAN overview.txt file, writing a dataframe per tool.

    :param overview_file_path: path, path to the output 'overview.txt' file

    Return 4 dataframes, one each for consensus dbCAN result, HMMER, Hotpep and DIAMOND
    """
    with open(overview_file_path, "r") as fh:
        overview_file = fh.read().splitlines()

    # build empty dataframes to add the repsective prediction tools output to
    hmmer_df = pd.DataFrame({}, columns=[
        "protein_accession",
        "cazy_family",
        "cazy_subfamily",
        "domain_range"
    ])
    hotpep_df = pd.DataFrame({}, columns=["protein_accession", "cazy_family", "cazy_subfamily"])
    diamond_df = pd.DataFrame({}, columns=["protein_accession", "cazy_family", "cazy_subfamily"])
    dbcan_df = pd.DataFrame({}, columns=["protein_accession", "cazy_family", "cazy_subfamily"])

    for line in overview_file[1:]:  # skip the first line becuase this is the head titles
        line = line.split("\t")

        # retrieve the data for HMMER
        hmmer_temp_df = parse_hmmer_output(line)
        hmmer_df = hmmer_df.append(hmmer_temp_df, ignore_index=True)

        # retrieve the data for Hotpep
        hotpep_temp_df = parse_hotpep_output(line)
        hotpep_df = hotpep_df.append(hotpep_temp_df, ignore_index=True)

        # retrieve the data for DIAMOND
        diamond_temp_df = parse_diamond_output(line)
        diamond_df = diamond_df.append(diamond_temp_df, ignore_index=True)

        # retrieve consensus results for dbCAN
        if line[-1] != "1":   # check #ofTools, same approach as other tools:
            # Do not include if non-CAZyme prediction
            consensus_cazy_fam = get_dbcan_consensus(
                hmmer_temp_df.iloc[0, 1],
                hotpep_temp_df.iloc[0, 1],
                diamond_temp_df.iloc[0, 1],
            )
            consensus_sub_fam = get_dbcan_consensus(
                hmmer_temp_df.iloc[0, 2],
                hotpep_temp_df.iloc[0, 2],
                diamond_temp_df.iloc[0, 2],
            )

            consensus_dict = {
                "protein_accession": [line[0]],
                "cazy_family": [consensus_cazy_fam],
                "cazy_subfamily": [consensus_sub_fam],
            }

            temp_consensus_df = pd.DataFrame(consensus_dict)

            dbcan_df = dbcan_df.append(temp_consensus_df, ignore_index=True)

    return dbcan_df, hmmer_df, hotpep_df, diamond_df


def add_hotpep_ec_predictions(hotpep_output_file, hotpep_df):
    """Retrieve predicated EC numbers from Hotpep output file and add to the Hotpep dataframe.

    :param hotpep_output_file: path, path to Hotpep.out file
    :param hotpep_df: pandas dataframe, containing Hotpep predicated CAZy families

    Return pandas dataframe (hotpep_df with EC number predictions)
    """
    with open(hotpep_output_file, "r") as fh:
        hotpep_file = fh.read().splitlines()

        ec_predictions = {}

        for line in hotpep_file[1:]:  # skip the first line which contains the titles
            line = line.split("\t")

            ec_numbers = line[-1].split(",")  # multiple EC numbers may be predicted

            # remove (":score") from each EC number, and convert "NA" to proper null value
            index = 0
            for index in range(len(ec_numbers)):
                ec = ec_numbers[index].split(":")[0]  # remove the (":score") from the EC number

                if ec == "NA":
                    ec_numbers[index] = np.nan

                else:
                    ec_numbers[index] = ec

            protein_accession = line[2]
            ec_predictions[protein_accession] = ec_numbers

    for protein_accession in ec_predictions:  # dict {protein_accession: [predicated EC#s]}
        # Check if there are EC numbers to add
        if type(ec_predictions[protein_accession][0]) == float:  # EC number is null value
            continue

        # find index of row in hotpep_df with the same protein accession
        row_index = hotpep_df.index[hotpep_df["protein_accession"] == protein_accession].tolist()[0]
        # add EC numbers to this row
        string = ""
        for item in ec_predictions[protein_accession][:-1]:
            string += f"{item},"
        string += ec_predictions[protein_accession][((len(ec_predictions[protein_accession])) - 1)]
        hotpep_df.iloc[row_index, 3] = string

    return hotpep_df


def parse_hmmer_output(line):
    """Parse the output from HMMER from the overview.txt file.

    Retrieve the protein accession, predicated CAZy families, accomanying subfamiles and domain
    ranges. Multiple domains can be predicated, and the families, subfamiles and ranges will be
    collected as as lists, with the index for each item in each list being associated with the
    matching items (by index) in the other lists. For example, the first item in every list will
    contain the prediciton for the first predicated domain, the second item the second prediction
    etc.

    :param line: str, line from the dbcan overview.txt file

    Return pandas dataframe."""
    # Retrieve predictions from HMMER (predicated domains and domain ranges)
    # separate out each of the domain predications
    hmmer_prediction = line[1].split("+")

    # create empty lists to store predicated domains and ranges
    cazy_fam = []
    cazy_subfam = []
    domain_ranges = []

    # separate out the name of the predicated domain (the CAZy family) and the domain's AA range
    for domain in hmmer_prediction:
        domain = domain.split("(")
        # standardise name if necessary, retrieve subfam if preciated, and produce null values
        domain_name = domain[0]
        if domain_name.startswith("GT2_Chitin_"):
            cazy_fam.append("GT2")
            cazy_subfam.append(np.nan)

        elif domain_name == '-':
            cazy_fam.append(np.nan)
            cazy_subfam.append(np.nan)

        elif domain_name.find("_") != -1:
            cutoff = domain_name.find("_")
            cazy_fam.append(domain_name[:cutoff])
            cazy_subfam.append(domain_name)

        else:
            cazy_fam.append(domain_name)
            cazy_subfam.append(np.nan)

        try:
            if domain[1] == "-":
                domain_ranges.append(np.nan)
            else:
                domain_ranges.append(domain[1][:-1])  # exlude the final ")"
        except IndexError:  # no domain predicated so no domain ranges
            domain_ranges.append(np.nan)

    prediction_dict = {
        "protein_accession": [line[0]],
        "cazy_family": [cazy_fam],
        "cazy_subfamily": [cazy_subfam],
        "domain_range": [domain_ranges],
    }

    prediction_df = pd.DataFrame(prediction_dict)

    return prediction_df


def parse_hotpep_output(line):
    """Parse the output from Hotpep from the overview.txt file.

    Retrieve the protein accession, predicated CAZy families, accomanying subfamiles.
    Multiple domains can be predicated, and the families and subfamiles will be collected as
    as lists, with the index for each item in each list being associated with the matching items
    (by index) in the other lists. For example, the first item in every list will contain the
    prediciton for the first predicated domain, the second item the second prediction etc.

    :param line: str, line from the dbcan overview.txt file

    Return pandas dataframe."""
    # Retrieve predictions from Hotpep
    hotpep_prediction = line[2].split("+")

    # create empty lists to separate the predicated CAZy fam and subfam
    cazy_fams = []
    cazy_subfams = []

    # remove the k-mer cluster number
    for prediction in hotpep_prediction:
        prediction = prediction.split("(")[0]

        if prediction == "-":
            cazy_fams.append(np.nan)
            cazy_subfams.append(np.nan)

        elif prediction.find("_") != -1:
            cutoff = prediction.find("_")
            cazy_fams.append(prediction[:cutoff])
            cazy_subfams.append(prediction)

        else:
            cazy_fams.append(prediction)
            cazy_subfams.append(np.nan)

    prediction_dict = {
        "protein_accession": [line[0]],
        "cazy_family": [cazy_fams],
        "cazy_subfamily": [cazy_subfams],
    }

    prediction_df = pd.DataFrame(prediction_dict)

    return prediction_df


def parse_diamond_output(line):
    """Parse the output from DIAMOND from the overview.txt file.

    Retrieve the protein accession, predicated CAZy families, accomanying subfamiles.
    Multiple domains can be predicated, and the families and subfamiles will be collected as
    as lists, with the index for each item in each list being associated with the matching items
    (by index) in the other lists. For example, the first item in every list will contain the
    prediciton for the first predicated domain, the second item the second prediction etc.

    :param line: str, line from the dbcan overview.txt file

    Return pandas dataframe."""

    # Retrieve domaind predications from DIAMOND
    diamond_prediction = line[3].split('+')

    # create empty lists to separate the CAZy family and subfamiles
    cazy_fams = []
    cazy_subfams = []

    for pred in diamond_prediction:
        if pred == '-':
            cazy_fams.append(np.nan)
            cazy_subfams.append(np.nan)

        elif pred.find("_") != -1:
            cutoff = pred.find("_")
            cazy_fams.append(pred[:cutoff])
            cazy_subfams.append(pred)

        else:
            cazy_fams.append(pred)
            cazy_subfams.append(np.nan)

    prediction_dict = {
        "protein_accession": [line[0]],
        "cazy_family": [cazy_fams],
        "cazy_subfamily": [cazy_subfams],
    }

    prediction_df = pd.DataFrame(prediction_dict)

    return prediction_df


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
    consensus = list(dict.fromkeys(consensus))

    if len(consensus) == 0:
        consensus = [np.nan]

    return consensus


# Functions for parsing the output from CUPP and eCAMI

def parse_cupp_output(log_file_path):
    """Parse the output from the output log file from CUPP and write out data to a dataframe.

    Retrieves the protein accession/name/identifier, predicated CAZy family, predicated CAZy
    subfamily, predicated EC number and predicated range of domain within the protein sequence
    (the index of the first and last residues of the domain).

    :param log_file_path: path, path to the output log file

    Return Pandas dataframe containing CUPP output
    """

    with open(log_file_path, "r") as lfh:
        log_file = lfh.read().splitlines()

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
        for item in prediction_output:
            if item.find("..") != -1:
                domain_range = item
                break
            else:
                domain_range = np.nan

        try:
            ec_number = prediction_output[-2].split(":")[1]
        except IndexError:
            ec_number = np.nan

        # build dict to enable easy building of df
        prediction = {
            "protein_accession": [prediction_output[0]],
            "cazy_family": [prediction_output[1]],
            "cazy_subfamily": [prediction_output[-1]],
            "ec_number": [ec_number],
            "domain_range": [domain_range],
        }

        # change empty string to null value
        if len(prediction["cazy_subfamily"][0]) == 0:
            prediction["cazy_subfamily"][0] = np.nan
        # If CUPP suggested subfamily, i.e. a subfamily no in CAZy, mark this with '*'
        else:
            if prediction["cazy_subfamily"][0].find(":") != -1:
                prediction["cazy_subfamily"][0] = "*" + prediction["cazy_subfamily"][0]
                print(prediction["cazy_subfamily"][0])

        prediction_df = pd.DataFrame(prediction)
        cupp_df = cupp_df.append(prediction_df)

    return cupp_df


def parse_ecami_output(txt_file_path):
    """Parse the output from the output text file from eCAMI and write out data to a dataframe.

    Retrieves the protein accession/name/identifier, predicated CAZy family, predicated CAZy
    subfamily, predicated EC number and additional/other domains predicated to also be within the
    protein sequence, indicating prediciton of a multiple module enzymes.

    :param text_file_path: path, path to the output text file

    Return Pandas dataframe containing eCAMI output
    """
    with open(txt_file_path, "r") as fh:
        ecami_file = fh.read().splitlines()

    # build an empty dataframe to add predication outputs to
    ecami_df = pd.DataFrame(columns=[
        "protein_accession",
        "cazy_family",
        "cazy_subfamily",
        "ec_number",
        "domain_range",
    ])

    # parse the outputs so in format suitable for final dataframe
    for line in ecami_file:
        if line.startswith(">"):  # identifies new protein
            prediction_output = line.split("\t")
            cazy_fam = prediction_output[1].split(":")[0]

            cazy_subfam = np.nan
            ec_number = np.nan

            # check for additional domains
            other_domains = prediction_output[2].split("|")
            additional_domains = []
            for domain in other_domains:
                domain = domain.split(":")[0]  # drop the eCAMI group number
                # check if subfamily is predicated for the main CAZy family
                if len(domain.split("_")) != 1:
                    cazy_subfam = domain
                # check if predicated EC number
                elif domain.find(".") != -1:
                    ec_number = domain
                else:
                    additional_domains.append(domain)

            if len(additional_domains) == 0:
                additional_domains = np.nan

            # build dict to enable easy building of df
            prediction = {
                "protein_accession": [prediction_output[0]],
                "cazy_family": [cazy_fam],
                "cazy_subfamily": [cazy_subfam],
                "ec_number": [ec_number],
                "additional_domains": [additional_domains],
            }

            prediction_df = pd.DataFrame(prediction)
            ecami_df = ecami_df.append(prediction_df)

    return ecami_df
