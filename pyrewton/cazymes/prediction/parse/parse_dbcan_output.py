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


def parse_dbcan_output(overview_file_path, logger):
    """Parse the output from the dbCAN overview.txt file, write a dataframe for each prediction tool

    :param overview_file_path: path, path to the output 'overview.txt' file
    :param logger: logger object

    Return 4 dataframes, one each for consensus dbCAN result, HMMER, Hotpep and DIAMOND
    """
    try:
        with open(overview_file_path, "r") as fh:
            overview_file = fh.read().splitlines()
    except FileNotFoundError:
        logger.error(
            "Could not open the dbCAN overview.txt file\n"
            f"{overview_file_path}\n"
            "Returning no dataframe for dbCAN, HMMER, Hotpep and DIAMOND for this output file"
        )
        return None, None, None, None

    # build empty dataframes to add the repsective prediction tools output to
    # list dfs store data as lists to build the dbcan_df
    hmmer_df = pd.DataFrame({}, columns=[
        "protein_accession",
        "cazyme_classification",
        "cazy_family",
        "cazy_subfamily",
        "domain_range"
    ])
    hotpep_df = pd.DataFrame({}, columns=[
        "protein_accession",
        "cazyme_classification",
        "cazy_family",
        "cazy_subfamily"
    ])
    diamond_df = pd.DataFrame({}, columns=[
        "protein_accession",
        "cazyme_classification",
        "cazy_family",
        "cazy_subfamily",
    ])
    dbcan_df = pd.DataFrame({}, columns=[
        "protein_accession",
        "cazyme_classification",
        "cazy_family",
        "cazy_subfamily",
    ])

    for line in overview_file[1:]:  # skip the first line becuase this is the head titles
        line = line.split("\t")

        # retrieve the data for HMMER
        hmmer_temp_df, hmmer_temp_list_df = parse_hmmer_output(line, logger)
        hmmer_df = hmmer_df.append(hmmer_temp_df, ignore_index=True)

        # retrieve the data for Hotpep
        hotpep_temp_df, hotpep_temp_list_df = parse_hotpep_output(line, logger)
        hotpep_df = hotpep_df.append(hotpep_temp_df, ignore_index=True)

        # retrieve the data for DIAMOND
        diamond_temp_df, diamond_temp_list_df = parse_diamond_output(line, logger)
        diamond_df = diamond_df.append(diamond_temp_df, ignore_index=True)

        # retrieve consensus results for dbCAN
        if line[-1] != "1":   # check ''#ofTools', don't include if non-CAZyme prediction
            consensus_cazy_fam = get_dbcan_consensus(
                hmmer_temp_list_df.iloc[0, 2],
                hotpep_temp_list_df.iloc[0, 2],
                diamond_temp_list_df.iloc[0, 2],
            )
            consensus_sub_fam = get_dbcan_consensus(
                hmmer_temp_list_df.iloc[0, 3],
                hotpep_temp_list_df.iloc[0, 3],
                diamond_temp_list_df.iloc[0, 3],
            )

            if consensus_cazy_fam == "-":  # consensus was the protein was not a CAZyme
                consensus_dict = {
                    "protein_accession": [line[0]],
                    "cazyme_classification": [0],
                    "cazy_family": [consensus_cazy_fam],
                    "cazy_subfamily": [consensus_sub_fam],
                }

            else:  # consensus was that the protein was a CAZyme
                consensus_dict = {
                    "protein_accession": [line[0]],
                    "cazyme_classification": [1],
                    "cazy_family": [consensus_cazy_fam],
                    "cazy_subfamily": [consensus_sub_fam],
                }

            temp_consensus_df = pd.DataFrame(consensus_dict)

            dbcan_df = dbcan_df.append(temp_consensus_df, ignore_index=True)

    return dbcan_df, hmmer_df, hotpep_df, diamond_df


def parse_hmmer_output(line, logger):
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

    :param line: str, line from the dbcan overview.txt file
    :param logger: logger object

    Return pandas dataframe.
    """
    # Retrieve predictions for HMMER and separate out each the predicated domains
    hmmer_predictions = line[1].split("+")

    # create empty lists to store predicated domains and ranges
    cazy_fams = []
    cazy_subfams = []
    domain_ranges = []

    # check if no CAZy prediction was made for protein represented in 'line'
    # if no family/subfamily (and thus domain ranges) predicted result will be '-' in 'line'
    if hmmer_predictions[0] != "-":  # CAZyme was predicted

        for domain in hmmer_predictions:
            # separate the predicted CAZy family/subfamily from the domain range
            domain = domain.split("(")

            # Get predicted CAZy family and subfamily if given, called domain_name
            domain_name = domain[0]
            # if domain_name == '-' no family/subfamily and domain ranges predicted
            # don't do anything, empty list converted to null value further down

            # Check for none standard CAZy family naming, and standardise
            # known culprits are GT2_Chitin.. and GT2_Glyco.. though there may be others
            try:
                non_standard_name = re.match(r"\D{2,3}\d+_\D", domain_name).group()
                # log the non-standard name incase the user wants it later
                domain_name = non_standard_name[:non_standard_name.find("_")]
            except AttributeError:
                pass

            # Check if a family or subfamily was predicted
            if domain_name.find("_") != -1:  # predicted CAZy subfamily
                cutoff = domain_name.find("_")
                fam = domain_name[:cutoff]
                subfam = domain_name
                # check family and subfamily formats are correct
                try:
                    re.match(r"\D{2,3}\d+", fam).group()
                    cazy_fams.append(fam)
                except AttributeError:
                    logger.warning(
                        "Non-standardised CAZy family name in HMMER for "
                        f"{line[0]}, {fam}"
                        f"Returning no CAZy family for {fam}"
                    )

                try:
                    re.match(r"\D{2,3}\d+_\d+", subfam).group()
                    cazy_subfams.append(subfam)
                except AttributeError:
                    logger.warning(
                        "Non-standardised CAZy subfamily name in HMMER for "
                        f"{line[0]}, {subfam}"
                        f"Returning no CAZy subfamily for {subfam}"
                    )

            else:  # only CAZy family predicted
                # check family format
                try:
                    re.match(r"\D{2,3}\d+", domain_name).group()
                    cazy_fams.append(domain_name)
                except AttributeError:
                    logger.warning(
                        "Non-standardised CAZy family name in HMMER for "
                        f"{line[0]}, {domain_name}"
                        f"Returning no CAZy family for {domain_name}"
                    )

            # Get predicted AA range for domain, called domain_range
            if len(domain) == 1:  # if no domain predicted then no domain ranges predicted
                continue

            else:
                domain_range = domain[1]
                try:
                    re.match(r"\d+-\d+\)", domain_range).group()
                    # remove terminal ')' and standardise range separation
                    domain_range = domain_range[:-1].replace("-", "..")
                    domain_ranges.append(domain_range)
                except AttributeError:
                    logger.warning(
                        "Non-standardised domain range in HMMER for "
                        f"{line[0]}, {domain_range}"
                        f"Returning no domain range for {domain_range}"
                    )

        # convert lists to strings for human readability
        # convert empty lists to null values
        # _hr suffix denotes variables for building human readable dataframe
        if len(cazy_fams) == 0:
            cazy_family_hr, cazy_fams = np.nan, [np.nan]
        else:
            cazy_family_hr = ", ".join(cazy_fams)

        if len(cazy_subfams) == 0:
            cazy_subfamily_hr, cazy_subfams = np.nan, [np.nan]
        else:
            cazy_subfamily_hr = ", ".join(cazy_subfams)

        if len(domain_ranges) == 0:
            domain_ranges_hr, domain_ranges = np.nan, [np.nan]
        else:
            domain_ranges_hr = ", ".join(domain_ranges)

        prediction_dict = {
            "protein_accession": [line[0]],
            "cazyme_classification": [1],
            "cazy_family": [cazy_family_hr],
            "cazy_subfamily": [cazy_subfamily_hr],
            "domain_range": [domain_ranges_hr],
        }

        prediction_list_dict = {
            "protein_accession": [line[0]],
            "cazyme_classification": [1],
            "cazy_family": [cazy_fams],
            "cazy_subfamily": [cazy_subfams],
            "domain_range": [domain_ranges],
        }

        prediction_df = pd.DataFrame(prediction_dict)
        prediction_list_df = pd.DataFrame(prediction_list_dict)

        return prediction_df, prediction_list_df

    else:  # HMMER did not predict the protein was a CAZyme
        cazy_family_hr, cazy_fams = np.nan, [np.nan]
        cazy_subfamily_hr, cazy_subfams = np.nan, [np.nan]
        domain_ranges_hr, domain_ranges = np.nan, [np.nan]

        prediction_dict = {
            "protein_accession": [line[0]],
            "cazyme_classification": [0],
            "cazy_family": [cazy_family_hr],
            "cazy_subfamily": [cazy_subfamily_hr],
            "domain_range": [domain_ranges_hr],
        }

        prediction_list_dict = {
            "protein_accession": [line[0]],
            "cazyme_classification": [0],
            "cazy_family": [cazy_fams],
            "cazy_subfamily": [cazy_subfams],
            "domain_range": [domain_ranges],
        }

        prediction_df = pd.DataFrame(prediction_dict)
        prediction_list_df = pd.DataFrame(prediction_list_dict)

        return prediction_df, prediction_list_df



def parse_hotpep_output(line, logger):
    """Parse the output from Hotpep from the overview.txt file.

    Retrieve the protein accession, predicated CAZy families, accomanying subfamiles.
    Multiple domains can be predicated, and the families and subfamiles will be collected as
    as lists, with the index for each item in each list being associated with the matching items
    (by index) in the other lists. For example, the first item in every list will contain the
    prediciton for the first predicated domain, the second item the second prediction etc.

    Two dataframes are produced in parallel. The first is to build a human readble dataframe for
    writing out to a .csv file. The second stores the data in lists to build the dbCAN consensus
    dataframe later on.

    :param line: str, line from the dbcan overview.txt file
    :param logger: logger object

    Return pandas dataframe."""
    # Separate out each of the predicted CAZy families/subfamiles
    hotpep_predictions = line[2].split("+")

    # create empty lists to separate the predicated CAZy fam and subfam
    cazy_fams = []
    cazy_subfams = []

    # remove the k-mer cluster number
    for prediction in hotpep_predictions:

        prediction = prediction.split("(")[0]  # drop the k-mer group number

        if prediction == "-":  # no CAZy family predicted
            cazy_fams_hr, cazy_fams = np.nan, [np.nan]
            cazy_subfams_hr, cazy_subfams = np.nan, [np.nan]

            prediction_dict = {
                "protein_accession": [line[0]],
                "cazyme_classification": [0],
                "cazy_family": [cazy_fams_hr],
                "cazy_subfamily": [cazy_subfams_hr],
            }

            prediction_list_dict = {
                "protein_accession": [line[0]],
                "cazyme_classification": [0],
                "cazy_family": [cazy_fams],
                "cazy_subfamily": [cazy_subfams],
            }

            prediction_df = pd.DataFrame(prediction_dict)
            prediction_list_df = pd.DataFrame(prediction_list_dict)

            return prediction_df, prediction_list_df
        
        else:  # Protein was predicated to be a CAZyme

            # check if a subfamily was predicated
            if prediction.find("_") != -1:  # predicted subfamily
                # separate predicted family and subfamily
                cutoff = prediction.find("_")
                fam = prediction[:cutoff]
                subfam = prediction

                # check family and subfamily formats are correct
                try:
                    re.match(r"\D{2,3}\d+", fam).group()
                    cazy_fams.append(fam)
                except AttributeError:
                    logger.warning(
                            "Non-standardised CAZy family name in overview.txt for Hotpep for "
                            f"{line[0]}, {fam}"
                            f"Returning no CAZy family for {fam}"
                        )

                try:
                    re.match(r"\D{2,3}\d+_\d+", subfam).group()
                    cazy_subfams.append(subfam)
                except AttributeError:
                    logger.warning(
                            "Non-standardised CAZy subfamily name in overview.txt for Hotpep for "
                            f"{line[0]}, {subfam}"
                            f"Returning no CAZy subfamily for {subfam}"
                    )

            else:  # only CAZy family predicted
                try:
                    re.match(r"\D{2,3}\d+", prediction).group()
                    cazy_fams.append(prediction)
                except AttributeError:
                    logger.warning(
                            "Non-standardised CAZy family name in overview.txt for Hotpep for "
                            f"{line[0]}, {prediction}"
                            f"Returning no CAZy family for {prediction}"
                    )
            
    # Convert lists to strings for human readability, identified by _hr suffix
    # Convert lists to null values
    if len(cazy_fams) == 0:
        cazy_fams_hr, cazy_fams = np.nan, [np.nan]
    else:
        cazy_fams_hr = ", ".join(cazy_fams)

    if len(cazy_subfams) == 0:
        cazy_subfams_hr, cazy_subfams = np.nan, [np.nan]
    else:
        cazy_subfams_hr = ", ".join(cazy_subfams)

    prediction_dict = {
        "protein_accession": [line[0]],
        "cazyme_classification": [1],
        "cazy_family": [cazy_fams_hr],
        "cazy_subfamily": [cazy_subfams_hr],
    }

    prediction_list_dict = {
        "protein_accession": [line[0]],
        "cazyme_classification": [1],
        "cazy_family": [cazy_fams],
        "cazy_subfamily": [cazy_subfams],
    }

    prediction_df = pd.DataFrame(prediction_dict)
    prediction_list_df = pd.DataFrame(prediction_list_dict)

    return prediction_df, prediction_list_df


def parse_diamond_output(line, logger):
    """Parse the output from DIAMOND from the overview.txt file.

    Retrieve the protein accession, predicated CAZy families, accomanying subfamiles.
    Multiple domains can be predicated, and the families and subfamiles will be collected as
    as lists, with the index for each item in each list being associated with the matching items
    (by index) in the other lists. For example, the first item in every list will contain the
    prediciton for the first predicated domain, the second item the second prediction etc.

    Two dataframes are produced in parallel. The first is to build a human readble dataframe for
    writing out to a .csv file. The second stores the data in lists to build the dbCAN consensus
    dataframe later on.

    :param line: str, line from the dbcan overview.txt file
    :param logger: logger object

    Return pandas dataframe."""

    # Separate out each of the CAZy family/subfamily predictions from DIAMOND
    diamond_prediction = line[3].split('+')

    # create empty lists to separate the CAZy family and subfamiles
    cazy_fams = []
    cazy_subfams = []

    for prediction in diamond_prediction:
        if prediction == '-':  # no CAZy family of subfamily predicted
            cazy_fams_hr, cazy_fams = np.nan, [np.nan]
            cazy_subfams_hr, cazy_subfams = np.nan, [np.nan]

            prediction_dict = {
                "protein_accession": [line[0]],
                "cazyme_classification": [0],
                "cazy_family": [cazy_fams_hr],
                "cazy_subfamily": [cazy_subfams_hr],
            }

            prediction_list_dict = {
                "protein_accession": [line[0]],
                "cazyme_classification": [0],
                "cazy_family": [cazy_fams],
                "cazy_subfamily": [cazy_subfams],
            }

            prediction_df = pd.DataFrame(prediction_dict)
            prediction_list_df = pd.DataFrame(prediction_list_dict)

            return prediction_df, prediction_list_df

        if prediction.find("_") != -1:  # CAZy subfamily predicted
            # separate predicted family and subfamily
            cutoff = prediction.find("_")
            fam = prediction[:cutoff]
            subfam = prediction

            # check family and subfamily formats are correct
            try:
                re.match(r"\D{2,3}\d+", fam).group()
                cazy_fams.append(fam)
            except AttributeError:
                logger.warning(
                        "Non-standardised CAZy family name in overview.txt for DIAMOND for "
                        f"{line[0]}, {fam}"
                        f"Returning no CAZy family for {fam}"
                    )

            try:
                re.match(r"\D{2,3}\d+_\d+", subfam).group()
                cazy_subfams.append(subfam)
            except AttributeError:
                logger.warning(
                        "Non-standardised CAZy subfamily name in overview.txt for DIAMOND for "
                        f"{line[0]}, {subfam}"
                        f"Returning no CAZy subfamily for {subfam}"
                    )

        else:  # only CAZy family predicted
            try:
                re.match(r"\D{2,3}\d+", prediction).group()
                cazy_fams.append(prediction)
            except AttributeError:
                try:
                    # check if an EC number was predicated
                    re.match(r"\d+\.(\d+|-)\.(\d+|-)\.(\d+|-)", prediction).group()
                    ec_numbers.append(prediction)
                except AttributeError:
                    logger.warning(
                            f"Could not determine data type of {line[0]}, {prediction} from "
                            "DIAMOND\nIn the overview.txt file.\n"
                            f"Returning no CAZy family for {prediction}"
                        )

    # Convert lists to strings for human readability, identified by _hr suffix
    # Convert lists to null values
    if len(cazy_fams) == 0:
        cazy_fams_hr, cazy_fams = np.nan, [np.nan]
    else:
        cazy_fams_hr = ", ".join(cazy_fams)

    if len(cazy_subfams) == 0:
        cazy_subfams_hr, cazy_subfams = np.nan, [np.nan]
    else:
        cazy_subfams_hr = ", ".join(cazy_subfams)

    prediction_dict = {
        "protein_accession": [line[0]],
        "cazyme_classification": [1],
        "cazy_family": [cazy_fams_hr],
        "cazy_subfamily": [cazy_subfams_hr],
    }

    prediction_list_dict = {
        "protein_accession": [line[0]],
        "cazyme_classification": [1],
        "cazy_family": [cazy_fams],
        "cazy_subfamily": [cazy_subfams],
    }

    prediction_df = pd.DataFrame(prediction_dict)
    prediction_list_df = pd.DataFrame(prediction_list_dict)

    return prediction_df, prediction_list_df


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

    # convert to strings for human readability or null value for empty list
    if (len(consensus) == 0) or (type(consensus[0]) == float):  # consensus was np.nan
        consensus = np.nan
    else:
        consensus = ", ".join(consensus)

    return consensus


def add_hotpep_ec_predictions(hotpep_output_file, hotpep_df, logger):
    """Retrieve predicated EC numbers from Hotpep output file and add to the Hotpep dataframe.

    :param hotpep_output_file: path, path to Hotpep.out file
    :param hotpep_df: pandas dataframe, containing Hotpep predicated CAZy families
    :param logger: logger object

    Return pandas dataframe (hotpep_df with EC number predictions)
    """
    # Add a column of null values to which EC numbers will be added
    hotpep_df["ec_number"] = np.nan

    # Retrieve all the predicted EC numbers from the Hotpep.out file
    try:
        with open(hotpep_output_file, "r") as fh:
            hotpep_file = fh.read().splitlines()
    except FileNotFoundError:
        logger.error(
            "Could no open Hotpep output file\n"
            f"{hotpep_output_file}"
            "Retunring no EC numbers for this file"
        )
        return hotpep_df

    # create empty dict to store EC numbers, keyed by protein accession, valued by predicted EC#
    ec_predictions = {}

    for line in hotpep_file[1:]:  # skip the first line which contains the titles
        # separate items in the line
        line = line.split("\t")

        # separate out the predicted EC numbers, multiple maybe predicted
        ec_numbers = line[-1].split(", ")

        # remove (":score") from each EC number, and convert "NA" to proper null value
        index = 0
        for index in range(len(ec_numbers)):
            ec = ec_numbers[index].split(":")[0].strip()  # remove the (":score") from the EC#

            if ec == "NA":  # used by Hotpep to show no predicated EC number
                ec_numbers[index] = np.nan

            else:  # Hotpep infers it is an EC number
                # check EC number is formated correctly
                try:
                    re.match(r"\d+?\.(\d+?|\-)\.(\d+?|\-)\.(\d+?|\-)", ec).group()
                    ec_numbers[index] = ec
                except AttributeError:
                    logger.warning(
                            "Non-standardised EC# formate in Hotpep.out for "
                            f"{line[2]}, {ec} in\n"
                            f"{hotpep_output_file}\n"
                            f"Returning no EC number for {ec}"
                        )

        protein_accession = line[2]
        ec_predictions[protein_accession] = ec_numbers

    # Add the predicted EC numbers for their respective protein in the hotpep_df (hotpep dataframe)

    for protein_accession in ec_predictions:  # dict {protein_accession: [predicated EC#s]}
        # Check if there are EC numbers to add
        if type(ec_predictions[protein_accession][0]) == float:  # EC numbers is stored as np.nan
            continue

        # find index of row in hotpep_df with the same protein accession
        row_index = hotpep_df.index[hotpep_df["protein_accession"] == protein_accession].tolist()[0]

        # compile all predicted EC numbers for current protein accession into a single string
        try:
            string = ", ".join(ec_predictions[protein_accession])
        except TypeError:  # raised in np.nan is in predicted EC list
            string = ""
            # add all but last EC number to string
            for item in ec_predictions[protein_accession][:-1]:
                if type(item) != float:  # if item is not np.nan
                    string += f"{item}, "
            # add the last predicted EC number
            string += ec_predictions[protein_accession][-1]

        # Add the string of predicted EC numbers to the dataframe
        hotpep_df.iloc[row_index, 3] = string

    return hotpep_df
