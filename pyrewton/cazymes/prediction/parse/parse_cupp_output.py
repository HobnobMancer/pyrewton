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
"""Contains function to parse the output from CUPP. Create standarised output.

:func parse_cupp_output: Coordinate parsing the CUPP output log file
:func get_cupp_domain_range: Retrieve the domain ranges of predicted domains
:func get_cupp_ec_number: Retrieve the predicted EC numbers of predicted domains
"""

import logging
import re

import numpy as np

from tqdm import tqdm

from pyrewton.cazymes.prediction.parse import CazymeDomain, CazymeProteinPrediction


def parse_cupp_output(log_file_path, fasta_path):
    """Parse the output log file from CUPP and write out a dataframe.

    Retrieves the protein accession, predicted CAZy families, CAZy families EC
    numbers and domain ranges.

    :param log_file_path: Path, path to the output log file
    :fasta path: Path, path to the FASTA containing the query sequences for CUPP

    Return a list of CUPPprediction instances, where each instance represented one protein.
    """
    # open the CUPP output log file
    with open(log_file_path, "r") as lfh:
        log_file = lfh.read().splitlines()

    cupp_predictions = {}  # stores proteins {protein_accession:CUPPprediction_instance}

    for line in tqdm(log_file, desc="parsing cupp output"):
        # separate the data fields
        prediction_output = line.split("\t")

        # retrieve the predicted domain range if given
        domain_range = get_cupp_domain_range(prediction_output)  # list of strs

        # retrieve predicted EC number if given
        ec_numbers = get_cupp_ec_number(prediction_output)  # list of strs

        # retrieve predicted CAZy subfamily if given
        cazy_subfamily = prediction_output[-1] 
        # if no CAZy subfamily was predicted, change empty string to null value
        if len(cazy_subfamily) == 0:
            cazy_subfamily = np.nan
        else:
            # write subfamilies using CAZy standard format
            cazy_subfamily = cazy_subfamily.replace(":","_")
            cazy_subfamily = cazy_subfamily.replace("+", ", ")  # separate multiple subfams with ', ' not '+'

        # retrieve the protein accession and check if it already has a corresponding CUPPprediction instance
        protein_accession = prediction_output[0]

        # retrieve the CAZyme classification and predicted CAZy family
        cazyme_classification = 1  # all proteins included in the file are identified as CAZymes
        cazy_family = prediction_output[1]  # This is selected from the CUPP 'Best function'

        # check if a CUPPprediction instance already exists for the protein
        try:
            existing_prediction = cupp_predictions[protein_accession]

            # check if the CAZyme domain has been been parsed before
            existing_cazyme_domains = existing_prediction.cazyme_domains

            existance = False
            for domain in existing_cazyme_domains:
                if (domain.cazy_family == cazy_family) and (domain.cazy_subfamily == cazy_subfamily):
                    for ec in ec_numbers:
                        domain.ec_numbers.append(ec)
                    for drange in domain_range:
                        domain.domain_range.append(drange)
                    existance = True

            if existance is False:
                # create new CAZyme domain
                new_cazyme_domain = CazymeDomain(
                    protein_accession,
                    cazy_family,
                    cazy_subfamily,
                    ec_numbers,
                    domain_range,
                )
                existing_prediction.cazyme_domains.append(new_cazyme_domain)

        except KeyError:  # raised if there is not instance for the protein

            new_cazyme_domain = CazymeDomain(
                protein_accession,
                cazy_family,
                cazy_subfamily,
                ec_numbers,
                domain_range,
            )

            new_protein = CazymeProteinPrediction(
                "CUPP",
                protein_accession,
                cazyme_classification,
                [new_cazyme_domain],
            )

            cupp_predictions[protein_accession] = new_protein

    # add non-CAZymes
    cupp_predictions = add_non_cazymes(fasta_path, cupp_predictions)

    cupp_predictions = list(cupp_predictions.values())

    return cupp_predictions


def get_cupp_domain_range(prediction_output):
    """Retrieve the predicted amino acid range of predicted CAZyme domain from CUPP log file.

    :param prediction_output: list of items from log file line.

    List of predicted domain ranges or null value in a list if not.
    """
    logger = logging.getLogger(__name__)

    domain_range = []  # store as a list in case multiple domain ranges are given

    for item in prediction_output:
        if item.find("..") != -1:
            domain_range.append(item)

    # check retrieved items are definetly the domain ranges
    for item in domain_range:
        try:
            re.match(r"\d+\.\.\d+", item).group()
        except AttributeError:
            # write as logger in pyrewton
            logger.warning("{item} misidentified as domain range")
            domain_range.remove(item)

    if len(domain_range) == 0:
        domain_range = [np.nan]

    return domain_range


def get_cupp_ec_number(prediction_output):
    """Retrieve predicted EC numbers from CUPP log file.

    EC numbers are represented as "CAZy_fam:EC_number" in the log file.
    If multiple EC numbers are predicated and the CAZy families of each are
    given then these are separated by '-'. If multiple EC numbers are predicated
    for the same CAZy family, these are separated by '&'.

    :param prediciton_output: list of items from log file line.

    List of predicted EC numbers or null value in a list if not.
    """
    # retrieve the data from the CUPP 'Best function' prediction
    ec_data = prediction_output[-2].split(":")

    if len(ec_data) == 1:  # Result of ":" not being present, caused by no EC number predictions
        ec_number = np.nan
        return ec_number

    # create empty list to store all predicted EC numbers
    all_ec_numbers = []

    for item in ec_data:
        # check if the item may be an EC number, which start with a digit
        try:
            re.match(r"\d.+", item).group()
        except AttributeError:  # not an EC number
            continue

        # if multiple best function and/or EC number predictions were made they may
        # be separated by a dash '-'
        dash_separated_data = item.split("-")

        for data in dash_separated_data:
            # if multiple EC numbers are predicated they are separate by '&
            split_data = data.split("&")

            for string in split_data:
                # check if the string is an EC number
                try:
                    re.match(r"\d+?\.(\d+?|\*)\.(\d+?|\*)\.(\d+?|\*)", string). group()
                    all_ec_numbers.append(string)
                except AttributeError:  # not an EC number
                    continue

    # standardise missing digits in the EC numbers from '*' to '-'                
    index = 0
    for index in range(len(all_ec_numbers)):
        all_ec_numbers[index] = all_ec_numbers[index].replace("*","-")
        # Sometimes 'Unknown' is written out by CUPP, ensure this is removed
        if all_ec_numbers[index].find("Unknown") != -1:
            all_ec_numbers.remove(all_ec_numbers[index])

    # write out the EC numbers in easy human readable format or create null value if none were predicted
    if len(all_ec_numbers) == 0:
        all_ec_numbers = [np.nan]

    return all_ec_numbers


def add_non_cazymes(fasta_path, cupp_predictions):
    """Add proteins that CUPP identified as non-CAZymes to the collection of CUPPprediction instances.

    :param fasta_path: Path, FASTA file used as input for CUPP.
    :param cupp_predictions: dict, key=protein_accession, value=CUPPprediction instance

    Return a dictionary valued by protein accessions and keyed by their respective CUPPprediction instance.
    """
    # open the FASTA path
    with open(fasta_path) as fh:
        fasta = fh.read().splitlines()

    for line in tqdm(fasta, desc="Adding non-CAZymes"):
        if line.startswith(">"):
            protein_accession = line[1:].strip()

            # check if the protein has been listed as a CAZyme by CUPP
            try:
                cupp_predictions[protein_accession]
            except KeyError:
                # raised if protein not in cupp_predictions, inferring proten was not labelled as CAZyme
                cazyme_classification = 0
                cupp_predictions[protein_accession] = CazymeProteinPrediction(
                    "CUPP",
                    protein_accession,
                    cazyme_classification,
                )

    return cupp_predictions
