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
"""Contains function to parse the output from eCAMI, and create a standarised output.

:func parse_ecami_output: Coordinate parsing the eCAMI output text file
:func get_ecami_additional_info: Retrieve predicted EC numbers and additional domains
"""

import logging
import re

from Bio import SeqIO
from tqdm import tqdm

from pyrewton.cazymes.evaluate_tools.parse import CazymeDomain, CazymeProteinPrediction


def parse_ecami_output(ecami_output_path, fasta_path):
    """Parse the output from eCAMI, retrieving predicted CAZyme domains, CAZy (sub)families and EC numbers.

    :param ecami_output_path: Path, output text file from eCAMI
    :param fasta_path: Path, FASTA file used as input by eCAMI

    Returns a list of ECAMIprediction instances. One instance per protein in the FASTA file used as eCAMI input"""
    # Read the output file from eCAMI
    with open(ecami_output_path, "r") as fh:
        ecami_file = fh.read().splitlines()

    ecami_predictions = {}  # stores proteins {protein_accession:ECAMIprediction_instance}

    for line in tqdm(ecami_file, desc="Parsing eCAMI output"):
        if not line.startswith('>'):  # only want the first line for a protien not the listed k-mers
            continue

        prediction_output = line.split("\t")

        # retrieve protein accession, removing the '>' prefix and stripping white space
        protein_accession = prediction_output[0][1:].strip()
        protein_accession = protein_accession.split(' ')[0]

        # retrieve the CAZy family
        cazy_family = prediction_output[1].split(":")[0].strip()

        # retrieve the children CAZy subfamily of the CAZy family and EC number annotations
        cazy_subfamilies, ec_numbers = get_subfamily_ec_numbers(prediction_output[2], cazy_family)

        # add ECAMIprediction instance to ecami_prediction
        cazyme_classification = 1  # all proteins included in eCAMI output file are identified as CAZymes

        try:
            # add new predicted CAZyme domain to an existing ECAMIprediction instance
            existing_prediction = ecami_predictions[protein_accession]

            # check if the CAZyme domain has already been passed
            existing_cazyme_domains = existing_prediction.cazyme_domains
            existance = False
            while existance is False:
                for domain in existing_cazyme_domains:
                    if (domain.cazy_family == cazy_family) and (domain.cazy_subfamily == cazy_subfamilies):
                        # Domain has been parsed previously, check if additional EC numbers were retrieved
                        for ec in ec_numbers:
                            if ec not in domain.ec_numbers:
                                domain.ec_numbers.append(ec)
                        existance = True
                # come to the end of existing domains and none match the new CAZyme domain
                break

            if existance is False:
                # create new CAZyme domain
                new_cazyme_domain = CazymeDomain(
                    prediction_tool="eCAMI",
                    protein_accession=protein_accession,
                    cazy_family=cazy_family,
                    cazy_subfamily=cazy_subfamilies,
                    ec_numbers=ec_numbers,
                )
                existing_prediction.cazyme_domains.append(new_cazyme_domain)

        except KeyError:  # raised if no corresponding ECAMIprediction instance was found
            new_cazyme_domain = CazymeDomain(
                prediction_tool="eCAMI",
                protein_accession=protein_accession,
                cazy_family=cazy_family,
                cazy_subfamily=cazy_subfamilies,
                ec_numbers=ec_numbers,
            )

            new_protein = CazymeProteinPrediction(
                prediction_tool="eCAMI",
                protein_accession=protein_accession,
                cazyme_classification=cazyme_classification,
            )

            new_protein.cazyme_domains.append(new_cazyme_domain)

            ecami_predictions[protein_accession] = new_protein

    ecami_predictions = add_non_cazymes(ecami_predictions, fasta_path)

    return list(ecami_predictions.values())


def get_subfamily_ec_numbers(subfam_group, cazy_family):
    """Retrieve the predicted CAZy subfamily and associated EC numbers.

    Retrieves only the child CAZy subfamilies for the CAZy familiy of the current working
    CAZyme domain in the protein. Returns the CAZy subfamilies as a string, and EC numbers
    as a list of string, with each string containing a unique EC number. The 'subfam_group'
    refers to the 'subfam_name_of_the_group:subfam_name_count' in the eCAMI output file.

    :param subfam_group: string, the 'subfam_name_of_the_group:subfam_name_count'
    :param cazy_family: string, CAZy family for the current working CAZyme domain

    Return the CAZy subfamily for the CAZyme domain (str) and list of associated EC numbers.
    """
    cazy_subfamilies = []  # store all listed predicted CAZy subfamilies
    ec_numbers = []  # store all predicted EC numbers

    # individual items are separated by "|"
    subfam_group = subfam_group.split("|")

    # check if the subfam_group contains a predicted CAZy subfamily
    for item in subfam_group:
        # remove the group number of the predicted item
        item = item.split(":")[0]

        # check if the item contains a CAZy subfamily
        try:
            re.match(r"\D{2,3}\d+_\d+", item).group()
            # check if the subfamily belongs to the CAZy family listed for the current CAZyme domain
            if item[:item.find("_")] == cazy_family:
                cazy_subfamilies.append(item)

        except AttributeError:  # raised if the item is not a CAZy subfamily
            pass

        # check if item contains an EC number
        try:
            re.match(r"\d+\.(\d+|-)\.(\d+|-)\.(\d+|-)", item).group()
            ec_numbers.append(item)

        except AttributeError:  # raised if the item is not an EC number
            pass

    cazy_subfamilies = ", ".join(cazy_subfamilies)  # convert to string

    return cazy_subfamilies, ec_numbers


def add_non_cazymes(ecami_predictions, fasta_path):
    """Add non-CAZymes to the parsed eCAMI output.

    The eCAMI output only includes the protein it predicts are CAZymes. Therefore, this function
    goes through the FASTA file that was used as input by eCAMI and adds the proteins not
    classified as CAZymes to the parsed proteins.

    :param ecami_predictions: dict, keyed by protein accession and valued by ECAMIprediction instance
    :param fasta_path: Path, FASTA file used as input by eCAMI

    Return a dictionary keyed by protein accession and valued by corresponding ECAMIprediction instance.
    """
    logger = logging.getLogger(__name__)

    # Add non-CAZymes
    # open the FASTA file containing the input protein sequences
    with open(fasta_path, "r") as fh:
        fasta = fh.read().splitlines()

    testset_proteins = 0

    for line in tqdm(fasta, desc="Adding non-CAZymes"):
        if line.startswith(">"):

            # remove '>' prefix and white space
            protein_accession = line[1:].strip()
            protein_accession = protein_accession.split(' ')[0]

            # check if the protein is already listed in the eCAMI predictions
            try:
                ecami_predictions[protein_accession]

            except KeyError:  # raised of protein not in ecami_predictions
                cazyme_classification = 0
                ecami_predictions[protein_accession] = CazymeProteinPrediction(
                    "eCAMI",
                    protein_accession,
                    cazyme_classification,
                )

    logger.info(
        f"Found {testset_proteins} proteins in test set FASTA:\n{fasta_path}\n"
        f"Parsing eCAMI output found {len(list(ecami_predictions.keys()))}"
    )

    return ecami_predictions
