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
"""Evaluate accuracy of dbCAN, CUPP and eCAMI prediction of CAZyme class and family.

:cmd_args:

:func --:
"""

import sys

import pandas as pd

from sklearn.metrics import adjusted_rand_score


def main():
    """Setup parser, logger and hyperparameters for ARI."""

    # Set up parser

    # Build logger

    # Set parameters for ARI
    X_dbcan, X_diamond, X_hmmer, X_hotpep, X_cupp, X_ecami = get_matrix_X(
        args, logger
    )  # contains dataset
    y = get_feature_array_y  # contains known lables (truths)

    # Calculate Rand index per CAZyme prediction tool
    rand_indexes = []  # empty tuple to store calculated results
    rand_indexes.append(calculate_ri(X_dbcan, y, logger))
    rand_indexes.append(calculate_ri(X_cupp, y, logger))
    rand_indexes.append(calculate_ri(X_ecami, y, logger))

    # Generate dataframe to display results
    result_df = pd.DataFrame(
        rand_indexes, columns=["dbCAN", "DIAMOND", "HMMER", "Hotpep", "CUPP", "eCAMI"]
    )

    # Write out results to file if enabled
    if args.output is not sys.stdout:
        write_output(result_df, args, logger)


def get_matrix_X(args, logger):
    result_df = calculate_ri(X, y, logger)  # dataframe containing Rand indexes
    """Generate one matrix containing predicated CAZyme classes per prediciton tool.

    :param args: parser object
    :param logger: logger object

    Returns pandas dataframe.
    """
    args.dbcan
    args.cupp
    args.ecami
    # args contains paths to the output files
    # parse each output file as required
    # create matrix for each output file containing a unique protein per row
    # and each row containing protein ID, host taxonomy ID and predicated CAZy family
    # sort each dataframe by protein ID (ensure will be in same order as feature array)
    # Create a single arrayF containing all dbCAN results where two or more tools agreed
    # Create a single matrix per dbcan tool (diamond, hmmer, hotpep, ecami)
    return X_dbcan, X_diamond, X_hmmer, X_hotpep, X_cupp, X_ecami


def get_feature_array_y(args, logger):
    """Generate one dimenstrional array of none CAZyme class/families from CAZy.

    :param args: parser object
    :param logger: logger object

    Return list.
    """
    args.cazy
    # args.cazy contains path to dataframe containing proteins IDs and CAZy
    # open dataframe
    # sort by protein ID
    # parse the dataframe to create a list containing the CAZy families
    # return list of CAZy families in same order as proteins in matrix X
    return cazy_list


def calculate_ri(X, y, logger):
    """Calculate the rand index.

    :param X: dataframe, ARI matrix X
    :param y: list, ARI feature array
    :param logger: logger object

    Return list of length 1.
    """


def write_output(result, args, logger):
    """Write out calculated Rand indexes to file"""
