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
"""Module for statistical evaluation of CAZyme prediction tools"""


import pandas as pd

from datetime import datetime

from tqdm import tqdm

from pyrewton.cazymes.prediction.parse import (
    parse_dbcan_output,
    parse_cupp_output,
    parse_ecami_output,
)
from pyrewton.cazymes.prediction.stats import (
    binary_cazyme_noncazyme_classification as binary_classification
)
from pyrewton.cazymes.prediction.stats import (
    multilabel_family_classification as multilabel_classification
)


class ClassificationDF:
    """Represents a CAZyme/non-CAZyme annotation/classification df for a test set"""

    def __init__(self, genome_id, df):
        self.genome_id = genome_id
        self.df = df

    def __str__(self):
        return f"<CAZyme/non-CAZyme classification df for test set {self.genome_id}>"

    def __repr__(self):
        return f"<CAZyme/non-CAZyme classification df for test set {self.genome_id}>"


def evaluate_performance(predictions, cazy_dict, args):
    """Evaluate the performance of the CAZymes prediction tools.

    Binary classification: CAZyme/non-CAZyme differentiation
    Multilabel classification (MLC): CAZy family predictions

    :param predictions: list of TestSet class instances
    :param cazy_dict: dict keyed by GenBank protein accession, valued by CAZy family classifications
    :param args: cmd-line args parser

    Return nothing.
    """
    time_stamp = datetime.now().strftime("%Y_%m_%d")

    fam_dict = multilabel_classification.foundation_dict()
    # fam_dict used to form foundation of multilabel classification

    all_binary_dfs, summary_binary_data, mlc_predictions, mlc_ground_truths = parse_predictions(
        predictions, time_stamp
    )

    binary_stat_df = pd.DataFrame(
        summary_binary_data,
        columns=['Statistic parameter', 'Genomic assembly', 'Prediction tool', 'Statistic value'],
    )
    binary_stat_df.to_csv(f'binary_classification_evaluation_{time_stamp}.csv')  # USED FOR EVALUATION IN R

    binary_classification.bootstrap_binary_classifications(all_binary_dfs, args)

    column_names = ["Genomic_accession", "Protein_accession", "Prediction_tool"]
    column_names += list(fam_dict.keys())
    mlc_predictions = pd.DataFrame(mlc_predictions, columns=column_names)
    mlc_predictions.to_csv(f"mlc_evaluation_{time_stamp}.csv")   # USED FOR EVALUATION IN R

    mlc_ground_truths = pd.DataFrame(mlc_ground_truths, columns=column_names)
    mlc_ground_truths.to_csv(f"mlc_ground_truths_{time_stamp}.csv")   # USED FOR EVALUATION IN R

    # evaluate the performance of predicting the correct CAZy family
    multilabel_classification.calc_fam_fbeta_score(
        mlc_predictions,
        mlc_ground_truths,
        time_stamp,
        args,
    )  # USED FOR EVALUATION IN R
    multilabel_classification.calc_fam_fbeta_score_per_testset(
        mlc_predictions,
        mlc_ground_truths,
        time_stamp,
        predictions,
        args,
    )  # USED FOR EVALUATION IN R
    multilabel_classification.calc_fam_stats(
        mlc_predictions,
        mlc_ground_truths,
        time_stamp,
        args,
    )  # USED FOR EVALUATION IN R

    # evaluate the performance of predicting the correct CAZy class
    multilabel_classification.build_class_dataframes(
        mlc_predictions,
        mlc_ground_truths,
        args,
        time_stamp,
    )  # USED FOR EVALUATION IN R

    return


def parse_predictions(predictions, time_stamp, cazy_dict):
    """Parse prediction outputs from prediction tools and create dataframes for the stats eval.

    :param predictions: list of TestSet class instances
    :param time_stamp: str, data and time when evaluation started, used for naming files
    :param cazy_dict: dict keyed by GenBank protein accession, valued by CAZy family classifications

    Return:
    all_binary_dfs - list of pandas dataframes of binary predicitons
    summary_binary_data - dictionary containing all binary predictions
    mlc_predictions - list of CAZy family predictions
    mlc_ground_truths - list of known CAZy family
    """
    all_binary_dfs = []  # used for bootstrapping accuracy
    summary_binary_data = []  # [[stat, genomic assembly, prediction tool, value]]

    mlc_predictions = []  # [[assembly, protein, prediction tool, fam1, fam2, fam3...]]
    mlc_ground_truths = []  # [[assembly, protein, prediction tool, fam1, fam2, fam3...]]

    for prediction in tqdm(predictions, desc="Standardising outputs and adding to cumulative df"):
        # standardise the output
        dir_path = prediction.prediction_paths["dir"]
        prediction.prediction_paths["dbcan"] = dir_path / "overview.txt"
        prediction.prediction_paths["hotpep"] = dir_path / "Hotpep.out"
        prediction.prediction_paths["cupp"] = dir_path / "cupp_output.fasta.log"
        prediction.prediction_paths["ecami"] = dir_path / "ecami_output.txt"

        (
            hmmer_predictions,
            hotpep_predictions,
            diamond_predictions,
            dbcan_predictions,
        ) = parse_dbcan_output.parse_dbcan_output(
            prediction.prediction_paths["dbcan"],
            prediction.prediction_paths["hotpep"],
            prediction.fasta,
        )

        cupp_predictions = parse_cupp_output.parse_cupp_output(
            prediction.prediction_paths["cupp"], prediction.fasta,
        )

        ecami_predictions = parse_ecami_output.parse_ecami_output(
            prediction.prediction_paths["ecami"], prediction.fasta,
        )

        classifications = []  # store all classifications

        tools = [
            [dbcan_predictions, "dbCAN"],
            [hmmer_predictions, "HMMER"],
            [hotpep_predictions, "Hotpep"],
            [diamond_predictions, "DIAMOND"],
            [cupp_predictions, "CUPP"],
            [ecami_predictions, "eCAMI"],
        ]

        # add all predictions from all prediction tools to a single dataframe
        for tool in tqdm(tools, "adding tools output to cumlative dfs"):

            # parse binary predictions
            classifications = binary_classification.get_prediction_classifications(
                tool[0],
                tool[1],
                classifications,
            )

            # parse multilabel classficiation (mlc) of CAZy family predictions
            (
                mlc_prediction, mlc_ground_truth,
            ) = multilabel_classification.get_multilable_classifications(
                tool[0],
                tool[1],
                cazy_dict,
                prediction.source,
            )

            mlc_predictions += mlc_prediction
            mlc_ground_truths += mlc_ground_truth

        # build a datafame of all binary classifications
        classifications_df = binary_classification.build_classification_df(classifications, "DF")

        # add ground truths to classifications.df
        classifications_df = binary_classification.retrieve_ground_truths(
            classifications_df, cazy_dict,
        )

        all_binary_dfs.append(ClassificationDF(prediction.source, classifications_df))
        # all_binary_dfs used for bootstrapping accuracy of binary predictions

        summary_lists = binary_classification.evaluate_binary_classification(
            classifications_df, prediction.source,
        )

        classifications_df.to_csv(f"binary_classifications_{time_stamp}_{prediction.source}.csv")

        for lst in summary_lists:
            summary_binary_data.append(lst)

    return all_binary_dfs, summary_binary_data, mlc_predictions, mlc_ground_truths
