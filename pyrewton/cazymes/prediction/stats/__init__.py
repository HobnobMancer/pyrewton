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

from pyrewton.cazymes.prediction.parse import parse_dbcan_output, parse_cupp_output, parse_ecami_output
from pyrewton.cazymes.prediction.stats import binary_classification, multilabel_classification


class ClassificationDF:
    """Represents a CAZyme/non-CAZyme annotation/classification df for a test set"""

    def __init__(self, genome_id, df):
        self.genome_id = genome_id
        self.df = df
    
    def __str__(self):
        return f"<CAZyme/non-CAZyme classification df for test set {self.genome_id}>"
    
    def __repr__(self):
        return f"<CAZyme/non-CAZyme classification df for test set {self.genome_id}>"


def evaluate_performance(predictions, cazy_dict):
    """Evaluate the performance of the CAZymes prediction tools.

    Binary classification: CAZyme/non-CAZyme differentiation
    Multilabel classification (MLC): CAZy family predictions

    Return nothing.
    """
    all_binary_dfs = []  # used for bootstrapping accuracy
    summary_binary_data = []  # [[stat, genomic assembly, prediction tool, value]]

    mlc_predictions = []  # [[assembly, protein, prediction tool, fam1, fam2, fam3...]]
    mlc_ground_truths = []  # [[assembly, protein, prediction tool, fam1, fam2, fam3...]]

    fam_dict = multilabel_classification.foundation_dict()  # used to form foundation of multilabel classification

    for prediction in tqdm(predictions, desc="Standardising outputs and adding to cumulative df"):
        # standardise the output
        dir_path = prediction.prediction_paths["dir"]
        prediction.prediction_paths["dbcan"] = dir_path / "overview.txt"
        prediction.prediction_paths["hotpep"] = dir_path / "Hotpep.out"
        prediction.prediction_paths["cupp"] = dir_path / "cupp_output.fasta.log"
        prediction.prediction_paths["ecami"] = dir_path / "ecami_output.txt"

        hmmer_predictions, hotpep_predictions, diamond_predictions, dbcan_predictions = parse_dbcan_output.parse_dbcan_output(
            prediction.prediction_paths["dbcan"], prediction.prediction_paths["hotpep"], prediction.fasta,
        )
        cupp_predictions = parse_cupp_output.parse_cupp_output(prediction.prediction_paths["cupp"], prediction.fasta)
        ecami_predictions = parse_ecami_output.parse_ecami_output(prediction.prediction_paths["ecami"], prediction.fasta)

        # BINARY EVALUATION
        print("build binary ground truths  dataframe")
        binary_ground_truths_df = binary_classification.retrieve_ground_truths(dbcan_predictions, cazy_dict)

        tools = [
            [dbcan_predictions, "dbCAN"],
            [hmmer_predictions, "HMMER"],
            [hotpep_predictions, "Hotpep"],
            [diamond_predictions, "DIAMOND"],
            [cupp_predictions, "CUPP"],
            [ecami_predictions, "eCAMI"],
        ]
        # add all predictions from all prediction tools to a single dict
        for tool in tqdm(tools, "adding tools output to cumlative dfs"):
            classifications = binary_classification.get_prediction_classifications(tool[0], tool[1], binary_ground_truths_df)
            mlc_prediction, mlc_ground_truth = multilabel_classification.get_multilable_classifications(tool[0], tool[1], cazy_dict, prediction.source)
            mlc_predictions += mlc_prediction
            mlc_ground_truths += mlc_ground_truth

        classifications_df = binary_classification.build_classification_df(classifications)
        all_binary_dfs.append(ClassificationDF(prediction.source, classifications_df))  # used for bootstrapping accuracy

        summary_lists = binary_classification.evaluate_binary_classification(classifications_df, prediction.source)
        for lst in summary_lists:
            summary_binary_data.append(lst)

    time_stamp = datetime.now().strftime("%Y_%m_%d")

    binary_stat_df = pd.DataFrame(summary_binary_data, columns=['Statistic parameter', 'Genomic assembly', 'Prediction tool', 'Statistic value'])
    binary_stat_df.to_csv(f'binary_classification_evaluation_{time_stamp}.csv')  # USED FOR EVALUATION IN R

    column_names = ["Genomic_accession", "Protein_accession", "Prediction_tool"]
    column_names += list(fam_dict.keys())

    mlc_predictions = pd.DataFrame(mlc_predictions, columns=column_names)
    mlc_predictions.to_csv(f"mlc_evaluation_{time_stamp}.csv")   # USED FOR EVALUATION IN R

    mlc_ground_truths = pd.DataFrame(mlc_ground_truths, columns=column_names)
    mlc_ground_truths.to_csv(f"mlc_ground_truths_{time_stamp}.csv")   # USED FOR EVALUATION IN R

    bootstrap_results = {}
    # boot strap accuracy
    for tool in ["dbCAN", "HMMER", "Hotpep", "DIAMOND", "CUPP", "eCAMI"]:
        bootstrap_results[tool] = {}

        for test_set in tqdm(all_binary_dfs, desc=f"Bootstrapping accuracies for {tool}"):
            # calculate the accuracy
            acc = binary_classification.calc_acc(test_set.df['CAZy'], test_set.df[tool])

            # bootstrap the accuracy
            bs_acc = binary_classification.bootstrap_acc(test_set.df['CAZy'], test_set.df[tool], 1000)

            bootstrap_results[tool][test_set.genome_id] = (acc, bs_acc)
    
    # build bootstrap data into a df that can be parsed by R
    data = []  # [tool, genomic_accessions, acc, lower, median, upper]
    for tool in bootstrap_results:
        for genome in bootstrap_results[tool]:
            new_data = [
                tool,
                genome,
                bootstrap_results[tool][genome][0],
                bootstrap_results[tool][genome][1][0],
                bootstrap_results[tool][genome][1][1],
                bootstrap_results[tool][genome][1][2],
            ]
            data.append(new_data)
    bootstrap_results = pd.DataFrame(data, columns=['Prediction_tool', 'Genome', 'Accuracy', 'Lower', 'Median', 'Upper'])

    bootstrap_results.to_csv(f'bootstrap_accuracy_evaluation_{time_stamp}.csv')  # USED FOR EVALUATION IN R

    return
