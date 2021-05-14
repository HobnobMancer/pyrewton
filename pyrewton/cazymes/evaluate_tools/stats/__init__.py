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
"""Module for statistical evaluation of CAZyme prediction tools

:class ClassificationDF: Represents a CAZyme/non-CAZyme annotation/classification df for a test set

:func evaluate_performance: Co-ordinates overall evaluation of prediction tools
:func parse_predictions: Parses the raw output files from the CAZyme prediction tools

"""


import pandas as pd

from datetime import datetime

from tqdm import tqdm

from pyrewton.cazymes.evaluate_tools.parse import (
    parse_dbcan_output,
    parse_cupp_output,
    parse_ecami_output,
)
from pyrewton.cazymes.evaluate_tools.stats import (
    binary_cazyme_noncazyme_classification as binary_classification
)
from pyrewton.cazymes.evaluate_tools.stats import (
    multilabel_family_classification as family_classifications
)
from pyrewton.cazymes.evaluate_tools.stats import (
    multilabel_class_classifications as class_classifications
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

    # build dataframes of predictions
    (
        all_binary_c_nc_dfs,
        binary_c_nc_statistics,
        class_predictions_df,
        class_ground_truths_df,
        all_family_predictions,
        all_family_ground_truths,
    ) = build_prediction_dataframes(
        predictions,
        time_stamp,
        cazy_dict,
    )

    binary_c_nc_statistics.to_csv(f'binary_classification_evaluation_{time_stamp}.csv')  # USED FOR EVALUATION IN R

    # bootstrap resample binary CAZyme/non-CAZyme classifications, to evaluate performance range
    binary_classification.bootstrap_binary_c_nc_classifications(all_binary_c_nc_dfs, args)

    # write out ground truth CAZy class annotations to disk
    all_class_ground_truths.to_csv(f"cazy_class_ground_truths_{time_stamp}.csv")

    # write out ground truth CAZy family annotations to disk
    all_family_ground_truths.to_csv(f"cazy_family_ground_truths_{time_stamp}.csv")

    # Evaluate the performance of CAZy class predictions
    # Calculate ARI and RI for multilabel evaluation, and add to CAZy class prediction dataframe
    class_predictions_df = class_classifications.calculate_class_ari_ri(
        class_ground_truths_df,
        class_predictions_df,
        time_stamp)
    # write out predicted CAZy class annotations and RI and ARI to disk
    class_predictions_df.to_csv(f"cazy_class_predictions_{time_stamp}.csv")   # USED FOR EVALUATION IN R

    # Calculate Fbeta, Sens, Spec and Acc for predicting each CAZy class
    # calculate the fbeta_score, sensitivity(recall), specificity and accuracy per CAZy class
    calculate_class_stats(class_ground_truths_df, class_predictions_df, time_stamp, args)

    calculate_class_stats_by_testsets(
        class_ground_truths_df,
        class_predictions_df,
        time_stamp,
        args,
    )

    # Evaluate the performance of CAZy Family predictions
    # Calculate ARI and RI for multilabel evaluation, and add to CAZy family prediction dataframe
    class_predictions_df = calculate_family_ari(
        class_ground_truths_df,
        class_predictions_df,
        time_stamp,
    )
    # write out predicted CAZy class annotations and RI and ARI to disk
    all_class_predictions.to_csv(f"cazy_class_predictions_{time_stamp}.csv")   # USED FOR EVALUATION IN R

    # Calculate Fbeta, Sens, Spec and Acc for predicting each CAZy class
    # calculate the fbeta_score, sensitivity(recall), specificity and accuracy per CAZy class
    calculate_class_stats(class_ground_truths_df, class_predictions_df, time_stamp, args)

    calculate_class_stats_by_testsets(
        class_ground_truths_df,
        class_predictions_df,
        time_stamp,
        args,
    )




    # evaluate the performance of predicting the correct CAZy family
    family_classifications.calc_fam_fbeta_score(
        all_family_predictions,
        all_family_ground_truths,
        time_stamp,
        args,
    )
      # USED FOR EVALUATION IN R
    family_classifications.calc_fam_fbeta_score_per_testset(
        all_family_predictions,
        all_family_ground_truths,
        time_stamp,
        predictions,
        args,
    )
      # USED FOR EVALUATION IN R
    family_classifications.calc_fam_stats(
        all_family_predictions,
        all_family_ground_truths,
        time_stamp,
        args,
    )
      # USED FOR EVALUATION IN R

    # evaluate the performance of predicting the correct CAZy class

      # USED FOR EVALUATION IN R

    return


def build_prediction_dataframes(predictions, time_stamp, cazy_dict):
    """Parse prediction outputs from prediction tools and create dataframes for the stats eval.

    :param predictions: list of TestSet class instances
    :param time_stamp: str, data and time when evaluation started, used for naming files
    :param cazy_dict: dict keyed by GenBank protein accession, valued by CAZy family classifications

    Return:
    all_binary_c_nc_dfs - list of pandas dfs of binary CAZyme/non-CAZyme predicitons, one df is
        created per test set containing CAZmye/non-CAZyme predictions from all prediction tools,
        and the ground truths from CAZy. This is used for bootstrapping the accuracy.
    binary_c_nc_statistics - df of statistical evaluations of binary CAZyme/non-CAZyme predictions
    all_class_predictions - df of predicted CAZy class classifications, including Rand Index
        and Adjusted Rand Index for every protein, except true negative non-CAZyme predictions
    all_class_ground_truths - df of CAZy determined CAZy class classifications
    all_family_predictions - 
    all_family_ground_truths - 
    """
    all_binary_c_nc_dfs = []  # used for bootstrapping accuracy
    binary_c_nc_statistics = []  # [[stat, genomic assembly, prediction tool, value]]

    all_family_predictions = []  # [[assembly, protein, prediction tool, fam1, fam2, fam3...]]
    all_family_ground_truths = []  # [[assembly, protein, prediction tool, fam1, fam2, fam3...]]

    for test_set in tqdm(predictions, desc="Standardising outputs"):
        # create paths to the raw prediction tool output files to be parsed
        dir_path = test_set.prediction_paths["dir"]
        test_set.prediction_paths["dbcan"] = dir_path / "overview.txt"
        test_set.prediction_paths["hotpep"] = dir_path / "Hotpep.out"
        test_set.prediction_paths["cupp"] = dir_path / "cupp_output.fasta.log"
        test_set.prediction_paths["ecami"] = dir_path / "ecami_output.txt"

        # standardise the raw outputs from the prediction tools
        # tool_predictions is a list of CazymeProteinPrediction instances
        # (see evaluate_tools.parse.__init__.py)
        (
            hmmer_predictions,
            hotpep_predictions,
            diamond_predictions,
            dbcan_predictions,
        ) = parse_dbcan_output.parse_dbcan_output(
            test_set.prediction_paths["dbcan"],
            test_set.prediction_paths["hotpep"],
            test_set.fasta,
        )

        cupp_predictions = parse_cupp_output.parse_cupp_output(
            test_set.prediction_paths["cupp"],
            test_set.fasta,
        )

        ecami_predictions = parse_ecami_output.parse_ecami_output(
            test_set.prediction_paths["ecami"],
            test_set.fasta,
        )

        standardised_outputs = {
            "dbCAN": dbcan_predictions,
            "HMMER": hmmer_predictions,
            "Hotpep": hotpep_predictions,
            "DIAMOND": diamond_predictions,
            "CUPP": cupp_predictions,
            "eCAMI": ecami_predictions,
        }

        binary_classifications = []  # store all binary C/NC classifications for this test set

        for tool in tqdm(standardised_outputs, "Adding raw predictions to cumlative dfs"):

            # parse binary CAZyme/non-CAZyme predictions for the prediction tool
            # retrieves list of nested lists [accession, tool, prediction (0 or 1)]
            all_binary_classifications = binary_classification.get_cazyme_noncazyme_predictions(
                standardised_outputs[tool],
                tool,
                all_binary_classifications,
            )

            # parse multilabel CAZy family predictions, and build lists of nested lists of all
            # predicted and ground truth CAZy family annotations
            (
                family_predictions,
                family_ground_truths,
            ) = family_classifications.get_family_classifications(
                tool[0],
                tool[1],
                cazy_dict,
                test_set.source,
            )

            all_family_predictions += family_predictions
            all_family_ground_truths += family_ground_truths

        # build a datafame containing all binary C/NC classifications from all prediction
        # tools for this test set
        classifications_df = binary_classification.build_classification_df(
            all_binary_classifications,
        )

        # add ground truth (CAZy) binary CAZyme/non-CAZyme classifications to classifications.df
        classifications_df = binary_classification.add_ground_truths(
            classifications_df, cazy_dict,
        )

        # write out binary CAZyme/non-CAZyme predictions and ground truth annotations for test set
        # for documentation
        classifications_df.to_csv(f"binary_classifications_{time_stamp}_{test_set.source}.csv")

        all_binary_c_nc_dfs.append(ClassificationDF(test_set.source, classifications_df))
        # all_binary_c_nc_dfs used for bootstrapping accuracy of binary predictions

        # calculate statistics parameters for each tool performance, evaluating the binary
        # classification of CAZymes and non-CAZymes for the current working test set
        binary_c_nc_statistics += binary_classification.evaluate_binary_cazyme_noncazyme_predictions(
            classifications_df,
            test_set.source,
            args,
        )

    # build dataframes of predicted and ground truth CAZy family annotations
    # define the columns names for the CAZy family known and predicted annotations dataframes
    fam_dict = family_classifications.foundation_dict() # for cataloguing family predictions
    column_names = ["Genomic_accession", "Protein_accession", "Prediction_tool"]
    column_names += list(fam_dict.keys())

    all_family_predictions = pd.DataFrame(all_family_predictions, columns=column_names)
    all_family_ground_truths = pd.DataFrame(all_family_ground_truths, columns=column_names)

    # from the CAZy family prediction data, build dataframes of the CAZy class annotations
    (
        all_class_predictions,
        all_class_ground_truths,
    ) = class_classifications.build_class_annotation_dataframes(
        all_family_predictions,
        all_family_ground_truths,
        args,
        time_stamp,
    )

    # build dataframe of statisticl evaluations of binary CAZyme/non-CAZyme classifications
    binary_c_nc_statistics = pd.DataFrame(
        binary_c_nc_statistics,
        columns=[
            'Statistic parameter',
            'Genomic assembly',
            'Prediction tool',
        ],
    )

    return(
        all_binary_c_nc_dfs,
        binary_c_nc_statistics,
        all_class_predictions,
        all_class_ground_truths,
        all_family_predictions,
        all_family_ground_truths,
    )
