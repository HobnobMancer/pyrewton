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
"""This script is for evaluating the prediction of CAZy class annotations."""

import logging

import numpy as np
import pandas as pd

from tqdm import tqdm

from sklearn.metrics import fbeta_score, confusion_matrix, recall_score, precision_score
from sklearn.metrics.cluster import rand_score, adjusted_rand_score

from pyrewton.utilities import build_logger


def build_class_annotation_dataframes(fam_predictions_df, fam_ground_truths_df, args, time_stamp):
    """Build dataframe of CAZy class annotation predictions and another for ground truths.

    :param fam_predictions_df: df of predicted CAZy family annotations from all prediction tools
    :param fam_ground_truths_df: df of CAZy annotations of proteins
    :param args: cmd-line args parser
    :param time_stramp: str, time evaluation was started

    Return nothing.
    """
    class_ground_truths = []  # [[genomic_accession, protein_accession, prediction_tool, cazy_class]]
    class_predictions = []  # [[genomic_accession, protein_accession, prediction_tool, cazy_class]]

    dbcan_tool_predictions = {}  # {protein_accession: {HMMER:set(), DIAMOND:set(), Hotpep:set(): genomic_accession:str}}

    row_index = 0
    for row_index in tqdm(
        range(len(fam_ground_truths_df["Prediction_tool"])), desc="Getting CAZy class predictions",
    ):
        if row_predictions["Prediction_tool"] == "dbCAN":
            continue  # otherwise would ignore when different families are predicted for the same class
            # HMMER GH1, DIAMOND GH2 and Hotpep GH3 would result in no dbCAN consensus class annotation
            # so parser separately to find the consensus class annotation, irrespective the family annotation

        row_ground_truths = fam_ground_truths_df.iloc[row_index]
        row_predictions = fam_predictions_df.iloc[row_index]

        ground_truth_row_data = [
            row_ground_truths["Genomic_accession"],
            row_ground_truths["Protein_accession"],
            row_ground_truths["Prediction_tool"],
        ]
        pred_row_data = [
            row_predictions["Genomic_accession"],
            row_predictions["Protein_accession"],
            row_predictions["Prediction_tool"],
        ]

        if row_predictions["Prediction_tool"] in ['HMMER', 'DIAMOND', 'Hotpep']:
            try:
                dbcan_tool_predictions[row_predictions["Protein_accession"]]
            
            except KeyError:
                dbcan_tool_predictions[row_predictions["Protein_accession"]] = {
                    'genomic_accession': row_predictions["Genomic_accession"]
                }

        # get the prediction and ground truth for each CAZy class
        for cazy_class in ["GH", "GT", "PL", "CE", "AA", "CBM"]:
            # retrieve names of CAZy families in the CAZy class
            fam_names = [col_name for col_name in fam_ground_truths_df.columns if col_name.startswith(cazy_class)]

            class_fam_ground_truths = row_ground_truths[fam_names]
            class_fam_predictions = row_predictions[fam_names]

            # were any of the child families of the parent class predicted for the protein?
            if (class_fam_ground_truths == 1).any():
                ground_truth_class_classification = 1
            else:
                ground_truth_class_classification = 0

            ground_truth_row_data.append(ground_truth_class_classification)

            if (class_fam_predictions == 1).any():
                predicted_class_classification = 1
            else:
                predicted_class_classification = 0

            pred_row_data.append(predicted_class_classification)

            if row_predictions["Prediction_tool"] in ['HMMER', 'DIAMOND', 'Hotpep']:
                try:
                    dbcan_tool_predictions[row_predictions["Protein_accession"]][row_predictions["Prediction_tool"]]
                    dbcan_tool_predictions[row_predictions["Protein_accession"]][row_predictions["Prediction_tool"]][cazy_class] = predicted_class_classification
                except KeyError:
                    dbcan_tool_predictions[row_predictions["Protein_accession"]][row_predictions["Prediction_tool"]] = {cazy_class: predicted_class_classification}
                
                try:
                    dbcan_tool_predictions[row_predictions["Protein_accession"]]['CAZy']
                    dbcan_tool_predictions[row_predictions["Protein_accession"]]['CAZy'][cazy_class] = ground_truth_class_classification
                except KeyError:
                    dbcan_tool_predictions[row_predictions["Protein_accession"]]['CAZy'] = {cazy_class: ground_truth_class_classification}

        class_ground_truths.append(ground_truth_row_data)
        class_predictions.append(pred_row_data)

    protein_accessions = list(dbcan_tool_predictions.keys())
    for protein_accession in tqdm(protein_accessions, desc="Getting consensus dbCAN Class predictions"):
        ground_truth_row_data = [
            dbcan_tool_predictions[protein_accession]["genomic_accession"],
            protein_accession,
            'dbCAN',
        ]
        pred_row_data = [
            dbcan_tool_predictions[protein_accession]["genomic_accession"],
            protein_accession,
            'dbCAN',
        ]
        for cazy_class in ["GH", "GT", "PL", "CE", "AA", "CBM"]:
            # retrieve the CAZy class annotations from CAZy and the classifiers within dbCAN
            hmmer = dbcan_tool_predictions[protein_accession]['HMMER'][cazy_class]
            hotpep = dbcan_tool_predictions[protein_accession]['Hotpep'][cazy_class]
            diamond = dbcan_tool_predictions[protein_accession]['DIAMOND'][cazy_class]
            cazy = dbcan_tool_predictions[protein_accession]['CAZy'][cazy_class]

            # get the consensus CAZy class annotations from HMMER, DIAMOND and Hotpep
            # annotations are stored as int, therefore the sum is the numnber of predictions for the class
            total = hmmer + hotpep + diamond
            if total >= 2:
                dbcan = 1
            else:
                dbcan = 0
                
            pred_row_data.append(dbcan)
            ground_truth_row_data.append(cazy)

        class_ground_truths.append(ground_truth_row_data)
        class_predictions.append(pred_row_data)

    # build dataframes
    class_ground_truths_df = pd.DataFrame(
        class_ground_truths,
        columns=[
            "Genomic_accession",
            "Protein_accession",
            "Prediction_tool",
            "GH",
            "GT",
            "PL",
            "CE",
            "AA",
            "CBM",
        ],
    )
    class_predictions_df = pd.DataFrame(
        class_predictions,
        columns=[
            "Genomic_accession",
            "Protein_accession",
            "Prediction_tool",
            "GH",
            "GT",
            "PL",
            "CE",
            "AA",
            "CBM",
        ],
    )

    return class_predictions_df, class_ground_truths_df


def calculate_class_ari_ri(
    ground_truth_df,
    prediction_df,
    time_stamp,
    class_list=["GH", "GT", "PL", "CE", "AA", "CBM"],
):
    """Calculate the Adjusted Rand Index (ARI) and Rand Index (RI) of CAZy class predictions.

    :param predictions_df: df of predicted CAZy family annotations from all prediction tools
    :param ground_truths_df: df of CAZy annotations of proteins
    :param time_stamp: str
    :param class_list: list of CAZy class names

    Return predictions_df with the calc ARI and RI added in as new columns (ARI and RI calc
    per protein, and each protein is a separate row).
    """
    ri_scores = []
    ari_scores = []

    row_index = 0
    for row_index in tqdm(
        range(len(prediction_df["Prediction_tool"])), desc="Calculating CAZy class ARI and RI",
    ):
        ground_truths_row = ground_truth_df.iloc[row_index]
        predictions_row = prediction_df.iloc[row_index]

        y_true = ground_truths_row[class_list]
        y_pred = predictions_row[class_list]

        ri = rand_score(y_true, y_pred)
        ri_scores.append(ri)

        ari = adjusted_rand_score(y_true, y_pred)
        ari_scores.append(ari)

    prediction_df["Rand_index"] = ri_scores
    prediction_df["Adjusted_Rand_index"] = ari_scores

    return prediction_df


def calculate_class_stats(
    ground_truth_df,
    prediction_df,
    args,
    class_list=["GH", "GT", "PL", "CE", "AA", "CBM"],
):
    """Calculate the adjusted rand index and rand index of CAZy class annotation predictions.

    :param predictions_df: df of predicted CAZy family annotations from all prediction tools
    :param ground_truths_df: df of CAZy annotations of proteins
    :param time_stamp: str
    :param args: cmd-line parser
    :param class_list: list of CAZy class names

    Return predictions_df including calculated CAZy class statistics.
    """
    logger = logging.getLogger(__name__)
    specific_logger = build_logger(args.output, "cazy_class_confusion_matrix_errors.log")

    across_all_test_sets_data = []

    # evaluate performance across all test sets and all proteins
    for tool in ["dbCAN", "HMMER", "Hotpep", "DIAMOND", "CUPP", "eCAMI"]:
        # retrieve the rows of interest from from predictions and ground truths df
        tool_class_ground_truths = ground_truth_df.loc[ground_truth_df["Prediction_tool"] == tool]
        tool_class_predictions = prediction_df.loc[prediction_df["Prediction_tool"] == tool]

        # build an empty dataframe to store all predictions EXCEPT true negative non-CAZymes
        tp_fp_fn_ground_truths = pd.DataFrame(columns=list(tool_class_ground_truths.columns))
        tp_fp_fn_predictions = pd.DataFrame(columns=list(tool_class_predictions.columns))

        index = 0
        for index in range(len(tool_class_ground_truths["Prediction_tool"])):
            y_true = tool_class_ground_truths.iloc[index]
            y_true = list(y_true[class_list])  # retrieve only the cazy class 0/1 annotations
            
            y_pred = tool_class_predictions.iloc[index]
            y_pred = list(y_pred[class_list])  # retrieve only the cazy class 0/1 annotations

            if (1 not in y_true) and (1 not in y_pred):
                # if y_true and y_pred are all 0s, this is a true negative non-CAZyme prediction
                # do not include true negative non-CAZyme predictions
                continue
            else:
                # add TP, FP or FN result to the dataframes
                tp_fp_fn_ground_truths = tp_fp_fn_ground_truths.append(
                    tool_class_ground_truths.iloc[index],
                )
                tp_fp_fn_predictions = tp_fp_fn_predictions.append(
                    tool_class_predictions.iloc[index],
                )

        # calculate statistics
        for cazy_class in tqdm(class_list, desc=f"Calc CAZy class stats for {tool}"):
            data = [tool, cazy_class]

            y_true = list(tp_fp_fn_ground_truths[cazy_class])
            y_pred = list(tp_fp_fn_predictions[cazy_class])

            # check if the CAZy class was included in predictions and ground truths
            # if not, exclude the class from the evaluation
            if (1 not in y_true) and (1 not in y_pred):
                # do not include in statistics
                logger.warning(
                    f"{cazy_class} not predicted by {tool} and not in known "
                    f"annotations\nExcluding {cazy_class} from evaluation by setting all "
                    "stats results as NaN"
                )
                specific_logger.warning(
                    f"{cazy_class} not predicted by {tool} and not in known "
                    f"annotations\nExcluding {cazy_class} from evaluation by setting all "
                    "stats results as NaN"
                )

                sensitivity = np.nan
                data.append(sensitivity)
                specificity = np.nan
                data.append(specificity)
                precision = np.nan
                data.append(precision)
                fbeta = np.nan
                data.append(fbeta)
                accuracy = np.nan
                data.append(accuracy)

                across_all_test_sets_data.append(data)

                continue

            recall = recall_score(y_true, y_pred)
            data.append(recall)

            precision = precision_score(y_true, y_pred)
            data.append(precision)

            fbeta = fbeta_score(y_true, y_pred, beta=args.beta)
            data.append(fbeta)

            cm = confusion_matrix(y_true, y_pred)
            try:
                tn = cm[0][0]
                fn = cm[1][0]
                tp = cm[1][1]
                fp = cm[0][1]

            except IndexError:  # cannot build confusion matrix if true negative, or 
                data.append(np.nan)  # add specificity
                data.append(np.nan)  # add acuracy

                specific_logger.warning(
                    f"Prediction Tool: {tool}\tCAZy class: {cazy_class}\n"
                    f"y_true: {y_true}\ny_pred: {y_pred}\n"
                )

                continue

            specificity = tn / (tn + fp)
            data.append(specificity)

            accuracy = (tp + tn)/(tp + fp + fn + tn)
            data.append(accuracy)

            across_all_test_sets_data.append(data)

    # build dataframe using across_all_test_sets_data
    class_stats_df = pd.DataFrame(
        across_all_test_sets_data,
        columns=[
            "Prediction_tool",
            "CAZy_class",
            "Specificity",
            "Sensitivity",
            "Precision",
            "Fbeta_score",
            "Accuracy",
        ],
    )

    return class_stats_df


def calculate_class_stats_by_testsets(
    ground_truth_df,
    prediction_df,
    args,
    class_list=["GH", "GT", "PL", "CE", "AA", "CBM"],
):
    """Calculate statistical parameters for evaluating performance of CAZy clas prediction for each
    test set and each prediction tool.

    :param predictions_df: df of predicted CAZy family annotations from all prediction tools
    :param ground_truths_df: df of CAZy annotations of proteins
    :param time_stamp: str
    :param args: cmd-line args parser
    :param class_list: list of CAZy class names

    Return predictions_df including calculated CAZy class statistics.
    """
    specific_logger = build_logger(
        args.output,
        "cazy_class_confusion_matrix_errors_per_testset.log",
    )
    logger = logging.getLogger(__name__)

    # get the names of the genomic accessions, one accession = one test set
    all_genomic_accessions = ground_truth_df["Genomic_accession"]
    all_genomic_accessions = set(all_genomic_accessions)

    class_stats_data = []  # Genomic_accession, Prediction_tool, Statistic_parameter, Stat_value

    for tool in ["dbCAN", "HMMER", "Hotpep", "DIAMOND", "CUPP", "eCAMI"]:
        # retrieve the rows of interest from dataframes with the current working prediction tool
        tool_ground_truths = ground_truth_df.loc[ground_truth_df["Prediction_tool"] == tool]
        tool_predictions = prediction_df.loc[prediction_df["Prediction_tool"] == tool]

        # one accession is one test set
        for accession in tqdm(
            all_genomic_accessions,
            desc=f"Evaluate Class prediction by test set for {tool}:",
        ):
            testset_tool_ground_truths = tool_ground_truths.loc[
                tool_ground_truths["Genomic_accession"] == accession
            ]
            testset_tool_predictions = tool_predictions.loc[
                tool_predictions["Genomic_accession"] == accession
            ]

            # build empty dataframes to store all results EXCEPT true negative non-CAZymes
            tp_fp_fn_ground_truths = pd.DataFrame(columns=list(testset_tool_ground_truths.columns))
            tp_fp_fn_predictions = pd.DataFrame(columns=list(testset_tool_predictions.columns))
           
            # exclude true negative non-CAZyme predictions
            index = 0
            for index in range(len(testset_tool_predictions["Genomic_accession"])):
                y_true = testset_tool_ground_truths.iloc[index]
                y_true = list(y_true[class_list])  # retrieve only the cazy class 0/1 annotations

                y_pred = testset_tool_predictions.iloc[index]
                y_pred = list(y_pred[class_list])  # retrieve only the cazy class 0/1 annotations

                if (1 not in y_true) and (1 not in y_pred):
                    # if y_true and y_pred are all 0s, this is a true negative non-CAZyme prediction
                    # do not include true negative non-CAZyme predictions
                    continue

                else:
                    tp_fp_fn_ground_truths = tp_fp_fn_ground_truths.append(
                        testset_tool_ground_truths.iloc[index],
                    )
                    tp_fp_fn_predictions = tp_fp_fn_predictions.append(
                        testset_tool_predictions.iloc[index],
                    )

            for cazy_class in class_list:
                y_true = list(tp_fp_fn_ground_truths[cazy_class])
                y_pred = list(tp_fp_fn_predictions[cazy_class])

                # check if the CAZy class was included in predictions and ground truths
                # if not exclude the class from the evaluation
                if (1 not in y_true) and (1 not in y_pred):
                    # do not include in statistics
                    logger.warning(
                        f"{cazy_class} not predicted by {tool} in {accession} and not in known "
                        f"annotations\nExcluding {cazy_class} from evaluation by setting all "
                        "stats results as NaN"
                    )
                    specific_logger.warning(
                        f"{cazy_class} not predicted by {tool} in {accession} and not in known "
                        f"annotations\nExcluding {cazy_class} from evaluation by setting all "
                        "stats results as NaN"
                    )

                    specificity = np.nan
                    class_stats_data.append(
                        [accession, tool, cazy_class, "Specificity", specificity],
                    )

                    sensitivity = np.nan
                    class_stats_data.append([accession, tool, cazy_class, "Sensitivity", sensitivity])
                    
                    precision = np.nan
                    class_stats_data.append([accession, tool, cazy_class, "Precision", precision])
                    
                    fbeta = np.nan
                    class_stats_data.append([accession, tool, cazy_class, "Fbeta_score", fbeta])
                    
                    accuracy = np.nan
                    class_stats_data.append([accession, tool, cazy_class, "Accuracy", accuracy])

                    continue

                recall = recall_score(y_true, y_pred)
                class_stats_data.append([accession, tool, cazy_class, "Sensitivity", recall])

                precision = precision_score(y_true, y_pred)
                class_stats_data.append([accession, tool, cazy_class, "Precision", precision])

                fbeta = fbeta_score(y_true, y_pred, beta=args.beta)
                class_stats_data.append([accession, tool, cazy_class, "Fbeta_score", fbeta])

                cm = confusion_matrix(y_true, y_pred)
                try:
                    tn = cm[0][0]
                    fn = cm[1][0]
                    tp = cm[1][1]
                    fp = cm[0][1]

                    accuracy = (tp + tn)/(tp + fp + fn + tn)
                    class_stats_data.append([accession, tool, cazy_class, "Accuracy", accuracy])

                except IndexError as e:
                    class_stats_data.append([accession, tool, cazy_class, "Accuracy", np.nan])
                    logger.warning(
                        f"Error raised when creating confusion matrix for protein {accession}, "
                        f"{tool}, {cazy_class}, when calculating accuracy.\nError raised:\n{e}"
                    )
                    specific_logger.warning(
                        f"Prediction Tool: {tool}\tClass: {cazy_class}\tAccession: {accession}\t"
                        f"Stat: Specificity\ny_true: {y_true}\ny_pred: {y_pred}\n"
                    )

                try:
                    tn = cm[0][0]
                    fp = cm[0][1]
                    specificity = tn / (tn + fp)

                    class_stats_data.append(
                        [
                            accession,
                            tool,
                            cazy_class,
                            "Specificity",
                            specificity,
                        ],
                    )

                except IndexError as e:
                    class_stats_data.append([accession, tool, cazy_class, "Specificity", np.nan])
                    logger.warning(
                        f"Error raised when creating confusion matrix for protein {accession}, "
                        f"{tool}, {cazy_class}, when calculating specificity.\nError raised:\n{e}"
                    )
                    specific_logger.warning(
                        f"Prediction Tool: {tool}\tClass: {cazy_class}\tAccession: {accession}\t"
                        f"Stat: Specificity\ny_true: {y_true}\ny_pred: {y_pred}\n"
                    )

    # build dataframe using class_stats_data
    class_stats_df = pd.DataFrame(
        class_stats_data,
        columns = [
            "Genomic_accession",
            "Prediction_tool",
            "CAZy_class",
            "Statistic_parameter",
            "Statistic_value",
        ],
    )

    return class_stats_df


def evaluate_taxa_performance(
    all_ground_truths,
    all_predictions,
    tax_dict,
    time_stamp,
    args,
):
    """Evaluate performance per CAZy family per tax group
    
    :param all_predictions: Pandas df, 
        columns: Genomic_accession, Protein_accession, Prediction_tool, one column per CAZy family
            denoting if annotation predicted (1) or not predicted (0), Rand_index and Adjusted_rand_index
    :param all_ground_truths: Pandas df, 
        columns: Genomic_accession, Protein_accession, Prediction_tool, one column per CAZy family
            denoting if annotation given by CAZy (1) or not  (0)
    :param class_predictions_df: Pandas df,
        columns: Genomic_accession, Protein_accession, Prediction_tool, one column per CAZy class
            denoting if annotation predicted (1) or not predicted (0), Rand_index and Adjusted_rand_index
    :param tax_dict: dict of tax groups, keyed by tax group name, values by list of genomic accessions
    :param time_stamp: str, date and time evaluation was invoked
    :param args: cmd-line args parser
    
    Return nothing.
    """
    # evaluate performance per tax group per family
    for tax_group in tqdm(tax_dict, 'Evaluating performance per tax group per CAZy class'):
        genomic_accessions = tax_dict[tax_group]

        # create empty dfs to store rows of interest
        gt_col_names = list(all_ground_truths.columns)
        gt_col_names.append('Tax_group')
        tax_ground_truths = pd.DataFrame(columns=gt_col_names)

        pred_col_names = list(all_predictions.columns)
        pred_col_names.append('Tax_group')
        tax_predictions = pd.DataFrame(columns=pred_col_names)

        for accession in genomic_accessions:
            grnd_trth_df = tax_ground_truths.loc[tax_ground_truths["Genomic_accession"] == accession]
            pred_df = tax_predictions.loc[tax_predictions["Genomic_accession"] == accession]

            # add tax_group column
            gt_tax_group_col = [tax_group] * len(grnd_trth_df['Genomic_accession'])
            pred_tax_group_col = [tax_group] * len(pred_df['Genomic_accession'])

            grnd_trth_df['Tax_group'] = gt_tax_group_col
            pred_df['Tax_group'] = pred_tax_group_col

            tax_ground_truths = tax_ground_truths.append(grnd_trth_df)
            tax_predictions = tax_predictions.append(pred_df)

        class_stats_df = calculate_class_stats(
            tax_ground_truths,
            tax_predictions,
            args,
        )

        output_path = args.output / f"class_stats_all_test_sets_tax_comparison_{tax_group}_{time_stamp}.csv"
        class_stats_df.to_csv(output_path)

        class_stats_df_testset = calculate_class_stats_by_testsets(
            tax_ground_truths,
            tax_predictions,
            args,
        )

        output_path = args.output / f"class_stats_per_test_set_tax_comparison_{tax_group}_{time_stamp}.csv"
        class_stats_df_testset.to_csv(output_path)

    return

