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


import logging

import numpy as np
import pandas as pd

from tqdm import tqdm

from sklearn.metrics import fbeta_score, confusion_matrix, recall_score, precision_score
from sklearn.metrics.cluster import rand_score, adjusted_rand_score


def get_family_classifications(predictions, prediciton_tool, cazy_dict, genomic_accession):
    """Retrieve the predicted  and ground truth (CAZy determined) CAZy family annotations.
    
    :param predictions:
    :param prediction_tools:
    :param cazy_dict:
    :param genomic_accession:

    Return two lists, one of predicted results, one of ground truths.
    Each list is a list of nested lists. Each nested list contains:
    [genomic_accession, protein_accession, prediction_tool, one column for each fam with 
    its 0 (not predicted) / 1 (predicted) prediction].
    """
    logger = logging.getLogger(__name__)

    all_predicted_annotations = []  # [[assembly, protein_accession, prediciton_tool, fam1, fam2..]]
    all_ground_truths = []  # [[assembly, protein_accession, prediciton_tool, fam1, fam2, ...]]

    for protein in predictions:
        # retrieve dicts of {cazy_family_name: 0} for storing CAZy family annotation predictions
        fam_predictions = foundation_dict()
        known_fams = foundation_dict()

        new_predictions = [genomic_accession, protein.protein_accession, prediciton_tool]
        new_ground_truths = [genomic_accession, protein.protein_accession, prediciton_tool]

        for domain in protein.cazyme_domains:
            try:
                fam_predictions[domain.cazy_family]
                fam_predictions[domain.cazy_family] = 1  # add one to cazy family 
            except KeyError:
                # check if cazy family annotation for the domain is a subfamily. If it is, mark the
                # parent family in the fam_predictions dict as being predicted
                if type(domain.cazy_family) is str:
                    if (domain.cazy_family).find("_") != -1:
                        fam = domain.cazy_family[:(domain.cazy_family).find("_")]
                        try:
                            fam_predictions[fam]
                            fam_predictions[fam] = 1
                        except KeyError:
                            logger.warning(
                                f"Did not recognise {fam} as a CAZy family "
                                f"from prediction {domain}, from {prediciton_tool}"
                            )
                    else:
                        logger.warning(
                            f"Did not recognise {fam} as a CAZy family "
                            f"from prediction {domain}, from {prediciton_tool}"
                        )

                elif np.isnan(domain.cazy_family):
                    # cazy family set as null value for the CAZy domain
                    logger.warning(
                        f"CAZy family set as null value for {domain}, from {prediciton_tool}"
                    )
                    continue

                else:
                    logger.warning(
                        f"Did not recognise {fam} as a CAZy family "
                        f"from prediction {domain}, from {prediciton_tool}"
                    )

        # add cazy family (0/1) predictions to list representing CAZy family predictions
        # for the current working protein
        new_prediction += list(fam_predictions.values())

        # retrieve the ground thruth (CAZy determined) CAZy family annotations for the protein
        try:
            cazy_annotations = cazy_dict[protein.protein_accession]
            try:
                for fam in cazy_annotations:
                    fam = fam.strip()
                    known_fams[fam]
                    known_fams[fam] = 1
            except KeyError as e:
                logger.warning(
                    f"KeyError raised for known CAZy family {fam}. Error raised:\n{e}"
                )
        except KeyError:
            pass

        # add cazy family (0/1) ground truth annotations to list representing all 
        # CAZy family ground truth annotations for the current working protein
        new_ground_truths += list(known_fams.values())

        all_predicted_annotations.append(new_prediction)
        all_ground_truths.append(new_ground_truths)

    return all_predicted_annotations, all_ground_truths


def calc_fam_fbeta_score(predictions_df, ground_truths_df, time_stamp, args):
    """Calculate the Fbeta score per CAZy family, per prediction tool and write out to csv file.

    :param predictions_df: df of predicted CAZy family annotations from all prediction tools
    :param ground_truths_df: df of CAZy annotations of proteins
    :param time_stamp: str, time evaluation was started
    :param args: cmd-line args parser

    Return nothing.
    """
    for tool in ["dbCAN", "HMMER", "Hotpep", "DIAMOND", "CUPP", "eCAMI"]:
        # retrieve the rows of interest from from predictions and ground truths df
        pred_df = predictions_df.loc[predictions_df["Prediction_tool"] == tool]
        grnd_trth_df = ground_truths_df.loc[ground_truths_df["Prediction_tool"] == tool]

        # list will form the row of calc fbeta scores that is later added to the dataframe
        fbeta_scores = [np.nan, np.nan, np.nan, np.nan]  # Nas reflect non-CAZy family cols in df

        # iterate through the columns
        # each column contains a unique CAZy family, except the first few columns
        for fam in tqdm(pred_df.columns[4:], desc=f"Calc fam Fbeta scores for {tool}"):

            # check if true negative, these are not taken into consideration by the Fbeta score
            pred = (pred_df[fam] == 0).all()
            known = (grnd_trth_df[fam] == 0).all()
            if pred and known:
                fbeta = np.nan
                fbeta_scores.append(fbeta)
                continue

            fbeta = fbeta_score(grnd_trth_df[fam], pred_df[fam], beta=args.beta)
            fbeta_scores.append(fbeta)

        # build df out of fbeta scores and add to the predicted df
        df = pd.DataFrame([fbeta_scores], columns=pred_df.columns)
        pred_df = pred_df.append(df)

        pred_df.to_csv(f"{tool}_fam_prediction_df_{time_stamp}.csv")

    return


def calc_fam_fbeta_score_per_testset(
    predictions_df,
    ground_truths_df,
    time_stamp,
    predictions,
    args,
):
    """Calculate the Fbeta score per CAZy family, per prediction tool, per test set and write out
    to csv file.

    :param predictions_df: df of predicted CAZy family annotations from all prediction tools
    :param ground_truths_df: df of CAZy annotations of proteins
    :param time_stramp: str, time statistical evaluation was initated, used for naming files
    :param predictions: list of Prediction instances, used for retrieving genomic accession of each
    test set
    :param args: cmd-line args parser

    Return nothing.
    """
    df_data = []  # [[CAZy_family, Prediction_tool, Genomic_accession, Fbeta_score]]

    for testset in tqdm(predictions, desc="Calc CAZy Fam Fbeta-score per test set"):
        # retrieve all rows containing proteins from the current working test set
        # test sets identified by genomic accession
        genomic_accession = testset.source
        pred_df = predictions_df.loc[predictions_df["Genomic_accession"] == genomic_accession]
        grnd_trth_df = ground_truths_df.loc[ground_truths_df["Genomic_accession"] == genomic_accession]

        for tool in ["dbCAN", "HMMER", "Hotpep", "DIAMOND", "CUPP", "eCAMI"]:
            # retrieve the rows of interest from from predictions and ground truths df
            pred_df = pred_df.loc[pred_df["Prediction_tool"] == tool]
            grnd_trth_df = grnd_trth_df.loc[grnd_trth_df["Prediction_tool"] == tool]

            # iterate through the columns, each column contains a unique CAZy family,
            # except the first few columns
            for fam in tqdm(pred_df.columns[4:], desc=f"Calc fam Fbeta scores for {tool}"):

                # check if true negative, these are not taken into consideration by the Fbeta score
                pred = (pred_df[fam] == 0).all()
                known = (grnd_trth_df[fam] == 0).all()
                if pred and known:
                    fbeta = np.nan
                    results = [fam, tool, genomic_accession, fbeta]
                    df_data.append(results)
                    continue

                # calculate the Fbeta score for the CAZy fam, for current working prediction tool,
                # for the current working test set (identified by its source genomic accession)
                fbeta = fbeta_score(grnd_trth_df[fam], pred_df[fam], beta=args.beta)

                results = [fam, tool, genomic_accession, fbeta]
                df_data.append(results)

    # build df out of fbeta scores and add to the predicted df
    column_names = ["CAZy_family", "Prediction_tool", "Genomic_accession", "Fbeta_score"]
    df = pd.DataFrame(df_data, columns=column_names)
    df.to_csv(f"cazy_fam_fbeta_scores_{time_stamp}.csv")

    return


def calc_fam_stats(predictions_df, ground_truths_df, time_stamp, args):
    """Calculate the Fbeta score, accuracy, recall and specificity per CAZy family, per prediction
    tool and write out to csv file.

    :param predictions_df: df of predicted CAZy family annotations from all prediction tools
    :param ground_truths_df: df of CAZy annotations of proteins
    :param time_stamp: str, time evaluation was started
    :param args: cmd_line args parser

    Return nothing.
    """
    logger = logging.getLogger(__name__)
    stats_orientated_df = []  # [[CAZyFam, Prediction_Tool, Specificity, Recall, Fbeta, Accuracy]]
    long_dataframe = []  # [[CAZyFam, StatParameter, PredicitonTool, StatValue]]

    for tool in ["dbCAN", "HMMER", "Hotpep", "DIAMOND", "CUPP", "eCAMI"]:
        # retrieve the rows of interest from from predictions and ground truths df
        pred_df = predictions_df.loc[predictions_df["Prediction_tool"] == tool]
        grnd_trth_df = ground_truths_df.loc[ground_truths_df["Prediction_tool"] == tool]

        # iterate through the columns, each column contains a unique CAZy family,
        # except the first few columns
        for fam in tqdm(pred_df.columns[4:], desc=f"Calc stats for fams for {tool}"):

            y_true = grnd_trth_df[fam]
            y_pred = pred_df[fam]

            # check if true negative, these are not taken into consideration by the Fbeta score
            pred = (pred_df[fam] == 0).all()
            known = (grnd_trth_df[fam] == 0).all()
            if pred and known:
                stats_orientated_df.append(
                    [
                        fam,
                        tool,
                        np.nan,
                        np.nan,
                        np.nan,
                        np.nan,
                        np.nan,
                        np.nan,
                        np.nan,
                        np.nan,
                        np.nan,
                    ],
                )
                long_dataframe.append([fam, "specificity", tool, np.nan])
                long_dataframe.append([fam, "recall", tool, np.nan])
                long_dataframe.append([fam, "fbeta", tool, np.nan])
                long_dataframe.append([fam, "accuracy", tool, np.nan])
                long_dataframe.append([fam, "TNs", tool, np.nan])
                long_dataframe.append([fam, "FNs", tool, np.nan])
                long_dataframe.append([fam, "TPs", tool, np.nan])
                long_dataframe.append([fam, "FPs", tool, np.nan])
                long_dataframe.append([fam, "recall_sample_size", tool, np.nan])
                continue

            cm = confusion_matrix(y_true, y_pred)
            try:
                tn = cm[0][0]
                fn = cm[1][0]
                tp = cm[1][1]
                fp = cm[0][1]

                specificity = tn / (tn + fp)
                long_dataframe.append([fam, "specificity", tool, specificity])

                recall = recall_score(y_true, y_pred)
                long_dataframe.append([fam, "recall", tool, recall])

                fbeta = fbeta_score(grnd_trth_df[fam], pred_df[fam], beta=args.beta)
                long_dataframe.append([fam, "fbeta", tool, fbeta])

                accuracy = (tp + tn)/(tp + fp + fn + tn)
                long_dataframe.append([fam, "accuracy", tool, accuracy])

                # find the sample size
                long_dataframe.append([fam, "TNs", tool, tn])
                long_dataframe.append([fam, "FNs", tool, fn])
                long_dataframe.append([fam, "TPs", tool, tp])
                long_dataframe.append([fam, "FPs", tool, fp])
                long_dataframe.append([fam, "recall_sample_size", tool, (tp+fn)])
                stats_orientated_df.append(
                    [
                        fam,
                        tool,
                        specificity,
                        recall,
                        fbeta,
                        accuracy,
                        tn,
                        fn,
                        tp,
                        fp,
                        (tp+fn),
                    ],
                )

            except IndexError:  # no idea why
                stats_orientated_df.append(
                    [
                        fam,
                        tool,
                        np.nan,
                        np.nan,
                        np.nan,
                        np.nan,
                        np.nan,
                        np.nan,
                        np.nan,
                        np.nan,
                        np.nan,
                    ],
                )
                long_dataframe.append([fam, "specificity", tool, np.nan])
                long_dataframe.append([fam, "recall", tool, np.nan])
                long_dataframe.append([fam, "fbeta", tool, np.nan])
                long_dataframe.append([fam, "accuracy", tool, np.nan])
                long_dataframe.append([fam, "TNs", tool, np.nan])
                long_dataframe.append([fam, "FNs", tool, np.nan])
                long_dataframe.append([fam, "TPs", tool, np.nan])
                long_dataframe.append([fam, "FPs", tool, np.nan])
                long_dataframe.append([fam, "recall_sample_size", tool, np.nan])
                logger.warning(
                    "Error raised when creating confusion matrix for:\n"
                    f"{tool}\t{fam}"
                )

    # build df out of fbeta scores and add to the predicted df
    df = pd.DataFrame(
        long_dataframe,
        columns=["CAZy_family", "Stat_parameter", "Prediction_tool", "Stat_value"],
    )
    df.to_csv(f"fam_stats_df_{time_stamp}.csv")

    df1 = pd.DataFrame(
        stats_orientated_df,
        columns=[
            "CAZy_family",
            "Prediction_tool",
            "Specificity",
            "Recall",
            "Fbeta_score",
            "Accuracy",
            "#TN",
            "#FN",
            "#TP",
            "#FP",
            "Recall_sample_size",
        ],
    )
    df1.to_csv(f"fam_stats_df_single_{time_stamp}.csv")

    return


def build_class_dataframes(predictions_df, ground_truths_df, args, time_stamp):
    """Build dataframe of CAZy class annotation predictions.

    :param predictions_df: df of predicted CAZy family annotations from all prediction tools
    :param ground_truths_df: df of CAZy annotations of proteins
    :param args: cmd-line args parser
    :param time_stramp: str, time evaluation was started

    Return nothing.
    """
    class_ground_truths = []  # genomic_accession, protein_accession, prediction_tool, cazy_class
    class_predictions = []  # genomic_accession, protein_accession, prediction_tool, cazy_class

    row_index = 0
    for row_index in tqdm(range(len(ground_truths_df["Prediction_tool"])), desc="Getting CAZy class predictions"):
        row_ground_truths = ground_truths_df.iloc[row_index]
        row_predictions = predictions_df.iloc[row_index]

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

        # get the prediction and ground truth for each CAZy class
        for cazy_class in ["GH", "GT", "PL", "CE", "AA", "CBM"]:
            # retrieve names of CAZy families in the CAZy class
            fam_names = [col_name for col_name in ground_truths_df.columns if col_name.startswith(cazy_class)]

            class_fam_ground_truths = row_ground_truths[fam_names]
            class_fam_predictions = row_predictions[fam_names]

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

    # write out dataframe
    class_ground_truths_df.to_csv(f"class_ground_truths_classifications_{time_stamp}.csv")

    # calculate ARI and RI for multilabel evaluation
    class_predictions_df = calculate_class_ari(
        class_ground_truths_df,
        class_predictions_df,
        time_stamp,
    )
    class_predictions_df.to_csv(f"class_predictions_classifications_{time_stamp}.csv")

    # calculate the fbeta_score, sensitivity(recall), specificity and accuracy per CAZy class
    calculate_class_stats(class_ground_truths_df, class_predictions_df, time_stamp, args)
    calculate_class_stats_by_testsets(
        class_ground_truths_df,
        class_predictions_df,
        time_stamp,
        args,
    )

    return


def calculate_class_ari(
    ground_truth_df,
    prediction_df,
    time_stamp,
    class_list=["GH", "GT", "PL", "CE", "AA", "CBM"],
):
    """Calculate the adjusted rand index and rand index of CAZy class annotation predictions.

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
    for row_index in tqdm(range(len(prediction_df["Prediction_tool"])), desc="Calculating CAZy class ARI and RI"):
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
    time_stamp,
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
    across_all_test_sets_data = []

    # evaluate performance across all test sets and all proteins
    for tool in ["dbCAN", "HMMER", "Hotpep", "DIAMOND", "CUPP", "eCAMI"]:
        # retrieve the rows of interest from from predictions and ground truths df
        tool_class_ground_truths = ground_truth_df.loc[ground_truth_df["Prediction_tool"] == tool]
        tool_class_predictions = prediction_df.loc[prediction_df["Prediction_tool"] == tool]

        for cazy_class in tqdm(class_list, desc=f"Calc CAZy class stats for {tool}"):
            data = [tool, cazy_class]

            y_true = tool_class_ground_truths[cazy_class]
            y_pred = tool_class_predictions[cazy_class]

            cm = confusion_matrix(y_true, y_pred)
            try:
                tn = cm[0][0]
                fn = cm[1][0]
                tp = cm[1][1]
                fp = cm[0][1]

                specificity = tn / (tn + fp)

                recall = recall_score(y_true, y_pred)

                precision = precision_score(y_true, y_pred)

                fbeta = fbeta_score(y_true, y_pred, beta=args.beta)

                accuracy = (tp + tn)/(tp + fp + fn + tn)

                data.append(specificity)
                data.append(recall)
                data.append(precision)
                data.append(fbeta)
                data.append(accuracy)

            except IndexError:
                data.append(np.nan)
                data.append(np.nan)
                data.append(np.nan)
                data.append(np.nan)
                data.append(np.nan)

            across_all_test_sets_data.append(data)

    # build dataframe using across_all_test_sets_data
    class_stats_df = pd.DataFrame(
        across_all_test_sets_data,
        columns=[
            "Prediction_tool",
            "CAZy_class",
            "Specificity",
            "Recall",
            "Precision",
            "Fbeta_score",
            "Accuracy",
        ],
    )
    class_stats_df.to_csv(f"class_stats_across_all_test_sets_{time_stamp}.csv")
    return


def calculate_class_stats_by_testsets(
    ground_truth_df,
    prediction_df,
    time_stamp,
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
    logger = logging.getLogger(__name__)

    # get the names of the genomic accessions, one accession = one test set
    all_genomic_accessions = ground_truth_df["Genomic_accession"]
    all_genomic_accessions = set(all_genomic_accessions)

    class_stats_data = []  # Genomic_accession, Prediction_tool, Statistic_parameter, Stat_value

    for tool in ["dbCAN", "HMMER", "Hotpep", "DIAMOND", "CUPP", "eCAMI"]:
        # retrieve the rows of interest from from predictions and ground truths df
        tool_class_ground_truths = ground_truth_df.loc[ground_truth_df["Prediction_tool"] == tool]
        tool_class_predictions = prediction_df.loc[prediction_df["Prediction_tool"] == tool]

        for accession in tqdm(all_genomic_accessions, desc="Evaluate Class prediction by test set"):
            test_set_tool_class_ground_truths = tool_class_ground_truths.loc[tool_class_ground_truths["Genomic_accession"] == accession]
            test_set_tool_class_predictions = tool_class_predictions.loc[tool_class_predictions["Genomic_accession"] == accession]

            for cazy_class in class_list:
                y_true = test_set_tool_class_ground_truths[cazy_class]
                y_pred = test_set_tool_class_predictions[cazy_class]

                # exclude true negatives
                # (where protein is predicted to be a non-CAZyme and the protein
                # is not included in CAZy)
                pred = (y_true == 0).all()
                known = (y_true == 0).all()
                if pred and known:
                    class_stats_data.append([accession, tool, cazy_class, "Specificity", np.nan])
                    class_stats_data.append([accession, tool, cazy_class, "Recall", np.nan])
                    class_stats_data.append([accession, tool, cazy_class, "Precision", np.nan])
                    class_stats_data.append([accession, tool, cazy_class, "Fbeta_score", np.nan])
                    class_stats_data.append([accession, tool, cazy_class, "Accuracy", np.nan])
                    continue

                recall = recall_score(y_true, y_pred)
                class_stats_data.append([accession, tool, cazy_class, "Recall", recall])

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
    class_stats_df.to_csv(f"class_stats_per_test_set_{time_stamp}.csv")

    return


def foundation_dict():
    """Dict use as foundation for CAZy family predictions."""
    foundation_dict = {
        'GH1': 0, 'GH2': 0, 'GH3': 0, 'GH4': 0, 'GH5': 0, 'GH6': 0, 'GH7': 0,
        'GH8': 0, 'GH9': 0, 'GH10': 0, 'GH11': 0, 'GH12': 0, 'GH13': 0, 'GH14': 0, 'GH15': 0,
        'GH16': 0, 'GH17': 0, 'GH18': 0, 'GH19': 0, 'GH20': 0, 'GH21': 0, 'GH22': 0, 'GH23': 0,
        'GH24': 0, 'GH25': 0, 'GH26': 0, 'GH27': 0, 'GH28': 0, 'GH29': 0, 'GH30': 0, 'GH31': 0,
        'GH32': 0, 'GH33': 0, 'GH34': 0, 'GH35': 0, 'GH36': 0, 'GH37': 0, 'GH38': 0, 'GH39': 0,
        'GH40': 0, 'GH41': 0, 'GH42': 0, 'GH43': 0, 'GH44': 0, 'GH45': 0, 'GH46': 0, 'GH47': 0,
        'GH48': 0, 'GH49': 0, 'GH50': 0, 'GH51': 0, 'GH52': 0, 'GH53': 0, 'GH54': 0, 'GH55': 0,
        'GH56': 0, 'GH57': 0, 'GH58': 0, 'GH59': 0, 'GH60': 0, 'GH61': 0, 'GH62': 0, 'GH63': 0,
        'GH64': 0, 'GH65': 0, 'GH66': 0, 'GH67': 0, 'GH68': 0, 'GH69': 0, 'GH70': 0, 'GH71': 0,
        'GH72': 0, 'GH73': 0, 'GH74': 0, 'GH75': 0, 'GH76': 0, 'GH77': 0, 'GH78': 0, 'GH79': 0,
        'GH80': 0, 'GH81': 0, 'GH82': 0, 'GH83': 0, 'GH84': 0, 'GH85': 0, 'GH86': 0, 'GH87': 0,
        'GH88': 0, 'GH89': 0, 'GH90': 0, 'GH91': 0, 'GH92': 0, 'GH93': 0, 'GH94': 0, 'GH95': 0,
        'GH96': 0, 'GH97': 0, 'GH98': 0, 'GH99': 0, 'GH100': 0, 'GH101': 0, 'GH102': 0, 'GH103': 0,
        'GH104': 0, 'GH105': 0, 'GH106': 0, 'GH107': 0, 'GH108': 0, 'GH109': 0, 'GH110': 0,'GH111': 0,
        'GH112': 0, 'GH113': 0, 'GH114': 0, 'GH115': 0, 'GH116': 0, 'GH117': 0, 'GH118': 0, 'GH119': 0,
        'GH120': 0, 'GH121': 0, 'GH122': 0, 'GH123': 0, 'GH124': 0, 'GH125': 0, 'GH126': 0, 'GH127': 0,
        'GH128': 0, 'GH129': 0, 'GH130': 0, 'GH131': 0, 'GH132': 0, 'GH133': 0, 'GH134': 0, 'GH135': 0, 'GH136': 0, 'GH137': 0, 'GH138': 0, 'GH139': 0,
        'GH140': 0, 'GH141': 0, 'GH142': 0, 'GH143': 0, 'GH144': 0, 'GH145': 0, 'GH146': 0, 'GH147': 0, 'GH148': 0, 'GH149': 0, 'GH150': 0,
        'GH151': 0, 'GH152': 0, 'GH153': 0, 'GH154': 0, 'GH155': 0, 'GH156': 0, 'GH157': 0, 'GH158': 0, 'GH159': 0, 'GH160': 0, 'GH161': 0,
        'GH162': 0, 'GH163': 0, 'GH164': 0, 'GH165': 0, 'GH166': 0, 'GH167': 0, 'GH168': 0, 'GH169': 0, 'GH170': 0, 'GH171': 0, 'GH0': 0,
        'GT1': 0, 'GT2': 0, 'GT3': 0, 'GT4': 0, 'GT5': 0, 'GT6': 0, 'GT7': 0, 'GT8': 0, 'GT9': 0, 'GT10': 0, 'GT11': 0, 'GT12': 0, 'GT13': 0,
        'GT14': 0, 'GT15': 0, 'GT16': 0, 'GT17': 0, 'GT18': 0, 'GT19': 0, 'GT20': 0, 'GT21': 0, 'GT22': 0, 'GT23': 0, 'GT24': 0, 'GT25': 0,
        'GT26': 0, 'GT27': 0, 'GT28': 0, 'GT29': 0, 'GT30': 0, 'GT31': 0, 'GT32': 0, 'GT33': 0, 'GT34': 0, 'GT35': 0, 'GT36': 0, 'GT37': 0,
        'GT38': 0, 'GT39': 0, 'GT40': 0, 'GT41': 0, 'GT42': 0, 'GT43': 0, 'GT44': 0, 'GT45': 0, 'GT46': 0, 'GT47': 0, 'GT48': 0, 'GT49': 0,
        'GT50': 0, 'GT51': 0, 'GT52': 0, 'GT53': 0, 'GT54': 0, 'GT55': 0, 'GT56': 0, 'GT57': 0, 'GT58': 0, 'GT59': 0, 'GT60': 0, 'GT61': 0,
        'GT62': 0, 'GT63': 0, 'GT64': 0, 'GT65': 0, 'GT66': 0, 'GT67': 0, 'GT68': 0, 'GT69': 0, 'GT70': 0, 'GT71': 0, 'GT72': 0, 'GT73': 0,
        'GT74': 0, 'GT75': 0, 'GT76': 0, 'GT77': 0, 'GT78': 0, 'GT79': 0, 'GT80': 0, 'GT81': 0, 'GT82': 0, 'GT83': 0, 'GT84': 0, 'GT85': 0,
        'GT86': 0, 'GT87': 0, 'GT88': 0, 'GT89': 0, 'GT90': 0, 'GT91': 0, 'GT92': 0, 'GT93': 0, 'GT94': 0, 'GT95': 0, 'GT96': 0, 'GT97': 0,
        'GT98': 0, 'GT99': 0, 'GT100': 0, 'GT101': 0, 'GT102': 0, 'GT103': 0, 'GT104': 0, 'GT105': 0, 'GT106': 0, 'GT107': 0, 'GT108': 0,
        'GT109': 0, 'GT110': 0, 'GT111': 0, 'GT112': 0, 'GT113': 0, 'GT114': 0, 'GT0': 0,
        'PL1': 0, 'PL2': 0, 'PL3': 0, 'PL4': 0, 'PL5': 0, 'PL6': 0, 'PL7': 0, 'PL8': 0, 'PL9': 0, 'PL10': 0, 'PL11': 0, 'PL12': 0, 'PL13': 0,
        'PL14': 0, 'PL15': 0,  'PL16': 0, 'PL17': 0, 'PL18': 0, 'PL19': 0, 'PL20': 0, 'PL21': 0, 'PL22': 0, 'PL23': 0, 'PL24': 0, 'PL25': 0,
        'PL26': 0, 'PL27': 0, 'PL28': 0, 'PL29': 0, 'PL30': 0, 'PL31': 0, 'PL32': 0, 'PL33': 0, 'PL34': 0, 'PL35': 0, 'PL36': 0, 'PL37': 0,
        'PL38': 0, 'PL39': 0, 'PL40': 0, 'PL41': 0, 'PL0': 0,
        'CE1': 0, 'CE2': 0, 'CE3': 0, 'CE4': 0, 'CE5': 0, 'CE6': 0, 'CE7': 0, 'CE8': 0, 'CE9': 0, 'CE10': 0, 'CE11': 0, 'CE12': 0, 'CE13': 0,
        'CE14': 0, 'CE15': 0, 'CE16': 0, 'CE17': 0, 'CE18': 0, 'CE0': 0,
        'AA1': 0, 'AA2': 0, 'AA3': 0, 'AA4': 0, 'AA5': 0, 'AA6': 0, 'AA7': 0, 'AA8': 0, 'AA9': 0, 'AA10': 0, 'AA11': 0, 'AA12': 0, 'AA13': 0,
        'AA14': 0, 'AA15': 0, 'AA16': 0, 'AA0': 0,
        'CBM1': 0, 'CBM2': 0, 'CBM3': 0, 'CBM4': 0, 'CBM5': 0, 'CBM6': 0, 'CBM7': 0, 'CBM8': 0, 'CBM9': 0, 'CBM10': 0, 'CBM11': 0, 'CBM12': 0,
        'CBM13': 0, 'CBM14': 0, 'CBM15': 0, 'CBM16': 0, 'CBM17': 0, 'CBM18': 0, 'CBM19': 0, 'CBM20': 0, 'CBM21': 0, 'CBM22': 0, 'CBM23': 0,
        'CBM24': 0, 'CBM25': 0, 'CBM26': 0, 'CBM27': 0, 'CBM28': 0, 'CBM29': 0, 'CBM30': 0, 'CBM31': 0, 'CBM32': 0, 'CBM33': 0, 'CBM34': 0,
        'CBM35': 0, 'CBM36': 0, 'CBM37': 0, 'CBM38': 0, 'CBM39': 0, 'CBM40': 0, 'CBM41': 0, 'CBM42': 0, 'CBM43': 0, 'CBM44': 0, 'CBM45': 0,
        'CBM46': 0, 'CBM47': 0, 'CBM48': 0, 'CBM49': 0, 'CBM50': 0, 'CBM51': 0, 'CBM52': 0, 'CBM53': 0, 'CBM54': 0, 'CBM55': 0, 'CBM56': 0,
        'CBM57': 0, 'CBM58': 0, 'CBM59': 0, 'CBM60': 0, 'CBM61': 0, 'CBM62': 0, 'CBM63': 0, 'CBM64': 0, 'CBM65': 0, 'CBM66': 0, 'CBM67': 0,
        'CBM68': 0, 'CBM69': 0, 'CBM70': 0, 'CBM71': 0, 'CBM72': 0, 'CBM73': 0, 'CBM74': 0, 'CBM75': 0, 'CBM76': 0, 'CBM77': 0, 'CBM78': 0,
        'CBM79': 0, 'CBM80': 0, 'CBM81': 0, 'CBM82': 0, 'CBM83': 0, 'CBM84': 0, 'CBM85': 0, 'CBM86': 0, 'CBM87': 0, 'CBM88': 0, 'CBM0': 0
    }
    return foundation_dict
