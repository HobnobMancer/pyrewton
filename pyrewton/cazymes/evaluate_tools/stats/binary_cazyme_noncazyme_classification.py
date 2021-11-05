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
"""This script is for evaluating the binary CAZyme/non-CAZyme prediction performances"""

import logging
import random

import pandas as pd
import numpy as np

from tqdm import tqdm

from scraper.sql.sql_orm import Genbank, Session
from sklearn.metrics import confusion_matrix, fbeta_score, precision_score, recall_score
from sklearn.utils import resample


def get_cazyme_noncazyme_predictions(tool_predictions, prediction_tool, classifications):
    """Add CAZyme/non-CAZyme predictions to list of all CAZyme/non-CAZyme predictions.

    :param tool_predictions: list of CazymeProteinPrediction instances
    :param prediction_tool: str, name of the prediction tool being passed
    :param classifications: dict of CAZyme/non-CAZyme predicted classifications

    Return list of CAZyme/non-CAZymes classifications, written in long form to match TidyData
    requirements and recommendations in R.
    """
    # classifications  = [[protien_accession, tool, cazyme classification],]
    if (prediction_tool == "eCAMI") or (prediction_tool == "CUPP"):
        for prediction in tqdm(tool_predictions, desc=f"Retreiving {prediction_tool} predictions"):
            accession = prediction.protein_accession.split(" ")[0]
            classifications.append(
                [
                    accession,
                    prediction_tool,
                    prediction.cazyme_classification,
                ],
            )
        return classifications

    for prediction in tqdm(tool_predictions, desc=f"Retreiving {prediction_tool} predictions"):
        classifications.append(
            [
                prediction.protein_accession,
                prediction_tool,
                prediction.cazyme_classification
            ],
        )
    return classifications


def build_classification_df(predicted_classifications):
    """Build dataframe of classifications. Unique protein per row, unqiue tool per column.

    :param classifications: list of lists, each list contains a protein accession, name of the
        prediction tool, and the CAZyme (represented by a value of 1) or non-CAZyme (represented
        by a value of 0) classification.

    Return pandas dataframe.
    """
    # create an empty dictionary to store the data for building the dataframe
    # this dictionary is used to then organise the data into a list of lists, for creating
    # a dataframe in long/tidy form below, with a unique protein per row, a the corresponding
    # prediction tool classification per column
    classification_dict = {}  # { protein_acession : {prediction_tool : [all_classiciations]} }

    for classification in predicted_classifications:
        protein_accession = classification[0]
        prediciton_tool = classification[1]
        cazyme_classification = classification[2]
        try:
            classification_dict[protein_accession]

            try:
                classification_dict[protein_accession][prediciton_tool]
                classification_dict[protein_accession][prediciton_tool].append(
                    cazyme_classification
                )
            except KeyError:
                classification_dict[protein_accession][prediciton_tool] = [cazyme_classification]

        except KeyError:
            classification_dict[protein_accession] = {}
            classification_dict[protein_accession][prediciton_tool] = [cazyme_classification]

    # build a series of lists
    data = []  # [Protein_accession, dbcan, hmmer, hotpep, diamond, cupp, ecami]
    for protein_accession in classification_dict:
        for i in range(len(classification_dict[protein_accession]["dbCAN"])): 
            # for loop is included in case of duplicate proteins in the test set
            # create a unique row for each instance of a protein in the test set
            new_data = [protein_accession]
            for tool in ["dbCAN", "HMMER", "Hotpep", "DIAMOND", "CUPP", "eCAMI"]:
                new_data.append(classification_dict[protein_accession][tool][i])
            data.append(new_data)

    df = pd.DataFrame(
        data,
        columns=[
            "Protein_accession",
            "dbCAN",
            "HMMER",
            "Hotpep",
            "DIAMOND",
            "CUPP",
            "eCAMI",
        ],
    )
    df = df.set_index("Protein_accession")

    return df


def add_ground_truths(classifications_df, cazy, data_source):
    """Retrieve ground truth CAZyme/non-CAZyme classifications and add to the df of predictions.

    :param classifications_df: pandas dataframe of prediction tool CAZyme/non-CAZyme classifications
    for each protein. Each unqiue protein is one row, each prediction tool is one column.
    :param cazy: dict or sqlalchemy connection to a local SQLite3 db engine

    Return df containing the prediction tool predictions and CAZy ground truths added as an
    additional column, called 'CAZy'.
    """
    cazy_c_nc_classification = []
    # tool_predictions = clasifications_df
    protein_accessions = classifications_df.index

    if data_source == 'db':
        for protein_accession in protein_accessions:
            with Session(bind=cazy) as session:
                gbk_query = session(Genbank).\
                    filter(Genbank.genbank_accession == protein_accession).\
                    all()
            
            if len(gbk_query) != 0:
                cazyme_classification = 1
            else:
                cazyme_classification = 0

            cazy_c_nc_classification.append(cazyme_classification)
    
    else:
        for protein_accession in protein_accessions:
            try:
                cazy[protein_accession]
                cazyme_classification = 1
            except KeyError:
                cazyme_classification = 0

            cazy_c_nc_classification.append(cazyme_classification)

    classifications_df["CAZy"] = cazy_c_nc_classification

    return classifications_df


def evaluate_binary_cazyme_noncazyme_predictions(classifications_df, genomic_accession, args):
    """Calculate the Fbeta-score, recall (sensitivity), precision and specificity for
    CAZyme/non-CAZyme prediction.

    The value of beta for the Fbeta-score is defined via the cmd-line

    :param classifications_df: pandas df of CAZyme/non-CAZyme classifications
    :param genomic_accession: str, accession of the genomic assembly from which the test set is from
    :param args: cmd-line args parser

    Return list of lists, each list contains the name of statistical parameter, the genomic
    accession of the source genome for the test set, the prediction tool and the value
    of the statistical parameter.
    """
    y_true = classifications_df['CAZy'].to_numpy()  # convert to numpy array

    stats_results = []  # [[stat, genome, tool, value]]

    # build a dict to store statistical parameters
    for tool in ["dbCAN", "HMMER", "Hotpep", "DIAMOND", "CUPP", "eCAMI"]:
        y_pred = classifications_df[tool].to_numpy()

        # calculations
        cm = confusion_matrix(y_true, y_pred)
        tn = cm[0][0]
        fn = cm[1][0]
        tp = cm[1][1]
        fp = cm[0][1]
        specificity_result = tn / (tn + fp)
        stats_results.append(["Specificity", genomic_accession, tool, specificity_result])

        recall_result = recall_score(y_true, y_pred)
        stats_results.append(["Recall", genomic_accession, tool, recall_result])

        precision_result = precision_score(y_true, y_pred)
        stats_results.append(["Precision", genomic_accession, tool, precision_result])

        fbeta_result = fbeta_score(y_true, y_pred, beta=args.beta)
        stats_results.append(["FBeta-score", genomic_accession, tool, fbeta_result])

        accuracy = (tp + tn)/(tp + fp + fn + tn)
        stats_results.append(["Accuracy", genomic_accession, tool, accuracy])

    return stats_results


def bootstrap_binary_c_nc_classifications(all_binary_dfs, time_stamp, args):
    """Bootstrap resample binary CAZyme/non-CAZyme classifications to evaluated the expected
    range in performance of each prediciton tool.

    :param all_binary_dfs: list of ClassificationDFs, each one is a pandas df with asscoiated meta
    data
    :param time_stamp: str, time when evaluation started
    :param args: cmd-line args parser

    Return nothing.
    """
    # Perform bootstrapping of CAZyme/non-CAZyme predictions
    bootstrap_results = []
    #  bootstrap_results = [[genomic_accession, prediction_tool, bootstrap_number, accuracy_score]]

    binary_dfs = random.sample(all_binary_dfs, args.bs_sample_size)

    for tool in ["dbCAN", "HMMER", "Hotpep", "DIAMOND", "CUPP", "eCAMI"]:
        for df in tqdm(binary_dfs, desc=f"Bootstrapping accuracies for {tool}"):
            temp_df = pd.concat([df.df[tool], df.df["CAZy"]], axis=1, keys=None)

            # bootstrap the accuracy
            bs_acc = bootstrap_acc(temp_df, args.bs_resampling)
            bs_number = 1
            for result in bs_acc:
                bootstrap_results.append([df.genome_accession, tool, bs_number, result])
                bs_number += 1

    bootstrap_results = pd.DataFrame(
        bootstrap_results,
        columns=['Genomic_accession', 'Prediction_tool', 'Bootstrap_number', 'accuracy'],
    )

    output_path = args.output / f'bootstrap_accuracy_evaluation_{time_stamp}.csv'
    bootstrap_results.to_csv(output_path) # USED FOR EVALUATION IN R

    return


def calc_acc(y_true, y_pred):
    """Calculate accuracy.

    :param y_true: ground truths, pandas series
    :param y_pred: predictions, pandas series

    Return float.
    """
    cm = confusion_matrix(y_true, y_pred)
    tn = cm[0][0]
    fn = cm[1][0]
    tp = cm[1][1]
    fp = cm[0][1]
    accuracy = (tp + tn)/(tp + fp + fn + tn)

    return accuracy


def bootstrap_acc(df, n):
    """Bootstrap accuracy.

    :param df: Pandas dataframe. 2 columns. First is predicted, second is ground truth
    :param n: int, number of times to resample

    Return tuple (bs.lower, bs.median, bs.upper).
    """
    bootstraps = []
    df_array = df.to_numpy()
    for _ in range(n):
        # create sample with replacement
        sample = resample(df_array)
        bootstraps.append(calc_acc(sample[:, 1], sample[:, 0]))

    return bootstraps


def get_f_pos_f_neg_predictions(all_binary_c_nc_dfs, time_stamp, args):
    """Identify all false positive and false negative CAZyme and non-CAZyme predictions.

    Builds one dataframe of false positive CAZyme predictions, and another of false negative
    non-CAZyme predictions. Both dataframes also include the number of tools that produced the
    false positive/negative prediction.

    :param all_binary_c_nc_dfs: list of ClassificationDF class instances
    :param time_stamp: str, date/time script was enabled
    :param args: cmd-line args

    Return nothing.
    """

    df_data = {
        "Protein_accession":[],
        "dbCAN":[],
        "HMMER":[],
        "Hotpep":[],
        "DIAMOND":[],
        "CUPP":[],
        "eCAMI":[],
        "CAZy":[],
        "Number_of_tools": [],
        "Blast_score_ratio": [],
    }

    false_positives_df_data = df_data
    false_negatives_df_data = df_data

    for classification_df_instance in all_binary_c_nc_dfs:
        classification_df = classification_df_instance.df

        # protein accessions are the row indexes
        classification_df.reset_index(inplace=True)
        classification_df = classification_df.rename(columns = {'index':'Protein_accession'})

        row_index = 0
        for row_index in tqdm(
            range(len(classification_df["CAZy"])),
            desc="Identifying FP and FN binary predictions",
        ):
            row = classification_df.iloc[row_index]

            if row["CAZy"] == 0:  # non-CAZyme (according to CAZy)
                tool_total = (
                    row["dbCAN"] +
                    row["HMMER"] +
                    row["DIAMOND"] +
                    row["Hotpep"] +
                    row["CUPP"] +
                    row["eCAMI"]
                )

                if tool_total != 0:  # found at least one false positive prediction
                    false_positives_df_data["Protein_accession"].append(row["Protein_accession"])
                    false_positives_df_data["dbCAN"].append(row["dbCAN"])
                    false_positives_df_data["HMMER"].append(row["HMMER"])
                    false_positives_df_data["DIAMOND"].append(row["DIAMOND"])
                    false_positives_df_data["Hotpep"].append(row["Hotpep"])
                    false_positives_df_data["CUPP"].append(row["CUPP"])
                    false_positives_df_data["eCAMI"].append(row["eCAMI"])
                    false_positives_df_data["CAZy"].append(row["CAZy"])
                    false_positives_df_data["Number_of_tools"].append(tool_total)

                    # retrieve the highest blast score ratio between the protein an a known CAZyme
                    false_positives_df_data["Blast_score_ratio"].append(
                        get_blast_score_ratio(
                            classification_df_instance,
                            row["Protein_accession"],
                            "FP"
                        ),
                    )

            else:  # CAZyme (according to CAZy)
                tool_total = (
                    row["dbCAN"] +
                    row["HMMER"] +
                    row["DIAMOND"] +
                    row["Hotpep"] +
                    row["CUPP"] +
                    row["eCAMI"]
                )

                if tool_total != 6:  # found at least one false negative prediction
                    false_negatives_df_data["Protein_accession"].append(row["Protein_accession"])
                    false_negatives_df_data["dbCAN"].append(row["dbCAN"])
                    false_negatives_df_data["HMMER"].append(row["HMMER"])
                    false_negatives_df_data["DIAMOND"].append(row["DIAMOND"])
                    false_negatives_df_data["Hotpep"].append(row["Hotpep"])
                    false_negatives_df_data["CUPP"].append(row["CUPP"])
                    false_negatives_df_data["eCAMI"].append(row["eCAMI"])
                    false_negatives_df_data["CAZy"].append(row["CAZy"])
                    false_negatives_df_data["Number_of_tools"].append(tool_total)

                    # retrieve the highest blast score ratio between the protein an a known CAZyme
                    false_negatives_df_data["Blast_score_ratio"].append(
                        get_blast_score_ratio(
                            classification_df_instance,
                            row["Protein_accession"],
                            "FN"
                        ),
                    )
    
    false_positives_df = pd.DataFrame(df_data)
    false_negatives_df = pd.DataFrame(df_data)

    # write out the dataframes
    output_path = args.output / f"binary_false_positive_predictions_{time_stamp}.csv"
    false_positives_df.to_csv(output_path)

    output_path = args.output / f"binary_false_negative_predictions_{time_stamp}.csv"
    false_negatives_df.to_csv(output_path)

    return


def get_blast_score_ratio(classification_df_instance, protein_accession, prediction_type):
    """Retrieve the highest BLAST score ratio between the given protein and a known CAZyme.

    :param classification_df_instance: ClassificationDF instance.
    :param protein_accession: str, GenBank protein accession
    :param prediction_type: str, 'FP' or 'FN'

    Return the BLAST score ratio (str).
    """
    logger = logging.getLogger(__name__)

    # compile the path to the alignment scores created when making the test set
    alignment_score_file = classification_df_instance.testset_path
    alignment_score_file = str(alignment_score_file).replace(
        "_test_set.fasta",
        "_alignment_scores.csv",
    )

    try:
        # load in the data
        alignment_score_df = pd.read_csv(alignment_score_file)

    except FileNotFoundError:
        logger.error(
            f"Could not find alignment score file for\n{classification_df_instance.df_path}\n"
            "Adding BLAST score ratio to False Positive/Negative dataframe as null/nan"
        )
        return np.nan

    # retrieve all rows containing data for the given protein
    if prediction_type is "FN":
        protein_blast_scores_df = alignment_score_df[alignment_score_df["subject (Cazyme)"] == protein_accession]
    else:
        protein_blast_scores_df = alignment_score_df[alignment_score_df["query (non-CAZyme)"] == protein_accession]

    protein_blast_scores_df = protein_blast_scores_df.sort_values(by = "blast_score_ratio")
    highest_score = protein_blast_scores_df.iloc[-1]
    highest_score = highest_score["blast_score_ratio"]

    return highest_score
