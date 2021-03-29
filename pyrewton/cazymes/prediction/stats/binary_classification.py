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


import numpy as np
import pandas as pd

from statistics import mean, stdev, median
from tqdm import tqdm

from sklearn.metrics import confusion_matrix, f1_score, precision_score, recall_score
from sklearn.utils import resample


def retrieve_ground_truths(tool_predictions, cazy_dict):
    """Retrieve dict of ground truth CAZyme/non-CAZyme predictions. {accession: {tool: prediction}}
    
    :param tool_predicitons: list of CazymeProteinPrediction instances.
    
    Return dict of ground truths {accession: {"CAZy": prediction(0/1)}}"""
    y_true = {}

    for prediction in tqdm(tool_predictions, desc="Retrieving ground truths"):
        # get the CAZy CAZyme/non-CAZyme classification
        try:
            cazy_dict[prediction.protein_accession]
            cazyme_classification = 1
        except KeyError:
            cazyme_classification = 0

        y_true[prediction.protein_accession] = {"CAZy": cazyme_classification}

    return y_true


def get_prediction_classifications(tool_predictions, prediction_tool, classifications):
    """Retrieve dict of CAZyme prediction tool predictions of CAZyme/non-CAZyme classification.
    
    :param tool_predictions: list of CazymeProteinPrediction instances
    :param prediction_tool: str, name of the prediction tool being passed
    :param classifications: dict of CAZyme/non-CAZyme predicted classifications
    
    Return dict of CAZyme(1)/non-CAZyme(0) classifications.
    """
    if (prediction_tool == "eCAMI") or (prediction_tool == "CUPP"):
        for prediction in tqdm(tool_predictions, desc=f"Retreiving {prediction_tool} predictions"):
            accession = prediction.protein_accession.split(" ")[0]
            try:
                classifications[accession]
                classifications[accession][prediction_tool] = prediction.cazyme_classification
            except KeyError:
                print(f"WARNING -- no ground truth for {prediction.protein_accession} from {prediction_tool}")
        return classifications

    for prediction in tqdm(tool_predictions, desc=f"Retreiving {prediction_tool} predictions"):
        try:
            classifications[prediction.protein_accession]
            classifications[prediction.protein_accession][prediction_tool] = prediction.cazyme_classification
        except KeyError:
            print(f"WARNING -- no ground truth for {prediction.protein_accession} from {prediction_tool}")
        
    return classifications


def build_classification_df(classifications):
    """Build dataframe of classifications. Unique protein per row, unqiue tool per column.
    
    :param classifications: dict of known and predicted classificaitons.
    
    Return pandas dataframe.
    """
    classification_df = pd.DataFrame.from_dict({(i): classifications[i]
                       for i in classifications.keys() 
                       for j in classifications[i].keys()},
                       orient='index')
    return classification_df


def evaluate_binary_classification(classifications_df, genomic_accession):
    """Calculate the F1-score, recall (sensitivity), precision and specificity for CAZyme/non-CAZyme prediction.
    
    :param classifications_df: pandas df of CAZyme/non-CAZyme classifications
    :param genomic_accession: str, accession of the genomic assembly from which the test set is from
    
    Return F1-score, recall (sensitivity), precision and specificity, and summary dict.
    """
    y_true = classifications_df['CAZy'].to_numpy()  # convert to numpy array
    
    summary_lists = []  # [[stat, genome, tool, value]]
    
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
        summary_lists.append(["Specificity", genomic_accession, tool, specificity_result])

        recall_result = recall_score(y_true, y_pred)
        summary_lists.append(["Recall", genomic_accession, tool, recall_result])

        precision_result = precision_score(y_true, y_pred)
        summary_lists.append(["Precision", genomic_accession, tool, precision_result])
    
        f1_result = f1_score(y_true, y_pred)
        summary_lists.append(["F1-score", genomic_accession, tool, f1_result])
        
        accuracy = (tp + tn)/(tp + fp + fn + tn)
        summary_lists.append(["Accuracy", genomic_accession, tool, accuracy])

    return summary_lists


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


def bootstrap_acc(y_true, y_pred, n):
    """Bootstrap accuracy.

    :param y_true: ground truths, pandas series
    :param y_pred: predictions, pandas series
    :param n: int, number of times to resample
    
    Return tuple (bs.lower, bs.median, bs.upper).
    """
    bootstraps = []
    pred_array = y_pred.to_numpy()
    true_array = y_true.to_numpy()
    for _ in range(n):
        sample = resample(pred_array)
        bootstraps.append(calc_acc(true_array, sample))

    return min(bootstraps), median(bootstraps), max(bootstraps)
