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
"""This script is for evaluating the prediction of CAZy family annotations."""


import logging

import numpy as np
import pandas as pd

from tqdm import tqdm

from scraper.sql.sql_orm import CazyFamily, Genbank, Session
from sklearn.metrics import fbeta_score, confusion_matrix, recall_score, precision_score
from sklearn.metrics.cluster import rand_score, adjusted_rand_score

from pyrewton.utilities import build_logger


def get_family_classifications(predictions, prediciton_tool, cazy, data_source, genomic_accession):
    """Retrieve the predicted  and ground truth (CAZy determined) CAZy family annotations.

    :param predictions: list of CazymeProteinPrediction instances 
    :param prediction_tools: str, name of the prediction tool being passed
    :param cazy: dict keyed by GenBank protein accession, valued by CAZy family classifications
        or open connection to a local CAZyme db engine
    :param data_source: str, 'dict' or 'db' depending if accessing CAZy annotations from a dict or db
    :param genomic_accession: str, genomic accession of genomic assembly from which proteins are sourced

    Return two lists, one of predicted results, one of ground truths.
    Each list is a list of nested lists. Each nested list contains:
    [genomic_accession, protein_accession, prediction_tool, one column for each fam with 
    its 0 (not predicted) or 1 (predicted) prediction score].
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

        # add cazy family (0/1) predictions to list representing known CAZy family annotations
        new_predictions += list(fam_predictions.values())

        if data_source == 'db':
            with Session(bind=cazy) as session:
                fam_query = session.query(CazyFamily).\
                    join(CazyFamily, Genbank.families).\
                    filter(Genbank.genbank_accession == protein.protein_accession).\
                    all()

            for fam in fam_query:
                try:
                    known_fams[fam.family] = 1
                except KeyError as e:
                    logger.error(
                        f"KeyError raised for known CAZy family {fam.family}. Error raised:\n{e}"
                    )

        else:
            # retrieve the ground thruth (CAZy determined) CAZy family annotations for the protein
            try:
                cazy_annotations = cazy_dict[protein.protein_accession]
                try:
                    for fam in cazy_annotations:
                        fam = fam.strip()
                        known_fams[fam] = 1
                except KeyError as e:
                    logger.error(
                        f"KeyError raised for known CAZy family {fam}. Error raised:\n{e}"
                    )
            except KeyError:
                pass

        # add cazy family (0/1) ground truth annotations to list representing all 
        # CAZy family ground truth annotations for the current working protein
        new_ground_truths += list(known_fams.values())

        all_predicted_annotations.append(new_predictions)
        all_ground_truths.append(new_ground_truths)

    return all_predicted_annotations, all_ground_truths


def calculate_family_ari_ri(prediction_df, ground_truth_df, time_stamp):
    """Calculate the adjusted rand index (ARI) and rand index (RI) for every protein.

    :param prediction_df: pandas df of predicted CAZy family annotations
    :param ground_truth_df: pandas df of CAZy determined CAZy family annotations
    :param time_stamp: str, used to create file names

    Return predictions_df with RI and ARI added in.
    """
    ri_scores = []
    ari_scores = []

    family_names = foundation_dict()
    family_names = list(family_names.keys())

    row_index = 0
    for row_index in tqdm(
        range(len(prediction_df["Prediction_tool"])), desc="Calculating CAZy family ARI and RI",
    ):

        ground_truths_row = ground_truth_df.iloc[row_index]
        predictions_row = prediction_df.iloc[row_index]

        y_true = list(ground_truths_row[family_names])
        y_pred = list(predictions_row[family_names])

        ri = rand_score(y_true, y_pred)
        ri_scores.append(ri)

        ari = adjusted_rand_score(y_true, y_pred)
        ari_scores.append(ari)

    prediction_df["Rand_index"] = ri_scores
    prediction_df["Adjusted_Rand_index"] = ari_scores

    return prediction_df


def calc_fam_stats(predictions_df, ground_truths_df, time_stamp, args):
    """Calculate the Specificity, Sensitivity, Precision, Fbeta, and accuracy for each CAZy family.

    :param predictions_df: df of predicted CAZy family annotations from all prediction tools
    :param ground_truths_df: df of CAZy annotations of proteins
    :param time_stamp: str, time evaluation was started
    :param args: cmd_line args parser

    Return nothing, instead write out the created dataframes to disk.
    """
    specific_logger = build_logger(args.output, "cazy_family_performance.log")
    logger = logging.getLogger(__name__)

    stats_orientated_df_data = []  # [[CAZyFam, Prediction_Tool, Specificity, Recall, Fbeta, Accuracy]]
    long_dataframe_data = []  # [[CAZyFam, StatParameter, PredicitonTool, StatValue]]

    family_names = foundation_dict()
    family_names = list(family_names.keys())

    for tool in ["dbCAN", "HMMER", "Hotpep", "DIAMOND", "CUPP", "eCAMI"]:
        # retrieve the relevant rows for the prediction tool
        grnd_trth_df = ground_truths_df.loc[ground_truths_df["Prediction_tool"] == tool]
        pred_df = predictions_df.loc[predictions_df["Prediction_tool"] == tool]

        # build an empty dataframe to store all predictions EXCEPT true negative non-CAZymes
        tp_fp_fn_ground_truths = pd.DataFrame(columns=list(grnd_trth_df.columns))
        tp_fp_fn_predictions = pd.DataFrame(columns=list(pred_df.columns))

        # exclude true negative non-CAZyme predictions
        index = 0
        for index in tqdm(
            range(len(grnd_trth_df["Prediction_tool"])),
            desc=f"Removing TN non-CAZyme predictions from CAZy fam predictions for {tool}",
        ):
            y_true = grnd_trth_df.iloc[index]
            y_true = list(y_true[family_names])  # retrieve only the cazy family 0/1 annotations
            
            y_pred = pred_df.iloc[index]
            y_pred = list(y_pred[family_names])  # retrieve only the cazy family 0/1 annotations

            if (1 not in y_true) and (1 not in y_pred):
                # if y_true and y_pred are all 0s, this is a true negative non-CAZyme prediction
                # do not include true negative non-CAZyme predictions
                continue
            else:
                # add TP, FP or FN result to the dataframes
                tp_fp_fn_ground_truths = tp_fp_fn_ground_truths.append(grnd_trth_df.iloc[index])
                tp_fp_fn_predictions = tp_fp_fn_predictions.append(pred_df.iloc[index])

        # iterate through the families and evaluate the CAZy family performance
        for fam in tqdm(family_names, desc=f"Evaluating performance per CAZy family for {tool}"):
            y_true = list(tp_fp_fn_ground_truths[fam])
            y_pred = list(tp_fp_fn_predictions[fam])

            # check if family was predicted and was included in test set as a known annotations
            if (1 not in y_true) and (1 not in y_pred):
                # do not include in statistics
                logger.warning(
                    f"{fam} not predicted by {tool} and not in known annotations\n"
                    f"Excluding {fam} from evaluation by setting all stats results as NaN"
                )
                specific_logger.warning(
                    f"{fam} not predicted by {tool} and not in known annotations. "
                    f"Excluding {fam} from evaluation by setting all stats results as NaN"
                )
                specificity = np.nan
                recall = np.nan
                precision = np.nan
                fbeta = np.nan
                accuracy = np.nan

                # [[CAZyFam, Prediction_Tool, Specificity, Recall, Fbeta, Accuracy]]
                stats_orientated_row = [
                    fam,
                    tool,
                    specificity,
                    recall,
                    precision,
                    fbeta,
                    accuracy,
                ]
                stats_orientated_df_data.append(stats_orientated_row)

                for stat in [
                    [specificity, "Specificity"],
                    [recall, "Recall"],
                    [precision, "Precision"],
                    [fbeta, "Fbeta_score"],
                    [accuracy, "Accuracy"],
                ]:
                    # [[CAZyFam, StatParameter, PredicitonTool, StatValue]]
                    longform_row = [fam, stat[1], tool, stat[0]]
                    long_dataframe_data.append(longform_row)

                continue

            # else: calculate performance for the CAZy family

            # recall aka sensitivity
            recall = recall_score(y_true, y_pred)
            long_dataframe_data.append([fam, "Recall", tool, recall])

            precision = precision_score(y_true, y_pred)
            long_dataframe_data.append([fam, "Precision", tool, recall])

            fbeta = fbeta_score(y_true, y_pred, beta=args.beta)
            long_dataframe_data.append([fam, "Fbeta_score", tool, fbeta])

            # create confusion matrix for calculating Specificty and Accuracy
            cm = confusion_matrix(y_true, y_pred)
            try:
                tn = cm[0][0]
                fn = cm[1][0]
                tp = cm[1][1]
                fp = cm[0][1]
            except IndexError as e:
                logger.warning(
                    f"Error raised when creating confusion matrix for {tool}: {fam},\n"
                    f"Error raised:\n{e}"
                )
                specific_logger.warning(
                    f"Error raised when creating confusion matrix for {tool}: {fam},\n"
                    f"Error raised:\n{e}\n"
                    f"Prediction Tool: {tool}\tFamily: {fam}\tStat: building confusion matrix\n"
                    f"y_true: {y_true}\ny_pred: {y_pred}\n"
                )

                specificity = np.nan
                long_dataframe_data.append([fam, "Specificity", tool, specificity])

                accuracy = np.nan
                long_dataframe_data.append([fam, "Accuracy", tool, accuracy])

                # [[CAZyFam, Prediction_Tool, Specificity, Recall, Fbeta, Accuracy]]
                stats_orientated_row = [
                    fam,
                    tool,
                    specificity,
                    recall,
                    precision,
                    fbeta,
                    accuracy,
                ]
                stats_orientated_df_data.append(stats_orientated_row)

            # calculate specificity and accuracy
            specificity = tn / (tn + fp)
            long_dataframe_data.append([fam, "Specificity", tool, specificity])

            accuracy = (tp + tn)/(tp + fp + fn + tn)
            long_dataframe_data.append([fam, "Accuracy", tool, accuracy])

            # [[CAZyFam, Prediction_Tool, Specificity, Recall, Fbeta, Accuracy]]
            stats_orientated_row = [
                fam,
                tool,
                specificity,
                recall,
                precision,
                fbeta,
                accuracy,
            ]
            stats_orientated_df_data.append(stats_orientated_row)

    # build statistics orientated datafrae
    stats_df = pd.DataFrame(
        stats_orientated_df_data,
        columns=[
            "CAZy_family",
            "Prediction_tool",
            "Specificity",
            "Recall",
            "Precision",
            "Fbeta_score",
            "Accuracy",
        ],
    )

    output_path = args.output / f"cazy_fam_stats_fam_per_row_{time_stamp}.csv"
    stats_df.to_csv(output_path)

    # build long form dataframe
    longform_df = pd.DataFrame(
        long_dataframe_data,
        columns=["CAZy_family", "Stat_parameter", "Prediction_tool", "Stat_value"],
    )

    output_path = args.output / f"cazy_fam_long_form_stats_df_{time_stamp}.csv"
    longform_df.to_csv(output_path)

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
