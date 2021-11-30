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

import json
import re
import logging
import sys

import pandas as pd

from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

from scraper.sql.sql_orm import CazyFamily, Genbank, Session
from saintBioutils.utilities.file_io import get_paths
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


@dataclass
class TestSet:
    """TestSet from a genomic assembly.

    :param fasta: Path, path to query fasta file
    :param alignment: Path, path to BLAST all vs all alignment score csv file
    :param tax_id: str, NCBI taxonomy ID of host species
    :param source: str, source of proteins within fasta file
    :param prediction_paths: path, path to which prediction tool output is written
    """

    fasta: Path  # path to FASTA file containing proteins for prediction
    alignment: Path
    tax_id: str  # NCBI taxonomy id, prefix NCBI:txid
    source: str  # source of protein sequences, genomic accession number
    prediction_paths: dict()  # contains path to outputs

    def __str__(self):
        """Representation of object"""
        return(
            (
                f"<tax-id={self.tax_id} source={self.source} "
                f"fasta={self.fasta} prediction_paths={self.prediction_paths}"
            )
        )


class ClassificationDF:
    """Represents a CAZyme/non-CAZyme annotation/classification df for a test set"""

    def __init__(self, genome_accession, df, df_path, testset_path, alignment_scores):
        self.genome_accession = genome_accession
        self.df = df
        self.df_path = df_path  # path to csv file
        self.testset_path = testset_path  # path to test set FASTA file
        self.alignment_scores = alignment_scores  # path to csv containing BLAST all vs. all scores

    def __str__(self):
        return f"<CAZyme/non-CAZyme classification df for test set {self.genome_accession}>"

    def __repr__(self):
        return f"<CAZyme/non-CAZyme classification df for test set {self.genome_accession}>"



def get_predictions(args):
    """Retrieve create class instances to represent output for each test set.

    :param args: cmd-line args parser

    Return list containing class instances, one instance per test set.
    """
    logger = logging.getLogger(__name__)

    # path_dict = {genomic_accession: {testset: Path, alingment: Path}}
    path_dict = get_testset_alignment_paths(args) 

    output_dirs = [f for f in (args.prediction_dir).iterdir() if f.is_dir()]  # store the paths to the output dirs

    predictions = []  # store TestSet instances

    for output_dir in tqdm(output_dirs, desc="Building TestSet instances"):
        # get the genomic assembly accession
        search_result = re.search(r"GCA_\d+?(_|.)\d+?\.\d+?.", str(output_dir), re.IGNORECASE)
        try:
            genomic_accession = search_result.group()[:-1]
        except AttributeError as e:
            logger.error(
                f"Could not get genomic asseccion from {output_dir}.\n"
                "Not including test set in the evaluation. The following error was raised:\n"
                f"{e}"
            )
            continue

        # get the taxonomy id
        search_result = re.search(r"txid\d+?_", str(output_dir), re.IGNORECASE)
        try:
            tax_id = search_result.group()[:-1]
        except AttributeError as e:
            logger.error(
                f"Could not get taxonomy ID from {output_dir}.\n"
                "Not including test set in the evaluation. The following error was raised:\n"
                f"{e}"
            )
            continue

        # build path to the FASTA file
        try:
            fasta_path = path_dict[genomic_accession]["testset"]
        except KeyError:
            logger.error(
                f"No test set FASTA file found for {genomic_accession}\n"
                "Excluding test set from evaluation"
            )
            continue
        try:
            alignment_path = path_dict[genomic_accession]["alignment"]
        except KeyError:
            logger.error(
                f"No BLAST all vs all alignment score csv file found for {genomic_accession}"
            )
            alignment_path = None

        # build dict of output paths
        output_paths = {}
        output_paths["dir"] = output_dir

        # build TestSet instance
        test_set = TestSet(fasta_path, alignment_path, tax_id, genomic_accession, output_paths)
        predictions.append(test_set)

    return predictions


def get_testset_alignment_paths(args):
    logger = logging.getLogger(__name__)

    # retrieve paths to all test set FASTA files
    testset_dir = args.testset_dir / "test_sets"
    all_test_set_path = get_paths.get_file_paths(
        testset_dir,
        prefixes=['genbank_proteins_'], 
        suffixes=['_test_set.fasta'],
    )

    # retrieve paths to all test set BLAST all versus all alignment score .csv files
    alignment_dir = args.testset_dir / "alignment_scores"
    all_align_score_paths = get_paths.get_file_paths(
        alignment_dir,
        prefixes=['genbank_proteins_'], 
        suffixes=['_alignment_scores.csv'],
    )

    path_dict = {}  # {genomic_accession: {testset: Path, alingment: Path}}

    for fasta_path in all_test_set_path:
        filename = fasta_path.name
        filename_parts = filename.split("_")
        genomic_accession = f"{filename_parts[3]}_{filename_parts[4]}"
        
        try:
            path_dict[genomic_accession]
            path_dict[genomic_accession]["testset"] = fasta_path
            logger.warning(f"Multiple test sets found with genomic accession: {genomic_accession}")
        except KeyError:
            path_dict[genomic_accession] = {}
            path_dict[genomic_accession]["testset"] = fasta_path
    
    for alignment_path in all_align_score_paths:
        filename = alignment_path.name
        filename_parts = filename.split("_")
        genomic_accession = f"{filename_parts[3]}_{filename_parts[4]}"

        try:
            path_dict[genomic_accession]
            path_dict[genomic_accession]["alignment"] = alignment_path
        except KeyError:
            logger.warning(f"No test set FASTA file found for {genomic_accession}")
            path_dict[genomic_accession] = {}
            path_dict[genomic_accession]["alignment"] = alignment_path
    
    return path_dict


def evaluate_performance(predictions, cazy, data_source, args):
    """Evaluate the performance of the CAZymes prediction tools.

    Binary classification: CAZyme/non-CAZyme differentiation
    Multilabel classification (MLC): CAZy family predictions

    :param predictions: list of TestSet class instances
    :param cazy: dict keyed by GenBank protein accession, valued by CAZy family classifications
        or open connection to a local CAZyme db engine
    :param data_source: str, 'dict' or 'db' depending if accessing CAZy annotations from a dict or db
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
        cazy,
        data_source,
        args,
    )

    output_path = args.output / f'binary_classification_evaluation_{time_stamp}.csv'
    binary_c_nc_statistics.to_csv(output_path)  # USED FOR EVALUATION IN R
    # Written in long for, with the columns: Statistic_parameters, Genomic_assembly, 
    # Prediciton_tool, Statistic_value. Statistical parameters calcualted are specificity, 
    # recall, precition, F1-score and accuracy

    # bootstrap resample binary CAZyme/non-CAZyme classifications, to evaluate performance range
    binary_classification.bootstrap_binary_c_nc_classifications(
        all_binary_c_nc_dfs,
        time_stamp,
        args,
    )

    binary_classification.get_f_pos_f_neg_predictions(
        all_binary_c_nc_dfs,
        time_stamp,
        args,
    )

    # write out ground truth CAZy class annotations to disk
    output_path = args.output / f"class_ground_truths_classifications_{time_stamp}.csv"
    class_ground_truths_df.to_csv(output_path)

    # write out ground truth CAZy family annotations to disk
    output_path = args.output / f"family_ground_truths_classifications_{time_stamp}.csv"
    all_family_ground_truths.to_csv(output_path)

    # Evaluate the performance of CAZy class predictions
    # Calculate ARI and RI for multilabel evaluation, and add to CAZy class prediction dataframe
    class_predictions_df = class_classifications.calculate_class_ari_ri(
        class_ground_truths_df,
        class_predictions_df,
        time_stamp,
    )

    # write out predicted CAZy class annotations and RI and ARI to disk
    output_path = args.output / f"class_predicted_classifications_{time_stamp}.csv"
    class_predictions_df.to_csv(output_path)   # USED FOR EVALUATION IN R

    # Calculate Fbeta, Sens, Spec and Acc for predicting each CAZy class
    # calculate the fbeta_score, sensitivity(recall), specificity and accuracy per CAZy class
    class_classifications.calculate_class_stats(
        class_ground_truths_df,
        class_predictions_df,
        time_stamp,
        args,
    )

    class_classifications.calculate_class_stats_by_testsets(
        class_ground_truths_df,
        class_predictions_df,
        time_stamp,
        args,
    )

    # Evaluate the performance of CAZy Family predictions
    # Calculate ARI and RI for multilabel evaluation, and add to CAZy family prediction dataframe
    all_family_predictions = family_classifications.calculate_family_ari_ri(
        all_family_predictions,
        all_family_ground_truths,
        time_stamp,
    )
    # write out predicted CAZy class annotations and RI and ARI to disk
    output_path = args.output / f"family_predicted_classifications_{time_stamp}.csv"
    all_family_predictions.to_csv(output_path)   # USED FOR EVALUATION IN R

    # evaluate the performance of predicting the correct CAZy family ACROSS ALL test sets
    fam_stats_df, fams_longform_df = family_classifications.calc_fam_stats(
        all_family_predictions,
        all_family_ground_truths,
        args,
    )  # creates a dataframe USED FOR EVALUATION IN R

    output_path = args.output / f"family_per_row_stats_{time_stamp}.csv"
    fam_stats_df.to_csv(output_path)

    output_path = args.output / f"family_long_form_stats_df_{time_stamp}.csv"
    fams_longform_df.to_csv(output_path)

    if args.tax_groups is not None:  # compare perforamnce between taxonomy groups
        evaluate_tax_group_performance(
            binary_c_nc_statistics,
            all_family_predictions,
            all_family_ground_truths,
            class_predictions_df,
            time_stamp,
            args,
        )

    return


def build_prediction_dataframes(predictions, time_stamp, cazy, data_source, args):
    """Parse prediction outputs from prediction tools and create dataframes for the stats eval.

    :param predictions: list of TestSet class instances
    :param time_stamp: str, data and time when evaluation started, used for naming files
    :param cazy: dict keyed by GenBank protein accession, valued by CAZy family classifications
        or open connection to a local CAZyme db engine
    :param data_source: str, 'dict' or 'db' depending if accessing CAZy annotations from a dict or db
    :param args: cmd-line arguments parser

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

        all_binary_classifications = []  # store all binary C/NC classifications for this test set
        # all_binary_classifications = [[protien_accession, tool, cazyme classification],]

        for tool in tqdm(standardised_outputs, "Adding raw predictions to cumlative dfs"):

            # parse binary CAZyme/non-CAZyme predictions for the prediction tool
            # retrieves list of nested lists [accession, tool, prediction (0 or 1)]
            all_binary_classifications = binary_classification.get_cazyme_noncazyme_predictions(
                test_set.source,
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
                predictions=standardised_outputs[tool],
                prediciton_tool=tool,
                cazy=cazy,
                data_source=data_source,
                genomic_accession=test_set.source,
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
            classifications_df,
            cazy,
            data_source,
        )

        # write out binary CAZyme/non-CAZyme predictions and ground truth annotations for test set
        # for documentation
        output_path = args.output / "binary_classifications" / f"binary_classifications_{time_stamp}_{test_set.source}.csv"
        classifications_df.to_csv(output_path)

        all_binary_c_nc_dfs.append(ClassificationDF(
            genome_accession=test_set.source,
            df=classifications_df,
            df_path=output_path,
            testset_path=test_set.fasta,
            alignment_scores=test_set.alignment,
        ))
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
            'Statistic_parameter',
            'Genomic_assembly',
            'Prediction_tool',
            'Statistic_value',
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


def get_fam_freq(args, cazy, timestamp, data_source):
    """Retrieve the frequencies of all CAZy families total occurences across all test sets.

    :param args: cmd-line args parser
    :param cazy: dict of CAZy family annotations of proteins from CAZy.
    :param timestamp: str, date-time script was invoked
    :param data_source: str, source of CAZy annotations, db or dict

    Return nothing.
    """
    logger = logging.getLogger(__name__)

    # get the paths to all test sets
    testset_dir = args.testset_dir / "test_sets"
    all_test_sets = get_paths.get_file_paths(
        testset_dir,
        prefixes=['genbank_proteins_'], 
        suffixes=['_test_set.fasta'],
    )

    logger.warning(f"Found {len(all_test_sets)} test sets in {args.testset_dir}")

    # build a dictionary to add frequency data to, key by CAZy fam, value by frequency
    freq_dict = family_classifications.foundation_dict()

    for testset in tqdm(all_test_sets, desc="Retrieving CAZy family freqs"):
        freq_dict = add_fam_freq(testset, freq_dict, cazy, data_source)

    # write out freq_dict
    output_path = args.output / f"CAZy_fam_testset_freq_{timestamp}.json"
    with open(output_path, "w") as fh:
        json.dump(freq_dict, fh)

    return


def add_fam_freq(testset, freq_dict, cazy, data_source):
    """Retrieve the frequency of CAZy family annotations in a test set.

    :param testset: Testset class instance
    :param freq_dict: dict of total CAZy family frequencies across all test sets
    :param cazy: dict of CAZy family annotations of proteins from CAZy.
    :param data_source: str, CAZy annotation source type, db or dict

    Return freq_dict.
    """
    logger = logging.getLogger(__name__)

    filename = str(testset).split("/")[-1]

    # open the FASTA file (FASTA file is the test set)
    with open(testset, "r") as fh:
        lines = fh.read().splitlines()
    
    for line in tqdm(lines, desc=f"Parsing {filename}"):
        if line.startswith(">"):
            # retrieve protein GenBank accession
            protein_accession = line.replace(">", "")
            
            # check if protein is a CAZy classified CAZyme
            if data_source == 'db':
                with Session(bind=cazy) as session:
                    fam_query = session.query(CazyFamily).\
                        join(CazyFamily, Genbank.families).\
                        filter(Genbank.genbank_accession == protein_accession).\
                        all()
                
                fam_annotations = [fam.family for fam in fam_query]
                if len(fam_annotations) == 0:  # non-CAZyme
                    continue

            else:
                try:
                    fam_annotations = cazy[protein_accession]
                except KeyError:
                    continue
                
            for fam in fam_annotations:
                try:
                    freq_dict[fam] += 1
                except KeyError:
                    logger.warning(
                        f"Retrieved '{fam}' from CAZy dict, but not included "
                        "in this modules foundation CAZy family dict"
                    )

    return freq_dict


def evaluate_tax_group_performance(
    binary_c_nc_statistics,
    all_family_predictions,
    all_family_ground_truths,
    class_predictions_df,
    time_stamp,
    args,
):
    """Evaluate the performance of the CAZyme classifiers between the taxonomy groups.
    
    :param binary_c_nc_statistics: Pandas df, of binary evaluation
        columns: Statistic_parameter, Genomic_assembly, Prediction_tool, Statistic_value
    :param all_family_predictions: Pandas df, 
        columns: Genomic_accession, Protein_accession, Prediction_tool, one column per CAZy family
            denoting if annotation predicted (1) or not predicted (0), Rand_index and Adjusted_rand_index
    :param all_family_ground_truths: Pandas df, 
        columns: Genomic_accession, Protein_accession, Prediction_tool, one column per CAZy family
            denoting if annotation given by CAZy (1) or not  (0)
    :param class_predictions_df: Pandas df,
        columns: Genomic_accession, Protein_accession, Prediction_tool, one column per CAZy class
            denoting if annotation predicted (1) or not predicted (0), Rand_index and Adjusted_rand_index
    :param time_stamp: str, date and time evaluation was invoked
    :param args: cmd-line args parser
    
    Return nothing.
    """
    logger = logging.getLogger(__name__)

    tax_dict = {}  # keyed by tax group name, valued by genomic accessions

    # add tax group column to binary classificaiton evaluation df
    binary_c_nc_statistics_tax = add_tax_group(binary_c_nc_statistics, tax_dict, 'Genomic_assembly')

    output_path = args.output / f'binary_classification_tax_comparison_{time_stamp}.csv'
    binary_c_nc_statistics_tax.to_csv(output_path)  # USED FOR EVALUATION IN R

    # add tax group column to family classification df
    all_family_predictions_tax = add_tax_group(all_family_predictions, tax_dict, 'Genomic_assembly')

    output_path = args.output / f"family_classification_tax_comparison_{time_stamp}.csv"
    all_family_predictions_tax.to_csv(output_path)  # USED FOR EVALUATION IN R



    # add tax group column to family classification df
    class_predictions_df_tax = add_tax_group(class_predictions_df, tax_dict, 'Genomic_assembly')

    output_path = args.output / f"class_classification_tax_comparison_{time_stamp}.csv"
    class_predictions_df_tax.to_csv(output_path)  # USED IN R EVALUATION


def add_tax_group(df, tax_dict, column_name):
    """Add tax_group column to and existing dataframe
    
    :param df: Pandas df
    :param tax_dict: dict, keyed by tax group name, values by list of genomic accessions
    :param column_name: str, name of column containing the genomic accessions
    
    return df with added tax_group column
    """
    logger = logging.getLogger(__name__)

    # add on taxonomy group column to the binary_c_nc_statistics df
    tax_groups = []  # new content for the tax_group column
    index = 0
    for index in tqdm(range(len(df[column_name])), desc='Adding tax group to df'):
        row = df.iloc[index]
        genomic_accession = row[column_name]
        tax_group = None
        for taxa in tax_dict:
            tax_genomic_accs = tax_dict[taxa]
            if genomic_accession in tax_genomic_accs:
                tax_group = taxa

        if tax_group is None:
            logger.warning(
                f'Accession {genomic_accession} retrieved from test sets but not included in tax group data\n'
                'Setting tax_group as NaN for test set'
            )
            tax_groups.append(tax_group)
        else:
            tax_groups.append(tax_group)
    
    df['Tax_group'] = tax_groups

    return df
