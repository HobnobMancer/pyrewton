#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) University of St Andrews 2022
# (c) University of Strathclyde 2022
# (c) James Hutton Institute 2022
#
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
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
"""..."""


import logging
import json

import pandas as pd

from typing import List, Optional

from Bio import SeqIO
from saintBioutils.utilities.file_io import make_output_directory
from saintBioutils.utilities.logger import config_logger
from tqdm import tqdm


def main(args: Optional[List[str]] = None, logger: Optional[logging.Logger] = None):
    """Set up parser and logger, and coordinate prediction of CAZymes."""
    if logger is None:
        config_logger(args)
    logger = logging.getLogger(__package__)

    if str(args.output.parent) != '.':
        make_output_directory(args.output.parent, args.force, args.nodelete)
    
    dbcan_records = []  # dbCAN used to refer to any CAZyme classifier

    for record in SeqIO.parse(args.pyrewton_fasta, "fasta"):
        dbcan_records.append(record)

    cazy_records = []

    for record in SeqIO.parse(args.cazy_fasta, "fasta"):
        cazy_records.append(record)

    logger.warning(f"Loaded {len(dbcan_records)} and {len(cazy_records)} records from the CAZyme classifiers and CAZy fasta files")

    mmseq_clusters = parse_mmseq(args.mmseq_tsv)

    logger.warning(f"MMSeq2 found {len(list(mmseq_clusters.keys()))} clusters")

    write_out_summary(dbcan_records, cazy_records, mmseq_clusters, args.output)


def parse_mmseq(mmseq_tsv):
    """Parse mmseq output into a dict.
    
    :param mmseq_output: pandas df containing mmseq output
    
    Return dict"""
    mmseq_output = pd.read_table(mmseq_tsv)

    clusters = {}

    index = 0
    for index in tqdm(range(len(mmseq_output)), desc="Parsing MMseq tsv file"):
        row = mmseq_output.iloc[index]

        cluster_acc = row[0]
        member_acc = row[1]

        try:
            clusters[cluster_acc].add(member_acc)
        except KeyError:
            clusters[cluster_acc] = {member_acc}
    
    # check if the genbank accession used for the cluster name is in the cluster members
    # if not add it

    for cluster in clusters:
        if cluster not in clusters[cluster]:
            clusters[cluster_acc].add(cluster)

    return clusters


def write_out_summary(dbcan_records, cazy_records, mmseq_clusters, output_path):
    """Compile data for output files and write output files
    
    Data:
    (i) cluster size
    (ii) number of CAZymes from dbCAN
    (iii) number of CAZymes from CAZy
    (iv) percetage of CAZymes in cluster from dbCAN
    (v) percentage of CAZymes in cluster from CAZy
    
    Return nothing
    """
    info_cvs_data = []
    dbcan_csv_data = []

    info = ""
    dbcan_info = ""

    restructured_mmseq_clusters = {}  # {cluster: {'cazy':set, 'dbcan':set()}}

    for cluster in tqdm(mmseq_clusters, desc="Getting cluster summary info"):
        restructured_mmseq_clusters[cluster] = {'CAZy': set(), 'CAzymeClassifier': set()}

        parsed_records = set()

        cluster_members = mmseq_clusters[cluster]

        cluster_size = len(cluster_members)
        info += f"{cluster} cluster size: {cluster_size}\n"

        num_in_cazy = 0
        num_in_dbcan = 0
        for record in cazy_records:
            if record.id in cluster_members:
                num_in_cazy += 1
                restructured_mmseq_clusters[cluster]['CAZy'].add(record.id)
                parsed_records.add(record.id)
        for record in dbcan_records:
            if record.id in cluster_members:
                num_in_dbcan += 1
                restructured_mmseq_clusters[cluster]['CAzymeClassifier'].add(record.id)
                parsed_records.add(record.id)
            
        for protein_id in cluster_members:
            if protein_id not in parsed_records:
                print("WARNING: Could not identify source of protein:", protein_id)
        
        info += f"{cluster} cluster proteins in CAZy: {num_in_cazy}\n"
        info += f"{cluster} cluster proteins in CAzymeClassifier: {num_in_dbcan}\n"

        cazy_precent = (num_in_cazy / cluster_size) * 100
        dbcan_precent = (num_in_dbcan / cluster_size) * 100

        info += f"{cluster} cluster proteins in CAZy (%): {cazy_precent}\n"
        info += f"{cluster} cluster proteins in CAzymeClassifier (%): {dbcan_precent}\n"

        info_cvs_data.append([cluster, cluster_size, num_in_cazy, num_in_dbcan, cazy_precent, dbcan_precent])

        if num_in_dbcan != 0:
            dbcan_info += f"{cluster} cluster size: {cluster_size}\n"
            dbcan_info += f"{cluster} cluster proteins in CAZy: {num_in_cazy}\n"
            dbcan_info += f"{cluster} cluster proteins in CAzymeClassifier: {num_in_dbcan}\n"
            dbcan_info += f"{cluster} cluster proteins in CAZy (%): {cazy_precent}\n"
            dbcan_info += f"{cluster} cluster proteins in CAzymeClassifier (%): {dbcan_precent}\n"
            dbcan_csv_data.append([cluster, cluster_size, num_in_cazy, num_in_dbcan, cazy_precent, dbcan_precent])
    
    info_output_path = str(output_path) + ".txt"
    with open(info_output_path, "w") as fh:
        fh.write(info)
    
    dbcan_output_path = str(output_path) + "_dbcan_clusters.txt"
    with open(dbcan_output_path, "w") as fh:
        fh.write(dbcan_info)

    # convert sets to lists
    serialisable_clusters = {}
    for cluster in tqdm(restructured_mmseq_clusters, desc='Preparing clusters dict for serialisation'):
        serialisable_clusters[cluster] = {'CAZy': set(), 'CAZymeClassifier': set()}
        serialisable_clusters[cluster]['CAZy'] = list(restructured_mmseq_clusters[cluster]['CAZy'])
        serialisable_clusters[cluster]['CAZymeClassifier'] = list(restructured_mmseq_clusters[cluster]['CAZymeClassifier'])

    json_output_path = str(output_path) + "_clusters.json"
    with open(json_output_path, "w") as fh:
        json.dump(serialisable_clusters, fh)

    column_names = ["Cluster", "Cluster_size", "Num_in_CAZy", "Num_in_CAZymeClassifier", "Percentage_CAZy", "Percentage_CAZymeClassifier"]

    info_csv_path = str(output_path) + "_info.csv"
    info_df = pd.DataFrame(info_cvs_data, columns=column_names)
    info_df.to_csv(info_csv_path)

    dbcan_csv_path = str(output_path) + "_CAZymeClassifier_info.csv"
    dbcan_df = pd.DataFrame(dbcan_csv_data, columns=column_names)
    dbcan_df.to_csv(dbcan_csv_path)


if __name__ == "__main__":
    main()
