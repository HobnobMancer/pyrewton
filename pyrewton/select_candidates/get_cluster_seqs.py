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
"""Gather seqs for each cluster, and write to multiseq FASTA file per cluster"""


import logging
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
    
    sequences = load_seqs(args, logger)

    clusters = parse_mmseq(args.mmseq_tsv)  # {cluster name: {members ids}}

    write_cluster_seqs(sequences, clusters, args, logger)


def load_seqs(args, logger):
    seqs = {}  # id : seq
    for record in SeqIO.parse(args.fasta, "fasta"):
        seqs[record.id] = record
    
    logger.warning(f"Loaded {len(list(seqs.keys()))} sequences into memory")

    return seqs


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


    seq_record_ids = set()
    seq_records = []
    i = 0
    for record in SeqIO.parse(args.fasta_1, "fasta"):
        record.description=""
        if record.id not in seq_record_ids:
            seq_records.append(record)
            seq_record_ids.add(record.id)
        i += 1
    logger.warning(f"Loaded {i} sequences from: {args.fasta_1}")
    i = 0
    for record in SeqIO.parse(args.fasta_2, "fasta"):
        record.description=""
        if record.id not in seq_record_ids:
            seq_records.append(record)
            seq_record_ids.add(record.id)
        i += 1
    logger.warning(f"Loaded {i} sequences from: {args.fasta_2}")

    logger.warning(f"Writing {len(seq_records)} to {output_path}")
    SeqIO.write(seq_records, output_path, "fasta")


def write_cluster_seqs(sequences, clusters, args, logger):
    """Gather and write out sequences for each cluster to a multiseq FASTA file
    
    :param sequences: dict, id : seq record
    :param clusters: dict, cluster name : {members ids}
    :params args: CLI args parser
    :param logger: logger object
    """
    for cluster in tqdm(clusters, desc="Writing out clusters"):
        if len(clusters[cluster] < args.min_size):
            continue
        
        cluster_seqs = []
        for member_id in clusters[cluster]:
            try:
                cluster_seqs.append(sequences[member_id])
            except KeyError:
                logger.warning(
                    f"Could not find sequence for {member_id} in {args.fasta}\n"
                    f"Not adding {member_id} seq to cluster sequences"
                )

        cluster_path = args.output / f"{cluster}.fasta"
        with open(cluster_path) as fh:
            SeqIO.write(cluster_seqs, fh, "fasta")


if __name__ == "__main__":
    main()
