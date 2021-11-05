#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) University of St Andrews 2020-2021
# (c) University of Strathclyde 2020-2021
# (c) James Hutton Institute 2020-2021
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
"""Functions for compiling and writing out test sets"""


import logging
import random

import pandas as pd

from pathlib import Path

from Bio.Blast.Applications import NcbiblastpCommandline
from tqdm import tqdm


def align_cazymes_and_noncazymes(cazyme_fasta, noncazyme_fasta, temp_alignment_dir):
    """Coordinate alignment of each selected CAZyme against the non-CAZymes.

    :param cazyme_fasta: Path to FASTA containing selected CAZymes from the proteome
    :param noncazyme_fasta: Path to FASTA containing all non-CAZymes from the proteome
    :param temp_alignment_dir: Path to temporary dir used for storing alignment i/o

    Return pandas df of alignment results.
    """
    logger = logging.getLogger(__name__)
    logger.warning("Starting alignment")
    alignment_output = temp_alignment_dir / "alignment_output.tab"
    # perform all-verus-all alignment
    all_v_all_blastp = NcbiblastpCommandline(
        query=noncazyme_fasta,
        subject=cazyme_fasta,
        out=alignment_output,
        outfmt="6 qseqid sseqid pident qcovs qlen slen length bitscore evalue",
        threshold=100,
    )
    stdout, stderr = all_v_all_blastp()

    # check if alignment was successful
    if len(stderr) != 0:
        logger.warning(
            f"Error raised by performing all-vs-all BLASTP\nstdout={stdout}\nstderr={stderr}"
        )
        return

    alignment_results = pd.read_csv(alignment_output, sep="\t", header=None)
    headers = ["query (non-CAZyme)", "subject (Cazyme)", "identity", "coverage",
               "qlength", "slength", "alength",
               "bitscore", "E-value"]
    alignment_results.columns = headers

    # Create a new column: blast_score_ratio a.k.a normalised bitscore
    alignment_results['blast_score_ratio'] = alignment_results.bitscore/alignment_results.qlength
    logger.warning("Alignment complete and Blast score ratio calculated")
    return alignment_results


def compile_output_file_path(fasta_path, args):
    """Compile paths write out the final test sets (FASTA files).

    :param fasta_path: str, path to fasta file containing protein sequences extracted from genome
    :param args: cmd-line args parser

    Return path to output FASTA file.
    """
    fasta_path = Path(fasta_path)
    genomic_assembly = fasta_path.name
    fasta = genomic_assembly.replace(".fasta", "_test_set.fasta")
    fasta = args.output / "test_sets" / fasta

    return fasta


def write_out_test_set(
    selected_cazymes,
    non_cazymes,
    alignment_df,
    final_fasta,
    genomic_acc,
    total_proteins,
    total_cazymes,
    coverage_log_path,
    composition_log_path,
    args,
):
    """Create the test set and write out a FASTA file.

    The FASTA file contains the protein sequences of the test set (in a random order).

    :param selected_cazymes: list of Protein instances of CAZymes for the test set
    :param non_cazymes: dict of all non-CAZymes from the input FASTA file {accession: Protein instance}
    :param alignment_df: Pandas df of Blast Score Ratios of non-CAZyme alignemnts against selected CAZymes
    :param final_fasta: path to output FASTA file of the final test set
    :param genomic_acc: str, genomic assembly accession
    :param total_proteins: int, total number of proteins extracted from the genomic assembly
    :param total_proteins: int, total number of CAZy annotated CAZymes extracted from the genomic assembly
    :param coverage_log_path: Path, path for logging test set genome coverage stats
    :param composition_log_path: Path, path for logging the composition of the test sets
    :param args: cmd-line args parser

    Return nothing.
    """
    logger = logging.getLogger(__name__)

    # sort results by BLAST score ratio
    alignment_results = alignment_df.sort_values(by=["blast_score_ratio"], ascending=False)

    # write alignment scores to csv for future reference and documentation
    csv_filename = (final_fasta.name).replace('_test_set.fasta', "_alignment_scores.csv")
    csv_fh = args.output / "alignment_scores" / csv_filename
    alignment_results.to_csv(csv_fh)

    selected_noncazyme_accessions = set()  # used to prevent introduction of duplicates

    # create empty df to store selected non-CAZymes
    headers = list(alignment_results.columns)
    selected_non_cazymes = pd.DataFrame(columns=headers)

    for row_index in tqdm(
        range(len(alignment_results["blast_score_ratio"])),
        desc="Selecting unique non-CAZymes",
    ):
        df_row = alignment_results.iloc[row_index]
        accession = df_row['query (non-CAZyme)']
        if accession not in selected_noncazyme_accessions:
            selected_noncazyme_accessions.add(accession)
            selected_non_cazymes = selected_non_cazymes.append(df_row)
        if len(selected_non_cazymes["blast_score_ratio"]) == 200:
            break
    
    if len(selected_non_cazymes["blast_score_ratio"]) != 200:
        noncazyme_count = len(selected_non_cazymes["blast_score_ratio"])
        logger.error(
            f"Could not retrieve 100 non-CAZymes from {genomic_acc}\n"
            f"Only retrieved {noncazyme_count} non-CAZymes\n"
            f"Test set for {genomic_acc} not created."
        )
        return

    all_selected_proteins = []  # list of Protein instances

    index = 0
    for index in tqdm(
        range(len(selected_non_cazymes["blast_score_ratio"])),
        desc="Retrieving selected non-CAZymes",
    ):
        df_row = alignment_results.iloc[index]
        protein_accession = df_row['query (non-CAZyme)']
        all_selected_proteins.append(non_cazymes[protein_accession])  # retrieve Protein instance from dict

    all_selected_proteins += selected_cazymes
    random.shuffle(all_selected_proteins)  # write out the test set in a random order

    composition_log_content = []

    for protein in all_selected_proteins:
        seq = str(protein.sequence)
        seq = "\n".join([seq[i : i + 60] for i in range(0, len(seq), 60)])
        file_content = f">{protein.accession}\n{seq}\n"
        with open(final_fasta, "a") as fh:
            fh.write(file_content)
    
    with open(composition_log_path, 'a') as fh:
        for protein in selected_cazymes:
            new_line = f"{genomic_acc}\t{protein.accession}\t1\n"
            fh.write(new_line)

        for accession in selected_noncazyme_accessions:
            new_line = f"{genomic_acc}\t{accession}\t0\n"
            fh.write(new_line)

    cazome_genome_coverage = ((total_cazymes / total_proteins) * 100)
    cazome_coverage = ((args.sample_size / total_cazymes) * 100)

    log_message = (
        f"{genomic_acc}\t"
        f"{total_proteins}\t"
        f"{total_cazymes}\t"
        f"{cazome_genome_coverage}\t"
        f"{cazome_coverage}\t"
        f"{args.sample_size}\n"
    )

    with open(coverage_log_path, 'a') as fh:
        fh.write(log_message)

    return
