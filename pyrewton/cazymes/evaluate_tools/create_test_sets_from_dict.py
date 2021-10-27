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
"""Retrieves datasets for evaluating the CAZyme prediction tools: dbCAN, CUPP and eCAMI.

Writes out a FASTA per candidate species, containing all the protein sequences to analysed by the
CAZymes prediction tools. The datasets contain an equal number of CAZymes to non-CAZymes.
"""

import logging
import random
import shutil

import pandas as pd
from typing import List, Optional

from Bio import Entrez, SeqIO
from Bio.Blast.Applications import NcbiblastpCommandline, NcbimakeblastdbCommandline
from tqdm import tqdm

from pyrewton.genbank.genomes import download_genomes, parse_genomes
from pyrewton.utilities import config_logger
from pyrewton.utilities.file_io import io_get_paths, make_output_directory, io_create_eval_testsets
from pyrewton.utilities.parsers.cmd_parser_get_evaluation_dataset_from_dict import build_parser


class Protein:
    """A single protein from a genomic assembly.
    
    :attr accession: str, unique GenBank accession
    :attr sequence: str, protein amino acid sequence
    :attr cazyme_classification: int, CAZyme = 1, non-CAZyme = 0.
    """

    def __init__(self, accession, sequence, cazyme_classification):
        self.accession = accession
        self.sequence = sequence
        self.cazyme_classification = cazyme_classification

    def __str__(self):
        return f"<{self.accession} <CAZyme classification = {self.cazyme_classification}>>"

    def __repr__(self):
        return f"<{self.accession} <CAZyme classification = {self.cazyme_classification}>>"


def main(argv: Optional[List[str]] = None, logger: Optional[logging.Logger] = None):
    if argv is None:
        parser = build_parser()
        args = parser.parse_args()
    else:
        parser = build_parser(argv)
        args = parser.parse_args()

    if logger is None:
        config_logger(args)
    logger = logging.getLogger(__name__)

    Entrez.email = args.email

    extract_seq_dir = args.output / "extract_protein_seqs"
    test_set_dir = args.output / "test_sets"
    make_output_directory(extract_seq_dir, args.force, args.nodelete)
    make_output_directory(test_set_dir, args.force, args.nodelete)

    # get the YAML file containing the genomic assemblies to be used for creating test sets
    assembly_dict = io_create_eval_testsets.retrieve_assemblies_dict(args.yaml)

    # get dict containing the genomic assemblies of all CAZymes in CAZy
    cazy_dict = io_create_eval_testsets.get_cazy_dict(args.cazy)

    temp_alignment_dir = args.output / "temp_alignment_dir"

    if args.genomes is not None:
        logger.info(f"Retrieving genomes from {args.genomes}")
        genomic_assembly_paths = io_get_paths.get_genomic_assembly_paths(args)
        logger.info(f"Retrieved {len(genomic_assembly_paths)} genomic assemblies")
        genomic_path_dict = {}
        for _path in genomic_assembly_paths:
            filename = (_path.name).split(".")[0]
            filename = filename.split("_")
            accession = f"{filename[0]}_{filename[1]}.{filename[2]}"
            genomic_path_dict[accession] = _path

    # create a test set for each genomic assembly
    for txid in tqdm(assembly_dict, desc="Parsing assemblies in config file"):
        for assembly in assembly_dict[txid]:
            # assembly = [genomic accession, assembly name]

            # whipe temp dir clean
            io_create_eval_testsets.prepare_output_dir(temp_alignment_dir)

            if args.genomes is not None:
                assembly_path = genomic_path_dict[assembly[0]]
            else:
                # download genomic assembly
                assembly_path = download_genomes.get_genomic_assembly(assembly)
                
            # create a FASTA file containing all proteins sequences in the genomic assembly
            fasta_path = parse_genomes.extract_protein_seqs(
                assembly_path,
                assembly[0],
                txid,
                args.output,
                subfolder="extract_protein_seqs",
            )

            # differentiate between CAZymes and non-CAZymes and get test set of 100 known CAZymes
            (
                selected_cazymes,
                cazyme_fasta,
                non_cazymes,
                noncazyme_fasta,
            ) = separate_cazymes_and_noncazymes(
                cazy_dict,
                fasta_path,
                temp_alignment_dir,
                assembly[0],
            )

            if selected_cazymes is None:
                continue

            alignment_df = align_cazymes_and_noncazymes(
                cazyme_fasta,
                noncazyme_fasta,
                temp_alignment_dir,
            )
            if alignment_df is None:
                continue

            final_fasta = compile_output_file_path(fasta_path)

            write_out_test_set(
                selected_cazymes,
                non_cazymes,
                alignment_df,
                final_fasta,
                assembly[0],
            )

    # delete the temporary alignment dir
    shutil.rmtree(temp_alignment_dir)


def separate_cazymes_and_noncazymes(cazy_dict, input_fasta, temp_alignment_dir, assembly):
    """Identify CAZymes and non-CAZymes in the input FASTA file.

    Writes non-CAZymes to a FASTA file and builds a non-CAZyme database needed for later alignments
    to identify the non-CAZymes that are most similar to the CAZymes (positive controls)
    selected for the test set.

    :param session: open SQL database session
    :param input_fasta: path to input FASTA file containing all protein seqs from a genomic assembly
    :param temp_alignment_dir: path to dir where alignment inputs/outputs are written
    :param assmebly: genomic assembly accessions retrieved from input yaml file

    Return list of CAZymes, path to temp FASTA of selected CAZymes, dict of non-CAZymes
    {protein_accession: Protein_instance}, and path to temporary non-CAZyme dir.
    """
    logger = logging.getLogger(__name__)
    cazyme_accessions = set()  # ott checking duplicate cazymes not selected
    cazymes = []
    non_cazymes = {}  # {accession: Protein_instance}

    # create temporary, empty dir for holding alignment input and output
    io_create_eval_testsets.prepare_output_dir(temp_alignment_dir)

    # alignment_output = temp_alignment_dir / "alignment_output.tab"
    noncazyme_fasta = temp_alignment_dir / "noncazymes.fasta"
    cazyme_fasta = temp_alignment_dir / "cazymes.fasta"

    with open(noncazyme_fasta, "a") as fh:
        for seq_record in tqdm(
            SeqIO.parse(input_fasta, "fasta"),
            total=len(list(SeqIO.parse(input_fasta, "fasta"))),
            desc="Differentiating CAZymes and non-CAZymes",
        ):
            accession = seq_record.id
            accession = accession.split(" ")[0]
            sequence = seq_record.seq

            # get the CAZyme classification, 1 = CAZyme, 0 = non-CAZyme
            try:
                cazy_dict[accession]
                cazyme_classification = 1
                if accession not in cazyme_accessions:
                    cazyme_accessions.add(accession)
                    cazymes.append(Protein(accession, sequence, cazyme_classification))

            except KeyError:
                cazyme_classification = 0
                non_cazymes[accession] = Protein(accession, sequence, cazyme_classification)
                # write to FASTA file
                seq = str(sequence)
                seq = "\n".join([seq[i : i + 60] for i in range(0, len(seq), 60)])
                file_content = f">{accession} \n{seq}\n"
                fh.write(file_content)

    # randomly selected 100 CAZymes
    try:
        selected_cazymes = random.sample(cazymes, 100)
        # write selected CAZymes to FASTA file in temp dir
        with open(cazyme_fasta, "w") as fh:
            for cazyme in tqdm(selected_cazymes, desc="Writing selected CAZymes to FASTA"):
                # write query CAZyme to FASTA and build a database
                seq = str(cazyme.sequence)
                seq = "\n".join([seq[i : i + 60] for i in range(0, len(seq), 60)])

                file_content = f">{cazyme.accession} \n{seq}\n"

                fh.write(file_content)

    except ValueError as e:
        logger.warning(f"No CAZymes selected for {assembly}, the following error was raised:\n{e}")
        selected_cazymes = None

    return selected_cazymes, cazyme_fasta, non_cazymes, noncazyme_fasta


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

    :param fasta_path: path to fasta file containing protein sequences extracted from genome
    :param args: cmd-line args parser

    Return path to output FASTA file.
    """
    genomic_assembly = fasta_path.name
    fasta = genomic_assembly.replace(".fasta", "_test_set.fasta")
    fasta = args.output / "test_sets" / fasta

    return fasta


def write_out_test_set(selected_cazymes, non_cazymes, alignment_df, final_fasta, genomic_acc):
    """Create the test set and write out a FASTA file.

    The FASTA file contains the protein sequences of the test set (in a random order).

    :param selected_cazymes: list of Protein instances of CAZymes for the test set
    :param non_cazymes: dict of all non-CAZymes from the input FASTA file {accession: Protein instance}
    :param alignment_df: Pandas df of Blast Score Ratios of non-CAZyme alignemnts against selected CAZymes
    :param final_fasta: path to output FASTA file of the final test set
    :param genomic_acc: str, genomic assembly accession

    Return nothing.
    """
    logger = logging.getLogger(__name__)

    # sort results by BLAST score ratio
    alignment_results = alignment_df.sort_values(by=["blast_score_ratio"], ascending=False)

    # write alignment scores to csv for future reference and documentation
    csv_fh = final_fasta.replace("_test_set.fasta", "_alignment_scores.csv")
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

    for protein in all_selected_proteins:
        seq = str(protein.sequence)
        seq = "\n".join([seq[i : i + 60] for i in range(0, len(seq), 60)])
        file_content = f">{protein.accession}\n{seq}\n"
        with open(final_fasta, "a") as fh:
            fh.write(file_content)

    return


if __name__ == "__main__":
    main()
