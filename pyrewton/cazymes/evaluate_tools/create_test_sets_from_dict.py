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

Uses a local CAZyme database to identify CAZy annotated CAZymes.
"""


import logging
import random
import shutil

from datetime import datetime
from typing import List, Optional

from Bio import Entrez, SeqIO
from saintBioutils.genbank import get_genomes, parse_genomes
from saintBioutils.utilities.file_io import get_paths, make_output_directory
from saintBioutils.utilities.logger import config_logger
from tqdm import tqdm

from pyrewton.cazymes.evaluate_tools.test_sets import (
    align_cazymes_and_noncazymes,
    compile_output_file_path,
    write_out_test_set,
)
from pyrewton.utilities.file_io import io_create_eval_testsets
from pyrewton.utilities.parsers.cmd_parser_create_eval_test_sets import build_parser_dict


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
        parser = build_parser_dict()
        args = parser.parse_args()
    else:
        parser = build_parser_dict(argv)
        args = parser.parse_args()

    if logger is None:
        config_logger(args)
    logger = logging.getLogger(__name__)

    Entrez.email = args.email

    time_stamp = datetime.now().strftime("%Y_%m_%d-%H_%M_%S")

    extract_seq_dir = args.output / f"extracted_protein_seqs_{time_stamp}"
    alignment_score_dir = args.output / f"alignment_scores_{time_stamp}"
    test_set_dir = args.output / f"test_sets_{time_stamp}"
    make_output_directory(extract_seq_dir, args.force, args.nodelete)
    make_output_directory(alignment_score_dir, args.force, args.nodelete)
    make_output_directory(test_set_dir, args.force, args.nodelete)
    if args.genomes is None:
        genome_dir = args.output / f"genomes_{time_stamp}"
        make_output_directory(genome_dir, args.force, args.nodelete)

    # get the YAML file containing the genomic assemblies to be used for creating test sets
    assembly_dict = io_create_eval_testsets.retrieve_assemblies_dict(args.yaml)

    # get dict containing the genomic assemblies of all CAZymes in CAZy
    cazy_dict = io_create_eval_testsets.get_cazy_dict(args.cazy)

    header = (
        "Genomic_accession\t"
        "Total_proteins\t"
        "Total_CAZymes\t"
        "Genome_CAZome_percentage\t"
        "CAZome_coverage_percentage\t"
        "CAZyme_sample_size\n"
    )
    coverage_log_path = args.output / f"cazome_coverage_{time_stamp}.txt"
    with open(coverage_log_path, 'a') as fh:
        fh.write(header)
    
    headers = (
        "Genomic_accession\t"
        "Protein_accession\t"
        "CAZyme_classification\n"
    )
    composition_log_path = args.output / f"test_set_composition_{time_stamp}.txt"
    with open(composition_log_path, 'a') as fh:
        fh.write(headers)

    temp_alignment_dir = args.output / "temp_alignment_dir"

    if args.genomes is not None:
        logger.info(f"Retrieving genomes from {args.genomes}")
        genomic_assembly_paths = get_paths.get_file_paths(args.genomes, suffixes=[".gbff.gz"])
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
                assembly_path = get_genomes.get_genomic_assembly(assembly, outdir=genome_dir)
                
            # create a FASTA file containing all proteins sequences in the genomic assembly
            fasta_path = parse_genomes.extract_protein_seqs(
                assembly_path,
                assembly[0],
                txid,
                extract_seq_dir,
            )

            # differentiate between CAZymes and non-CAZymes and get test set of 100 known CAZymes
            (
                selected_cazymes,
                cazyme_fasta,
                non_cazymes,
                noncazyme_fasta,
                total_proteins,
                total_cazymes,
            ) = separate_cazymes_and_noncazymes(
                cazy_dict,
                fasta_path,
                temp_alignment_dir,
                assembly[0],
                args,
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

            final_fasta = compile_output_file_path(fasta_path, args)

            write_out_test_set(
                selected_cazymes,
                non_cazymes,
                alignment_df,
                final_fasta,
                assembly[0],
                total_proteins,
                total_cazymes,
                coverage_log_path,
                composition_log_path,
                args,
            )

    # delete the temporary alignment dir
    shutil.rmtree(temp_alignment_dir)


def separate_cazymes_and_noncazymes(cazy_dict, input_fasta, temp_alignment_dir, assembly, args):
    """Identify CAZymes and non-CAZymes in the input FASTA file.

    Writes non-CAZymes to a FASTA file and builds a non-CAZyme database needed for later alignments
    to identify the non-CAZymes that are most similar to the CAZymes (positive controls)
    selected for the test set.

    :param session: open SQL database session
    :param input_fasta: path to input FASTA file containing all protein seqs from a genomic assembly
    :param temp_alignment_dir: path to dir where alignment inputs/outputs are written
    :param assmebly: genomic assembly accessions retrieved from input yaml file

    Return list of CAZymes, path to temp FASTA of selected CAZymes, dict of non-CAZymes
    {protein_accession: Protein_instance}, path to temporary non-CAZyme dir, and total number
    of proteins and CAZymes extracted from the genome.
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

    protein_count = 0

    with open(noncazyme_fasta, "a") as fh:
        for seq_record in tqdm(
            SeqIO.parse(input_fasta, "fasta"),
            total=len(list(SeqIO.parse(input_fasta, "fasta"))),
            desc="Differentiating CAZymes and non-CAZymes",
        ):
            protein_count += 1
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
    
    total_proteins = len(cazymes) + len(list(non_cazymes.keys()))
    total_cazymes = len(cazymes)

    # randomly selected n CAZymes
    try:
        selected_cazymes = random.sample(cazymes, args.sample_size)
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

    return selected_cazymes, cazyme_fasta, non_cazymes, noncazyme_fasta, protein_count, len(cazymes)


if __name__ == "__main__":
    main()
