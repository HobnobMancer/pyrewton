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
"""Retrieves datasets for evaluating the CAZyme prediction tools: dbCAN, CUPP and eCAMI.

Writes out a FASTA per candidate species, containing all the protein sequences to analysed by the
CAZymes prediction tools. The datasets contain an equal number of CAZymes to non-CAZymes.
"""

import gzip
import json
import logging
import random
import re
import yaml
import sys
import shutil
import time

import pandas as pd
from pathlib import Path
from socket import timeout
from typing import List, Optional
from urllib.error import HTTPError, URLError
from urllib.request import urlopen

from Bio import Entrez, SeqIO
from Bio.Blast.Applications import NcbiblastpCommandline, NcbimakeblastdbCommandline
from tqdm import tqdm

from pyrewton.utilities import config_logger
from pyrewton.utilities.file_io import make_output_directory
from pyrewton.utilities.parsers.cmd_parser_get_evaluation_dataset_from_dict import build_parser


class Protein:
    """A single protein from a genomic assembly."""

    def __init__(self, accession, sequence, cazyme_classification):
        """Initiate class instance.

        :attr accession: str, unique GenBank accession
        :attr sequence: str, protein amino acid sequence
        :attr cazyme_classification: int, CAZyme = 1, non-CAZyme = 0.
        """
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

    make_output_directory(args.output, args.force, args.nodelete)

    # get the YAML file containing the genomic assemblies to be used for creating test sets
    assembly_dict = retrieve_assemblies_dict(args.yaml)

    # get dict containing the genomic assemblies of all CAZymes in CAZy
    cazy_dict = get_cazy_dict(args.cazy)

    temp_alignment_dir = args.output / "temp_alignment_dir"

    # create a test set for each genomic assembly
    for txid in tqdm(assembly_dict, desc="Parsing assemblies in config file"):
        for assembly in assembly_dict[txid]:
            # whipe temp dir clean
            prepare_output_dir(temp_alignment_dir)

            # download genomic assembly
            assembly_path = get_genomic_assembly(assembly, txid)

            # create a FASTA file containing all proteins sequences in the genomic assembly
            fasta_path = get_protein_seqs(assembly_path, assembly, txid)

            # differentiate between CAZymes and non-CAZymes and get test set of 100 CAZymes
            selected_cazymes, cazyme_fasta, non_cazymes, noncazyme_fasta = differentiate_cazymes_and_noncazymes(
                cazy_dict,
                fasta_path,
                temp_alignment_dir,
            )

            if selected_cazymes is None:
                continue

            alignment_df = align_cazymes_and_noncazymes(cazyme_fasta, noncazyme_fasta, temp_alignment_dir)
            if alignment_df is None:
                continue

            final_fasta = compile_output_file_path(fasta_path)

            write_out_test_set(selected_cazymes, non_cazymes, alignment_df, final_fasta)

    # delete the temporary alignment dir
    shutil.rmtree(temp_alignment_dir)


def retrieve_assemblies_dict(yaml_path):
    """Retrieve dict of genomic assemblies to create test sets from.

    :param yaml_path: Path to yaml file, keyed by NCBI:txid and valued by genomic assemblies

    Return dict, keyed by NCBI:txid and valued by list genomic assemblies"""
    logger = logging.getLogger(__name__)
    # open the YAML file
    try:
        with open(yaml_path, "r") as fh:
            assemblies_dict = yaml.full_load(fh)
    except FileNotFoundError:
        logger.error(
            "Did not find the configuration file. Check the path is correct.\n"
            "Terminating programme"
        )
        sys.exit(1)
    return assemblies_dict


def get_cazy_dict(cazy_path):
    logger = logging.getLogger(__name__)
    try:
        with open(cazy_path, "r") as fh:
            cazy_dict = json.load(fh)
    except FileNotFoundError:
        logger.error(
            "Did not find the configuration file. Check the path is correct.\n"
            "Terminating programme"
        )
        sys.exit(1)
    return cazy_dict


def get_genomic_assembly(assembly_accession, txid, suffix="genomic.gbff.gz"):
    """Coordinate downloading Genomic assemmbly from the NCBI Assembly database.

    :param assembly_accession: str, accession of the Genomic assembly to be downloaded
    :param txid: str, NCBI Taxonomy ID of the host species of the genomic assembly
    :param suffix: str, suffix of file

    Return path to downloaded genomic assembly.
    """
    # compile url for download
    genbank_url, filestem = compile_url(assembly_accession, suffix)

    # create path to write downloaded Genomic assembly to
    out_file_path = "_".join([filestem.replace(".", "_"), suffix])
    out_file_path = Path(out_file_path)

    # download GenBank file
    download_file(genbank_url, out_file_path, assembly_accession, "GenBank file")

    return out_file_path


def compile_url(accession_number, suffix, ftpstem="ftp://ftp.ncbi.nlm.nih.gov/genomes/all"):
    """Retrieve URL for downloading the assembly from NCBI, and create filestem of output file path

    :param accession_number: str, asseccion number of genomic assembly

    Return str, url required for download and filestem for output file path for the downloaded assembly.
    """
    # search for the ID of the record
    with entrez_retry(
        10, Entrez.esearch, db="Assembly", term="GCA_000021645.1[Assembly Accession]", rettype='uilist',
    ) as handle:
        search_record = Entrez.read(handle)

    # retrieve record for genomic assembly
    with entrez_retry(10, Entrez.esummary, db="assembly", id=search_record['IdList'][0], report="full") as handle:
        record = Entrez.read(handle)

    assembly_name = record["DocumentSummarySet"]["DocumentSummary"][0]["AssemblyName"]

    escape_characters = re.compile(r"[\s/,#\(\)]")
    escape_name = re.sub(escape_characters, "_", assembly_name)

    # compile filstem
    filestem = "_".join([accession_number, escape_name])

    # separate out filesteam into GCstem, accession number intergers and discarded
    url_parts = tuple(filestem.split("_", 2))

    # separate identifying numbers from version number
    sub_directories = "/".join(
        [url_parts[1][i : i + 3] for i in range(0, len(url_parts[1].split(".")[0]), 3)]
    )

    # return url for downloading file
    return (
        "{0}/{1}/{2}/{3}/{3}_{4}".format(
            ftpstem, url_parts[0], sub_directories, filestem, suffix
        ),
        filestem,
    )


def entrez_retry(retries, entrez_func, *func_args, **func_kwargs):
    """Call to NCBI using Entrez.

    :param retries: int, maximum number of retries excepted if network error encountered
    :param entrez_func: function, call method to NCBI
    :param *func_args: tuple, arguments passed to Entrez function
    :param ** func_kwargs: dictionary, keyword arguments passed to Entrez function

    Returns record.
    """
    logger = logging.getLogger(__name__)
    record, retries, tries = None, retries, 0

    while record is None and tries < retries:
        try:
            record = entrez_func(*func_args, **func_kwargs)

        except IOError:
            # log retry attempt
            if tries < retries:
                logger.warning(
                    f"Network error encountered during try no.{tries}.\nRetrying in 10s",
                    exc_info=1,
                )
                time.sleep(10)
            tries += 1

    if record is None:
        logger.error(
            "Network error encountered too many times. Exiting attempt to call to NCBI"
        )
        return

    return record


def download_file(
    genbank_url, out_file_path, accession_number, file_type
):
    """Download file.

    :param genbank_url: str, url of file to be downloaded
    :param out_file_path: path, output directory for file to be written to
    :param accession_number: str, accession number of genome
    :param file_type: str, denotes in logger file type downloaded

    Return nothing.
    """
    logger = logging.getLogger(__name__)
    # Try URL connection
    try:
        response = urlopen(genbank_url, timeout=45)
    except (HTTPError, URLError, timeout) as e:
        logger.error(
            f"Failed to download {file_type} for {accession_number}", exc_info=1,
        )
        return

    if out_file_path.exists():
        logger.warning(f"Output file {out_file_path} exists, not downloading")
        return

    # Download file
    file_size = int(response.info().get("Content-length"))
    bsize = 1_048_576
    try:
        with open(out_file_path, "wb") as out_handle:
            # Using leave=False as this will be an internally-nested progress bar
            with tqdm(
                total=file_size,
                leave=False,
                desc=f"Downloading {accession_number} {file_type}",
            ) as pbar:
                while True:
                    buffer = response.read(bsize)
                    if not buffer:
                        break
                    pbar.update(len(buffer))
                    out_handle.write(buffer)
    except IOError:
        logger.error(f"Download failed for {accession_number}", exc_info=1)
        return

    return


def get_protein_seqs(assembly_path, accession, txid, filestem="genbank_proteins"):
    """Retrieve annoated protein sequences from genomic assembly and write to a single FASTA file.

    :param assemly_path: Path to genomic assembly
    :param accession: str, accession number of the genomic assembly
    :param txid: str, NCBI taxonomy id of the host species
    :param filestem: str, file name prefix

    Return path to FASTA file containing the protein sequences from the assembly.
    """
    logger = logging.getLogger(__name__)

    # build path to the output FASTA file
    fasta_path = f"{filestem}_{txid}_{accession}.fasta"

    protein_count = 0

    with open(fasta_path, "a") as fh:
        with gzip.open(assembly_path, "rt") as handle:  # unzip the genomic assembly
            # parse proteins in the genomic assembly
            for gb_record in SeqIO.parse(handle, "genbank"):
                for (index, feature) in enumerate(gb_record.features):
                    # Parse over only protein encoding features (type = 'CDS')
                    if feature.type == "CDS":
                        # retrieve data from protein feature record
                        protein_id = get_record_feature(feature, "protein_id", accession_number)
                        locus_tag = get_record_feature(feature, "locus_tag", accession_number)
                        # extract protein sequence
                        seq = get_record_feature(feature, "translation", accession_number)
                        if seq is None:
                            continue

                        # create file content for writing protein to fasta file
                        # FASTA sequences have 60 characters per line
                        seq = "\n".join([seq[i : i + 60] for i in range(0, len(seq), 60)])
                        protein_id = protein_id + " " + locus_tag

                        file_content = f">{protein_id} \n{seq}\n"

                        fh.write(file_content)

                        protein_count += 1

    logger.warning(f"{protein_count} proteins in genomic assembly {accession}")

    return fasta_path


def get_record_feature(feature, qualifier, accession_number):
    """Retrieve data from BioPython feature object.

    :param feature: BioPython feature object representing the curernt working protein
    :param qualifier: str, data to be retrieved
    :accession_number: str, accession of the protein being parsed.

    Return feature data.
    """
    logger = logging.getLogger(__name__)
    try:
        data = feature.qualifiers[qualifier][0]
        return data
    except KeyError:
        logger.warning(f"Failed to retrieve feature {qualifier}, returning None value, accession: {accession}")
        return None


def prepare_output_dir(output_dir):
    """Delete and make temporary output directory"""
    output_dir = Path(output_dir)
    try:
        if output_dir.exists:
            shutil.rmtree(output_dir)
    except (FileExistsError, FileNotFoundError) as e:
        pass
    try:
        output_dir.mkdir(exist_ok=True)
    except FileExistsError:
        print("output_dir_exists")
    return


def differentiate_cazymes_and_noncazymes(cazy_dict, input_fasta, temp_alignment_dir, assembly):
    """Identify CAZymes and non-CAZymes in the input FASTA file.

    Writes non-CAZymes to a FASTA fie and builds a non-CAZyme database needed for later alignments
    to identify the non-CAZymes that are most similar to the CAZymes (positive controls)
    selected for the test set.

    :param session: open SQL database session
    :param input_fasta: path to input FASTA file containing all protein seqs from a genomic assembly

    Return list of CAZymes, path to temp FASTA of selected CAZymes, dict of non-CAZymes
    {protein_accession: Protein_instance}, and path to temporary non-CAZyme dir.
    """
    logger = logging.getLogger(__name__)
    cazymes = []
    non_cazymes = {}  # {accession: Protein_instance}

    # create temporary, empty dir for holding alignment input and output
    prepare_output_dir(temp_alignment_dir)
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

    return alignment_results


def compile_output_file_path(input_fasta):
    """Compile paths to the FASTA file for the test set from the input FASTA file.

    :param input_fasta: path to the input fasta file

    Return path to output FASTA file.
    """
    input_fasta = input_fasta.split("/")
    genomic_assembly = input_fasta[0]

    fasta = genomic_assembly.replace("fasta", "_test_set.fasta")

    return fasta


def write_out_test_set(selected_cazymes, non_cazymes, alignment_df, final_fasta):
    """Create the test set and write out a FASTA file.

    The FASTA file contains the protein sequences of the test set (in a random order).

    :param selected_cazymes: list of Protein instances of CAZymes for the test set
    :param non_cazymes: dict of all non-CAZymes from the input FASTA file {accession: Protein instance}
    :param alignment_df: Pandas df of Blast Score Ratios of non-CAZyme alignemnts against selected CAZymes
    :param final_fasta: path to output FASTA file of the final test set

    Return nothing.
    """
    # sort results by BLAST score ratio
    alignment_results = alignment_df.sort_values(by=["blast_score_ratio"], ascending=False)

    # write alignment scores to csv for future reference and documentation
    csv_fh = final_fasta.replace("_test_set.fasta", "_alignment_scores.csv")
    alignment_results.to_csv(csv_fh)

    selected_non_cazymes = alignment_results[:200]

    all_selected_proteins = []  # list of Protein instances

    index = 0
    for index in tqdm(range(len(selected_non_cazymes["blast_score_ratio"])), desc="Retrieving selected non-CAZymes"):
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
