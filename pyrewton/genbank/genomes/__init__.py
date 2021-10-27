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
"""Module for handle parsing genomes and genomic assemblies"""


import gzip
import logging

from Bio import SeqIO


def extract_protein_seqs(assembly_path, accession, txid, filestem="genbank_proteins"):
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
                        protein_id = get_record_feature(feature, "protein_id", accession)
                        locus_tag = get_record_feature(feature, "locus_tag", accession)
                        # extract protein sequence
                        seq = get_record_feature(feature, "translation", accession)
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


def get_record_feature(feature, qualifier, accession):
    """Retrieve data from BioPython feature object.

    :param feature: BioPython feature object representing the curernt working protein
    :param qualifier: str, data to be retrieved
    :accession: str, accession of the protein being parsed.

    Return feature data.
    """
    logger = logging.getLogger(__name__)
    try:
        data = feature.qualifiers[qualifier][0]
        return data
    except KeyError:
        logger.warning(
            f"Failed to retrieve feature {qualifier}, returning None value, accession: {accession}"
        )
        return None

