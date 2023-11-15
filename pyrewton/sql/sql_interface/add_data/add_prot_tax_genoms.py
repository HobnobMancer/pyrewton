#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) University of St Andrews 2020-2021
# (c) University of Strathclyde 2020-2021
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
#
"""Script housing funcs for adding data from GenBank assemblies and proteins to the db"""


import logging
import sys

from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from saintBioutils.utilities.file_io import get_paths
from tqdm import tqdm

from pyrewton.sql.sql_interface.add_data import insert_data
from pyrewton.sql.sql_interface.load_data import load_genbank_data


def add_protein_data(connection, args, logger):
    """Load protein data from FASTA files and local CAZome database
    and add new protein, taxonomy and genomic assembly data to the local CAZome db.
    
    :param connection: open connection to sqlite3 db
    :param args: CLI args parser
    :param logger: logger object
    """
    # get proteins from FASTA files of proteins to be added to the local CAZome database
    protein_fasta_files = get_paths.get_file_paths(args.protein_dir, suffixes='fasta')

    if len(protein_fasta_files) == 0:
        logger.error(
            f"No FASTA files containing protein sequences extracted from the genomes retrieved from\n"
            f"{args.protein_dir}\nCheck the correct dir was provided\n"
            "Terminating program"
        )
        sys.exit(1)

    # retrieve protein data from the FASTA file parsed by dbCAN
    # tax_dict = {genomic_accession: {'genus': str, 'species': str, 'tax_id': str}}
    # protein_dict = {genomic_accession: {protein_accession: str(sequence)}}
    protein_dict, tax_dict = parse_protein_files(protein_fasta_files)

    # check if need to add any new taxs
    # test taxs using NCBI tax id, as genus and species name could change for the same NCBI tax id
    # complete_tax_table = {ncbi tax id: {'genus': genus, 'species': species, 'db_id': db_id}
    complete_tax_table = load_genbank_data.get_complete_tax_table(connection)
    taxs_to_add = []
    for genomic_acc in tax_dict:
        if tax_dict[genomic_acc]['tax_id'] not in complete_tax_table.keys():
            taxs_to_add.append( (
                tax_dict[genomic_acc]['tax_id'],
                tax_dict[genomic_acc]['genus'],
                tax_dict[genomic_acc]['species'],
            ) )
    if len(taxs_to_add) != 0:
        logger.warning(f"Adding {len(taxs_to_add)} new taxes to the Taxonomies table")
        insert_data(connection, "Taxonomies", ['ncbi_tax_id', 'genus', 'species'], taxs_to_add)

    # check if need to add new genomes
    # {genomic_accession: {taxonomy_id: int, db_id: int}}
    complete_assemblies_table = load_genbank_data.get_complete_assemblies_table(connection)
    assemblies_to_add = []
    for genomic_acc in tax_dict:
        if genomic_acc in complete_assemblies_table.keys():
            ncbi_tax_id = complete_tax_table[genomic_acc]['tax_id']
            assemblies_to_add.append( (ncbi_tax_id, genomic_acc) )
    if len(assemblies_to_add) != 0:
        logger.warning(f"Adding {len(assemblies_to_add)} new genomic assemblies to the Assemblies table")
        insert_data(connection, "Assemblies", ['taxonomy_id', 'assembly_accession'], assemblies_to_add)

    # test if need to add any new proteins
    existing_db_proteins = load_genbank_data.get_protein_db_ids(connection)
    # protein_dict = {genomic_accession: {protein_accession: str(sequence)}}
    proteins_to_add = []
    for genomic_acc in protein_dict:
        for protein_acc in protein_dict[genomic_acc]:
            if protein_acc not in existing_db_proteins.keys():
                assembly_id = complete_assemblies_table[genomic_acc]['db_id']

                sequence = protein_dict[genomic_acc][protein_acc]['sequence']
                processed_sequence = sequence.replace("X", "G")
                processed_sequence = processed_sequence.replace("Z", "G")
                processed_sequence = processed_sequence.replace("B", "G")
                processed_sequence = processed_sequence.replace("J", "G")
                analysed_seq = ProteinAnalysis(processed_sequence)
                mass = analysed_seq.molecular_weight()
                length = len(sequence)

                proteins_to_add.append( (
                    assembly_id,
                    protein_acc,
                    mass,
                    length,
                    processed_sequence,
                ) )
    if len(assemblies_to_add) != 0:
        logger.warning(f"Adding {len(proteins_to_add)} new proteins to the Proteins table")
        insert_data(connection, "Proteins", ['assembly_id', 'genbank_accession', 'mass', 'length', 'sequence'], assemblies_to_add)


def parse_protein_files(protein_fasta_files):
    """Retrieve protein, genomic and taxonomic data from the FASTA files parsed by dbCAN.
    
    :param protein_fasta_files: list of Paths, one path per FASTA file
    
    Return two dicts
    tax_dict = {genomic_accession: {'genus': str, 'species': str, 'tax_id': str}}
    protein_dict = {genomic_accession: {protein_accession: str(sequence)}}
    """
    logger = logging.getLogger(__name__)

    protein_dict = {}
    tax_dict = {}  # {genomes: {'genus': str, 'species': str, 'tax_id': str}}

    for fasta in tqdm(protein_fasta_files, desc="Retrieving protein data from FASTA files"):
        for record in SeqIO.parse(fasta, 'fasta'):
            protein_acc = record.id

            desc = record.description.split(" ")

            assembly_acc = desc[2]
            assembly_acc = assembly_acc.replace("_", ".")
            assembly_acc = assembly_acc.replace("GCA.", "GCA_")
            assembly_acc = assembly_acc.replace("GCF.", "GCF_")

            tax_id = desc[3].replace("txid", "")
            
            genus = desc[4]
            species = desc[5]

            # add taxonomy data
            try:
                tax_dict[assembly_acc]
            
            except KeyError:
                tax_dict[assembly_acc] = {
                    'genus': genus,
                    'species': species,
                    'txid': tax_id,
                }

            # add protein data
            try:
                protein_dict[assembly_acc]

                try:
                    protein_dict[assembly_acc][protein_acc]

                    if protein_dict[assembly_acc][protein_acc] != str(record.seq):
                        logger.warning(
                            f"Multiple unique sequences retrieved for {protein_acc}\n"
                            "Retrieving the first protein sequence only"
                        )

                except KeyError:
                    protein_dict[assembly_acc][protein_acc] = str(record.seq)

            except KeyError:
                protein_dict[assembly_acc] = {
                    protein_acc: str(record.seq),
                }

    return protein_dict, tax_dict
