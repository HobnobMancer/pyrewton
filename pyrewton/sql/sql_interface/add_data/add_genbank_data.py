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
"""Module for defining the db ORM and interacting with the db"""


import logging
from Bio.SeqUtils.ProtParam import ProteinAnalysis

from tqdm import tqdm

from pyrewton.sql.sql_interface.add_data import insert_data


def add_species_data(tax_dict, connection):
    species_to_insert = set()
    
    for genomic_accession in tqdm(tax_dict, desc="Adding species to db"):
        tax_id = tax_dict[genomic_accession]['txid']
        genus = tax_dict[genomic_accession]['genus']
        species = tax_dict[genomic_accession]['species']

        species_to_insert.add( (tax_id, genus, species) )

    if len(species_to_insert) != 0:
        insert_data(connection, 'Taxonomies', ['ncbi_tax_id', 'genus', 'species'], list(species_to_insert))

    return


def add_genomic_accessions(tax_dict, tax_table_dict, connection):
    """Add genomic accessions to the database.
    
    :param tax_dict: {genomic_accession: {'genus': str, 'species': str, 'tax_id': str}}
    :param tax_table_dict: {genus species: db_id}
    :param connection: open sqlalchemy connection to an SQLite3 engine
    
    Return nothing
    """
    assemblies_to_insert = set()
    for genomic_accession in tqdm(tax_dict, desc="Adding assemblies to db"):
        genus = tax_dict[genomic_accession]['genus']
        species = tax_dict[genomic_accession]['species']

        db_tax_id = tax_table_dict[f"{genus} {species}"]

        assemblies_to_insert.add( (db_tax_id, genomic_accession) )

    if len(assemblies_to_insert) != 0:
        insert_data(connection, 'Assemblies', ['taxonomy_id', 'assembly_accession'], list(assemblies_to_insert))

    return


def add_proteins(protein_dict, assembly_dict, connection):
    """Add proteins to the db
    
    :param protein_dict:
    :param assembly_dict:
    :param connection: open sqlalchemy connection to an SQLite3 engine
    
    Return nothing
    """
    proteins_to_insert = set()

    for genomic_accession in tqdm(protein_dict, desc="Adding proteins to db"):
        protein_accessions = list(protein_dict[genomic_accession].keys())
        
        for protein_accession in protein_accessions:
            assembly_db_id = assembly_dict[genomic_accession]

            sequence = protein_dict[genomic_accession][protein_accession]['sequence']
            processed_sequence = sequence.replace("X", "G")
            processed_sequence = processed_sequence.replace("Z", "G")
            processed_sequence = processed_sequence.replace("B", "G")
            processed_sequence = processed_sequence.replace("J", "G")
            analysed_seq = ProteinAnalysis(processed_sequence)
            mass = analysed_seq.molecular_weight()
            length = len(sequence)

            proteins_to_insert.add( (assembly_db_id, protein_accession, mass, length, sequence) )

    if len(proteins_to_insert) != 0:
        insert_data(
            connection,
            'Proteins',
            ['assembly_id', 'genbank_accession', 'mass', 'length', 'sequence'],
            list(proteins_to_insert),
        )

    return


def add_families(cazome_dict, connection):
    """Add CAZy families to db.
    
    :param cazome_dict:
    :param connection: open sqlalchemy connection to an SQLite3 engine
    
    Return nothing
    """
    families_to_insert = set()

    tools = ['dbcan', 'hmmer', 'hotpep', 'diamond', 'cazy']

    genomic_accessions = list(cazome_dict.keys())
    for genomic_accession in tqdm(genomic_accessions, desc="Adding families to db"):
        protein_accessions = list(cazome_dict[genomic_accession].keys())
        for protein in protein_accessions:
            for tool in tools:
                try:
                    fams = cazome_dict[genomic_accession][protein][tool]
                    if fams is not None:
                        for fam in fams:
                            families_to_insert.add( (fam,) )
                except KeyError:
                    pass  # raised if data not retrieved from CAZy

    if len(families_to_insert) != 0:
        insert_data(
            connection,
            'CazyFamilies',
            ['family'],
            list(families_to_insert),
        )
    
    return


def add_classifications(
    cazome_dict,
    domain_dict,
    protein_db_dict,
    classifer_db_dict,
    family_db_dict,
    connection,
    args,
):
    """Add domain classifications from dbCAN (and CAZy) to the db.
    
    :param cazome_dict:
    {genomic_accession: {protein_accession: {sequence: str, tool: set, #ofTools: int}}}
    :param connection: open sqlalchemy connection to an SQLite3 db engine
    :param domain_dict: {protein: {fam: set(ranges)}}
    :param protein_db_dict: {protein: db_id}
    :param classifier_db_dict: {classifier: db_id}
    :param family_db_dict: {fam: db_id}
    :param connection: open sqlalchemy connection to an SQLite3 engine
    :param args: cmd-line args parser

    Return nothing.
    """
    logger = logging.getLogger(__name__)

    if args.cazy is not None:
        tools = [
            ['hmmer', 'HMMER'],
            ['hotpep', 'Hotpep'],
            ['diamond', 'DIAMOND'],
            ['dbcan', 'dbCAN'],
            ['cazy', 'CAZy'],
        ]
    else:
        tools = [
            ['hmmer', 'HMMER'],
            ['hotpep', 'Hotpep'],
            ['diamond', 'DIAMOND'],
            ['dbcan', 'dbCAN'],
        ]

    domains_to_insert = set()

    for genomic_accession in tqdm(cazome_dict, desc="Adding domains to db"):
        protein_accessions = list(cazome_dict[genomic_accession].keys())
        for protein_accession in protein_accessions:
            protein_id = protein_db_dict[protein_accession]  # db protein record id

            for tool in tools:
                classifier_id = classifer_db_dict[tool[1]]

                domain_fams = cazome_dict[genomic_accession][protein_accession][tool[0]]

                if tool[0] == 'hmmer':  # includes adding domain ranges
                    for fam in domain_fams:
                        fam_id = family_db_dict[fam]
                        try:
                            domain_ranges = domain_dict[protein_accession][fam]
                        except KeyError:
                            logger.warning(
                                f"HMMER annotation {fam} retrieved for {genomic_accession}:{protein_accession}\n"
                                "but no domain range was retrieved. Adding domain without domain range"
                            )
                            domains_to_insert.add( (protein_id, classifier_id, fam_id, None) )

                        for drange in domain_ranges:
                             domains_to_insert.add( (protein_id, classifier_id, fam_id, drange) )

                else:  # no domain ranges
                    try:
                        for fam in domain_fams:
                            fam_id = family_db_dict[fam]
                            domains_to_insert.add( (protein_id, classifier_id, fam_id, None) )
                    except TypeError:  # raised if domain_fams is None
                        continue
    
    if len(domains_to_insert) != 0:
        insert_data(
            connection,
            'Domains',
            ['protein_id', 'classifier_id', 'family_id', 'domain_range'],
            list(domains_to_insert),
        )

    return
