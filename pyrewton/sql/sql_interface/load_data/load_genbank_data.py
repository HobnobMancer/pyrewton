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


from tqdm import tqdm

from pyrewton.sql.sql_orm import (
    Assembly,
    CazyFamily,
    Classifier,
    Domain,
    Protein,
    Taxonomy,
    Session,
)


def get_complete_domains_table(connection):
    """Parse Domains into dict
    
    :param connection: open sqlalchemy connection to an SQLite3 db engine
    
    Return dict {domain_id : {'protein_id': int, 'family_id': int, 'classifier_id': int, domain_range: str}}
    """
    with Session(bind=connection) as session:
        db_tax_records = session.query(Domain).all()

    domain_data = {}
    for record in tqdm(db_tax_records, desc="Retrieving db Tax records"):
        domain_data[record.domain_id] = {
            'protein_id': record.protein_id,
            'family_id': record.family_id,
            'classifier_id': record.classifier_id,
            'domain_range': record.domain_range,
        }

    return domain_data


def get_protein_annotations_table(connection):
    """Parse Domains into dict
    
    :param connection: open sqlalchemy connection to an SQLite3 db engine
    
    Return dict
    {protein_id : [ (classifier_id, family_id) ]}
    """
    with Session(bind=connection) as session:
        db_tax_records = session.query(Domain).all()

    domain_data = {}
    for record in tqdm(db_tax_records, desc="Retrieving db Tax records"):
        try:

            domain_data[record.protein_id].append( (record.classifier_id, record.family_id) )
        except KeyError:
            domain_data[record.protein_id] = [(record.classifier_id, record.family_id)]

    return domain_data



def get_complete_tax_table(connection):
    """Parse Taxonomies into dict
    
    :param connection: open sqlalchemy connection to an SQLite3 db engine
    
    Return dict {ncbi tax id: {'genus': genus, 'species': species, 'db_id': db_id}
    """
    with Session(bind=connection) as session:
        db_tax_records = session.query(Taxonomy).all()

    tax_data = {}  # {genus species: db_id}
    for record in tqdm(db_tax_records, desc="Retrieving db Tax records"):
        tax_data[record.ncbi_tax_id] = {
            'genus': record.genus,
            'species': record.species,
            'db_id': record.taxonomy_id,
        }

    return tax_data


def get_tax_table(connection):
    """Parse Taxonomies into dict
    
    :param connection: open sqlalchemy connection to an SQLite3 db engine
    
    Return dict {genus species: db_id}
    """
    with Session(bind=connection) as session:
        db_tax_records = session.query(Taxonomy).all()

    tax_data = {}  # {genus species: db_id}
    for record in tqdm(db_tax_records, desc="Retrieving db Tax records"):
        tax_data[f"{record.genus} {record.species}"] = record.taxonomy_id

    return tax_data


def get_complete_assemblies_table(connection):
    """Parse Assemblies into dict
    
    :param connection: open sqlalchemy connection to an SQLite3 db engine
    
    Return dict {genomic_accession: {taxonomy_id: int, db_id: int}}
    """
    with Session(bind=connection) as session:
        db_records = session.query(Assembly).all()

    assembly_dict = {}  # {genomic_accession: db_id}
    for record in tqdm(db_records, desc="Retrieving db Assembly records"):
        assembly_dict[record.assembly_accession] = {
            'taxonomy_id': record.taxonomy_id,
            'db_id': record.assembly_id,
        }

    return assembly_dict


def get_assemblies_table(connection):
    """Parse Assemblies into dict
    
    :param connection: open sqlalchemy connection to an SQLite3 db engine
    
    Return dict {genomic_accession: db_id}
    """
    with Session(bind=connection) as session:
        db_records = session.query(Assembly).all()

    assembly_dict = {}  # {genomic_accession: db_id}
    for record in tqdm(db_records, desc="Retrieving db Assembly records"):
        assembly_dict[record.assembly_accession] = record.assembly_id

    return assembly_dict


def get_protein_db_ids(connection):
    """Get Protein record Ids, load into a dict
    
    :param connection: open sqlalchemy connection to an SQLite3 db engine
    
    Return dict {genbank_accession: db_id}
    """
    with Session(bind=connection) as session:
        db_records = session.query(Protein).all()

    protein_dict = {}  # {genomic_accession: db_id}
    for record in tqdm(db_records, desc="Retrieving db Protein records"):
        protein_dict[record.genbank_accession] = record.protein_id

    return protein_dict


def get_classifier_db_ids(connection):
    """Get Classifier record Ids, load into a dict
    
    :param connection: open sqlalchemy connection to an SQLite3 db engine
    
    Return dict {classifier: db_id}
    """
    with Session(bind=connection) as session:
        db_records = session.query(Classifier).all()

    classifier_dict = {}  # {genomic_accession: db_id}
    for record in tqdm(db_records, desc="Retrieving db classifier records"):
        classifier_dict[record.classifier] = record.classifier_id

    return classifier_dict


def get_family_db_ids(connection):
    """Get CazyFamily record Ids, load into a dict
    
    :param connection: open sqlalchemy connection to an SQLite3 db engine
    
    Return dict {family: db_id}
    """
    with Session(bind=connection) as session:
        db_records = session.query(CazyFamily).all()

    fam_dict = {}  # {genomic_accession: db_id}
    for record in tqdm(db_records, desc="Retrieving db Protein records"):
        fam_dict[record.family] = record.family_id

    return fam_dict


def get_protein_records(connection, protein_accessions):
    """Retrieve records for proteins matching the user's criteria.
    
    :param connection: open sqlalchemy connection to an SQLite3 db engine
    :param protein_accessions: list of protein accessions
    
    Return list of Protein table records.
    """
    with Session(bind=connection) as session:
        db_records = session.query(Protein).all()

    selected_records = [r for r in db_records if r.genbank_accession in protein_accessions]
    
    return selected_records
