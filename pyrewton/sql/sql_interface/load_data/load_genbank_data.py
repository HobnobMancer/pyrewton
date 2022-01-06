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
    Protein,
    Taxonomy,
    Session,
)



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
    
    Return dict {genbank_accession: db_id}
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
    
    Return dict {genbank_accession: db_id}
    """
    with Session(bind=connection) as session:
        db_records = session.query(CazyFamily).all()

    fam_dict = {}  # {genomic_accession: db_id}
    for record in tqdm(db_records, desc="Retrieving db Protein records"):
        fam_dict[record.family] = record.family_id

    return fam_dict
