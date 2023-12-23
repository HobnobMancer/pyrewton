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

from tqdm import tqdm

from pyrewton.sql.sql_orm import (
    ActiveSiteType,
    AssociatedActivity,
    CofactorMolecule,
    Ec_number,
    Metal,
    Pdb,
    Protein,
    Session,
    Substrate,
)


def get_protein_db_ids(connection):
    """Parse Proteins table into dict:
    
    :param:
    
    Return {genbank_accession: db_id}
    """
    with Session(bind=connection) as session:
        db_records = session.query(Protein).all()

    data = {}  # {genbank_accession: db_id}
    for record in tqdm(db_records, desc="Retrieving db Protein records"):
        data[record.genbank_accession] = record.protein_id

    return data


def get_substrates_tables(connection):
    """Parse Substrates table into dict
    
    :param:
    
    Return dict {substrate: substrate_id}
    """
    with Session(bind=connection) as session:
        db_records = session.query(Substrate).all()

    data = {}  # {substrate: substrate_id}
    for record in tqdm(db_records, desc="Retrieving db Substrate records"):
        data[record.substrate] = record.substrate_id

    return data


def get_site_type_dict(connection):
    with Session(bind=connection) as session:
        db_records = session.query(ActiveSiteType).all()

    data = {}  # {site_type: db_id}
    for record in tqdm(db_records, desc="Retrieving db Site Types records"):
        data[record.site_type] = record.type_id

    return data


def get_activity_dict(connection):
    with Session(bind=connection) as session:
        db_records = session.query(AssociatedActivity).all()

    data = {}  # {site_type: db_id}
    for record in tqdm(db_records, desc="Retrieving db Acitivity records"):
        data[record.associated_activity] = record.activity_id

    return data


def load_metals_table(connection):
    """Parse Metals table into dict
    
    :param connection: open sqlalchemy connection to an SQLite3 db engine
    
    Return dict {metal: db_id}
    """
    with Session(bind=connection) as session:
        db_records = session.query(Metal).all()

    metal_data = {}  # {metal: db_id}
    for record in tqdm(db_records, desc="Retrieving db Metal Ion records"):
        metal_data[record.ion] = record.ion_id

    return metal_data


def get_molecules_dict(connection):
    """Parse CofactorMolecules table into dict
    
    :param:
    
    Return dict {molecule: molecule_id}
    """
    with Session(bind=connection) as session:
        db_records = session.query(CofactorMolecule).all()

    data = {}  # {molecule: molecule_id}
    for record in tqdm(db_records, desc="Retrieving db Cofactor Molecule records"):
        data[record.molecule] = record.molecule_id

    return data


def get_pdb_dict(connection):
    """Parse Pdbs table into dict
    
    :param:
    
    Return dict {pdb: pdb_id}
    """
    with Session(bind=connection) as session:
        db_records = session.query(Pdb).all()

    data = {}  # {pdb: molecule_id}
    for record in tqdm(db_records, desc="Retrieving db PDB records"):
        data[record.pdb_accession] = record.pdb_id

    return data


def get_ec_dict(connection):
    """Parse Ec_numbers table into dict
    
    :param:
    
    Return dict {ec: ec_id}
    """
    with Session(bind=connection) as session:
        db_records = session.query(Ec_number).all()

    data = {}  # {pdb: molecule_id}
    for record in tqdm(db_records, desc="Retrieving db EC number records"):
        data[record.ec_number] = record.ec_id

    return data
