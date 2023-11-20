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


from sqlalchemy import text
from tqdm import tqdm

from pyrewton.sql.sql_interface.load_data import load_uniprot_data
from pyrewton.sql.sql_interface.add_data.bulk_insert import insert_data


def insert_substrate_data(substrate_binding_inserts, connection):
    """Insert substrate data into the SubstrateBindingSites and Substrates tables
    
    :param substrate_binding_inserts: set of tuples (protein_id, position, substrate, evidence)
    :param connection:
    
    Return nothing.
    """
    substrates_to_insert = set()
    
    for substrate_tuple in substrate_binding_inserts:
        substrates_to_insert.add( (substrate_tuple[2],) )

    if len(substrates_to_insert) != 0:
        insert_data(connection, 'Substrates', ['substrate'], list(substrates_to_insert))
    
    substrates_table_dict = load_uniprot_data.get_substrates_tables(connection)

    binding_sites_to_insert = set()

    for substrate_tuple in substrate_binding_inserts:
        binding_sites_to_insert.add(
            (
                substrate_tuple[0],  # protein_id
                substrate_tuple[1],  # position
                substrates_table_dict[substrate_tuple[2]],  # substrate_id
                substrate_tuple[3],  # evidence
            )
        )

    if len(binding_sites_to_insert) != 0:
        insert_data(
            connection,
            'SubstrateBindingSites',
            ['protein_id', 'position', 'substrate_id', 'evidence'],
            list(binding_sites_to_insert),
        )

    return


def insert_active_site_data(
    active_sites_inserts,
    active_site_types_inserts,
    associated_activities_inserts,
    connection,
):
    """Insert data pertaining to active sites into the db
    
    :param:
    
    Return nothing
    """
    if len(active_site_types_inserts) != 0:
        insert_data(
            connection,
            'SiteTypes',
            ['site_type'],
            list(active_site_types_inserts),
        )

    if len(associated_activities_inserts) != 0:
        insert_data(
            connection,
            'AssociatedActivities',
            ['associated_activity'],
            list(associated_activities_inserts),
        )

    site_type_dict = load_uniprot_data.get_site_type_dict(connection)
    activity_dict = load_uniprot_data.get_activity_dict(connection)

    acitve_sites_to_insert = set()

    for site_tuple in active_sites_inserts:
        if site_tuple[2] is None:
            site_type = site_tuple[2]
        else:
            site_type = site_type_dict[site_tuple[2]]
        
        if site_tuple[3] is None:
            activity = site_tuple[3]
        else:
            activity = activity_dict[site_tuple[3]]

        acitve_sites_to_insert.add(
            (
                site_tuple[0],  # protein_id
                site_tuple[1],  # position
                site_type,  # type_id
                activity,  # activity_id
                site_tuple[4],  # note
                site_tuple[5],  # evidence
            )
        )

    if len(acitve_sites_to_insert) != 0:
        insert_data(
            connection,
            'ActiveSites',
            ['protein_id', 'position', 'type_id', 'activity_id', 'note', 'evidence'],
            list(acitve_sites_to_insert),
        )

    return


def insert_metal_binding_data(metal_binding_inserts, metals_inserts, connection):
    """Insert data pertaining to metal binding into the db
    
    :param:
    
    Return nothing
    """
    metals = [_ for _ in metals_inserts]
    if len(metals) != 0:
        insert_data(
            connection,
            'Metals',
            ['ion'],
            metals,
        )

    metals_table_dict = load_uniprot_data.load_metals_table(connection)

    metal_sites_to_insert = set()

    for site_tuple in metal_binding_inserts:
        if site_tuple[2] is None:
            ion = None
        else:
            ion = metals_table_dict[site_tuple[2]]

        metal_sites_to_insert.add(
            (
                site_tuple[0],  # protein_id
                site_tuple[1],  # position
                ion,  # ion_id
                site_tuple[3],  # ion_number,
                site_tuple[4],  # note
                site_tuple[5],  # evidence
            )
        )

    if len(metal_sites_to_insert) != 0:
        insert_data(
            connection,
            'MetalBindingSites',
            ['protein_id', 'position', 'ion_id', 'ion_number', 'note', 'evidence'],
            list(metal_sites_to_insert),
        )
    
    return


def insert_cofactor_data(
    cofactors_inserts,
    cofactor_molecules_inserts,
    connection,
):
    """Add data about cofactors to the db
    
    :param
    
    Return nothing
    """
    if len(cofactor_molecules_inserts) != 0:
        insert_data(
            connection,
            'CofactorMolecules',
            ['molecule'],
            list(cofactor_molecules_inserts),
        )

    molecules_dict = load_uniprot_data.get_molecules_dict(connection)

    cofactor_data_inserts = set()

    for data_tuple in cofactors_inserts:
        if data_tuple[1] is None:
            molecule = None
        else:
            molecule = molecules_dict[data_tuple[1]]

        cofactor_data_inserts.add(
            (
                data_tuple[0],  # protein_id
                molecule,  # molecule_id
                data_tuple[2],  # note
                data_tuple[3],  # evidence
            )
        )
    
    if len(cofactor_data_inserts) != 0:
        insert_data(
            connection,
            'Cofactors',
            ['protein_id', 'molecule_id', 'note', 'evidence'],
            list(cofactor_data_inserts),
        )

    return


def insert_pdb_data(
    protein_pdb_inserts,
    pdbs_inserts,
    connection,
):
    """Insert PDB accessions and relation to proteins to the db
    
    :param
    
    Return nothing
    """
    if len(pdbs_inserts) != 0:
        insert_data(
            connection,
            'Pdbs',
            ['pdb_accession'],
            list(pdbs_inserts),
        )

    pdb_table_dict = load_uniprot_data.get_pdb_dict(connection)

    pdb_relationships = set()

    for data_tuple in protein_pdb_inserts:
        if data_tuple[1] is None:
            pdb = None
        else:
            pdb = pdb_table_dict[data_tuple[1]]

        pdb_relationships.add(
            (
                data_tuple[0],  # protein_id
                pdb,  # pdb_id
            )
        )

    if len(pdb_relationships) != 0:
        insert_data(
            connection,
            'Proteins_Pdbs',
            ['protein_id', 'pdb_id'],
            list(pdb_relationships),
        )

    return


def insert_ec_data(
    protein_ec_inserts ,
    ec_inserts,
    connection,
):
    """Add EC numbers and relation to proteins to the db
    
    :param:
    
    Return nothing
    """
    if len(ec_inserts) != 0:
        insert_data(
            connection,
            'Ec_numbers',
            ['ec_number'],
            list(ec_inserts),
        )

    ec_table_dict = load_uniprot_data.get_ec_dict(connection)

    ec_relationships = set()

    for data_tuple in protein_ec_inserts:
        if data_tuple[1] is None:
            ec = None
        else:
            ec = ec_table_dict[data_tuple[1]]

        ec_relationships.add(
            (
                data_tuple[0],  # protein_id
                ec,  # ec_id
            )
        )

    if len(ec_relationships) != 0:
        insert_data(
            connection,
            'Proteins_Ecs',
            ['protein_id', 'ec_id'],
            list(ec_relationships),
        )
    
    return


def insert_protein_names(protein_table_updates, connection):
    """Add UniProt protein names and IDs to the Proteins table
    
    :param protein_table_updates:
    :param connection:
    
    Return nothing.
    """
    for record in tqdm(protein_table_updates, desc="Updating UniProt data in the Proteins table"):
        connection.execute(
            text(
                f"UPDATE Proteins SET uniprot_id = '{record[1]}' WHERE protein_id = '{record[0]}'"
            )
        )
        protein_name = record[2].replace('"', "")
        protein_name = protein_name.replace("'", "")
        protein_name = protein_name.replace("`", "")
        connection.execute(
            text(
                f"UPDATE Proteins SET protein_name = '{protein_name}' WHERE protein_id = '{record[0]}'"
            )
        )

    return