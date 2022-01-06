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
from re import sub
from Bio.SeqUtils.ProtParam import ProteinAnalysis

from tqdm import tqdm

from pyrewton.sql.sql_interface import load_data
from pyrewton.sql.sql_interface.add_data import SqlInterfaceException, insert_data


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
        insert_data(connection, 'Substrates', ['substrate'], substrates_to_insert)
    
    substrates_table_dict = load_data.get_substrates_tables(connection)

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
            binding_sites_to_insert,
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
            active_site_types_inserts
        )

    if len(associated_activities_inserts) != 0:
        insert_data(
            connection,
            'AssociatedActivities',
            ['associated_activity'],
            associated_activities_inserts
        )

    site_type_dict = load_data.get_site_type_dict(connection)
    activity_dict = load_data.get_activity_dict(connection)

    acitve_sites_to_insert = set()

    for site_tuple in active_sites_inserts:
        acitve_sites_to_insert.add(
            (
                site_tuple[0],  # protein_id
                site_tuple[1],  # position
                site_type_dict[site_tuple[2]],  # type_id
                activity_dict[site_tuple[3]],  # activity_id
                site_tuple[4],  # note
                site_tuple[5],  # evidence
            )
        )

    if len(acitve_sites_to_insert) != 0:
        insert_data(
            connection,
            'ActiveSites',
            ['protein_id', 'position', 'type_id', 'activity_id', 'note', 'evidence'],
            acitve_sites_to_insert
        )

    return
    