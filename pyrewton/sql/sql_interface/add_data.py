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

from pyrewton.sql.sql_interface import SqlInterfaceException


def insert_data(connection, table_name, column_names, insert_values):
    """Insert values into one or multiple rows in the database.
    
    :param connection: open connection to SQLite db engine
    :param table_name: str, name of table to be inserted into
    :param column_names: list of columns (str) to insert data into
    :param insert_values: list of tuples, one tuple per inserted row in the db
    
    Return nothing.
    """
    # set up series of ? to fill in the VALUES statement
    value_stmt = ''
    for name in range((len(column_names)) - 1):
        value_stmt += '?, '
    value_stmt += '?'  # statement should not end with a comma
    
    with connection.begin():
        try:
            connection.exec_driver_sql(
                f"INSERT INTO {table_name} ({', '.join(column_names)}) VALUES ({value_stmt})",
                insert_values,
            )
        except Exception as db_error:
            raise SqlInterfaceException(db_error)

    return


def add_species_data(tax_dict, ncbi_tax_dict, connection):
    species_to_insert = set()
    
    for genomic_accession in tqdm(tax_dict, desc="Adding species to db"):
        tax_id = tax_dict[genomic_accession]['tax_id']
        db_tax_id = ncbi_tax_dict[tax_id]

        genus = tax_dict[genomic_accession]['genus']
        species = tax_dict[genomic_accession]['species']

        species_to_insert.add( (db_tax_id, genus, species) )

    if len(species_to_insert) != 0:
        insert_data(connection, 'Taxonomies', ['ncbi_id', 'genus', 'species'], list(species_to_insert))

    return
