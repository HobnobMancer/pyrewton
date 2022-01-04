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
    Ncbi_Tax_Id,
    Taxonomy,
    Session,
)


def get_ncbi_tax_ids(connection):
    """Parse Ncbi_Tax_Ids into a dict
    
    :param connection: open sqlalchemy connection to an SQLite3 db engine
    
    Return dict {ncbi_tax_id: db_id}
    """
    with Session(bind=connection) as session:
        db_ncbi_records = session.query(Ncbi_Tax_Id).all()

    ncbi_tax_ids = {}  # {ncbi_tax_id: db_id}
    for record in tqdm(db_ncbi_records, desc="Retrieving NCBI db records"):
        ncbi_tax_ids[record.ncbi_tax_id] = record.ncbi_id

    return ncbi_tax_ids


def get_tax_table(connection):
    """Parse Taxonomies into dict
    
    :param connection: open sqlalchemy connection to an SQLite3 db engine
    
    Return dict {genus species: db_id}
    """
    with Session(bind=connection) as session:
        db_tax_records = session.query(Taxonomy).all()

    tax_data = {}  # {genus species: db_id}