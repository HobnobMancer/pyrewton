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


from tqdm import tqdm

from pyrewton.sql.sql_interface.add_data import insert_data
from pyrewton.sql.sql_interface.load_data import load_genbank_data


def add_families(cazy_families, connection, logger):
    """Add CAZy families to the database
    
    :param cazy_families: set of CAZy families
    :param connection: open connection to sqlite3 db
    :param logger: logger object
    """
    # load existing families = {family: db_id}
    family_table_dict = load_genbank_data.get_family_db_ids(connection)
    fams_in_db = set([fam for fam in family_table_dict])
    # what is in cazy_families but not in fams_in_db
    families_to_add = cazy_families.difference(fams_in_db)  
    families_to_add = [(fam,) for fam in families_to_add]
    if len(families_to_add) != 0:
        logger.warning(f"Adding {len(families_to_add)} CAZy families to the CazyFamilies table")
        insert_data(connection, "CazyFamilies", ['family'], families_to_add)


def add_classifiers(
        classifier,
        connection,
        logger,
        classifier_version=None,
        classifier_training_date=None,
    ):
    """Add classifier to database
    
    :param classifier: str, name of classifier / tool
    :param connection: open connection to SQL db
    :param logger: logger object
    :param classifier version: str, version number
    :param classifier training date: str, release data of CAZy db against which the tool was trained
    """
    db_classifiers = load_genbank_data.get_classifier_db_ids(connection) # {classifier : db id}
    if classifier not in db_classifiers.keys()
        insert_data(
            connection,
            "Classifiers",
            ['classifier', 'version', 'cazy_training_set_date'],
            [(classifier, classifier_version, classifier_training_date)]
        )


def add_classifications(cazy_dict, cazy_classifier_id, connection, logger):
    """Add cazy family classifications - proteins-family relationships
    
    :param cazy_dict: {prot acc : {families}}
    :param cazy_classifier_id: int, local db classifier id for CAZy
    :param connection: open connection to sqlite3 db
    :param logger: logger object
    """
    existing_db_proteins = load_genbank_data.get_protein_db_ids(connection)
    # load existing families = {family: db_id}
    # convert proteins to db protein ids
    family_table_dict = load_genbank_data.get_family_db_ids(connection)

    # existing_db_proteins = {genbank_accession: db_id}
    # {{protein_id : [ (classifier_id, family_id) ]}
    prot_annotation_table = load_genbank_data.get_protein_annotations_table(connection)
    domains_to_add = []  # one domain per cazy family classification for each protein
    for protein_acc in tqdm(cazy_dict, desc="Identifying new CAZy domains to add to db"):  # {protein acc: {families}}
        for family in cazy_dict[protein_acc]:
            family_id = family_table_dict[family]
            if (cazy_classifier_id, family_id) not in prot_annotation_table:
                domains_to_add.append((
                    existing_db_proteins[protein_acc], # protein_id
                    cazy_classifier_id,
                    family_id,
                    None,
                ))

    if len(domains_to_add) != 0:
        logger.warning(f"Adding {len(domains_to_add)} domains from CAZy to the Domains table")
        insert_data(
            connection,
            'Domains',
            ['protein_id','classifier_id','family_id','domain_range'],
            domains_to_add,
        )


