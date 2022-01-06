#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) University of St Andrews 2020-2021
# (c) University of Strathclyde 2020-2021
# (c) James Hutton Institute 2020-2021
#
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
"""Retrieve data from UniProt for CAZymes in a CAZome database"""


import logging
import urllib.parse
import urllib.request

from typing import List, Optional
from urllib.error import HTTPError

from bioservices import UniProt
from saintBioutils.utilities.logger import config_logger
from tqdm import tqdm

from pyrewton.sql.sql_orm import get_cazome_db_connection
from pyrewton.sql.sql_interface import load_data
from pyrewton.utilities.parsers.cmd_parser_add_uniprot import build_parser


def main(argv: Optional[List[str]] = None, logger: Optional[logging.Logger] = None):
    """Set up parser and logger, and coordinate prediction of CAZymes."""
    if argv is None:
        parser = build_parser()
        args = parser.parse_args()
    else:
        parser = build_parser(argv)
        args = parser.parse_args()

    if logger is None:
        config_logger(args)
    logger = logging.getLogger(__package__)

    connection = get_cazome_db_connection(args.database, args)

    # {genbank_accession: db_id}
    protein_db_dict = load_data.get_protein_db_ids(connection)

    # {uniprot_id: {'gbk_acc': str, 'db_id': int}}
    uniprot_id_dict = get_uniprot_ids(protein_db_dict)

    return


def get_uniprot_ids(protein_db_dict, args):
    """Batch query UniProt to get the UniProt record IDs

    Retrieve UniProt accessions for the GenBank accessions from UniProt REST API.
    
    UniProt requests batch queries of no larger than 20,000, athough queries longer than 500
    often raise HTTP 400 Error codes, especially in busy server times.
    
    :param protein_db_dict: {genbank_accession: db_id}
    :param args: cmd-line args parser

    Return dict {uniprot_id: {'gbk_acc': str, 'db_id': int}}
    Contains only proteins in UniProt and the local CAZome db
    """
    logger = logging.getLogger(__name__)
    uniprot_url = 'https://www.uniprot.org/uploadlists/'

    genbank_accessions = list(protein_db_dict.keys())

    uniprot_rest_queries = get_chunks_list(genbank_accessions, args.uniprot_batch_size)

    uniprot_gbk_dict = {}  # {uniprot_accession: gbk_accession}
    failed_queries = {}  # {query: tries}

    for query_chunk in tqdm(
        uniprot_rest_queries,
        desc='Batch retrieving UniProt accessions',
    ):
        if type(query) != str:
            # convert the set of gbk accessions into str format
            query = ' '.join(query_chunk)

        params = {
            'from': 'EMBL',
            'to': 'ACC',
            'format': 'tab',
            'query': query
        }

        # submit query data
        data = urllib.parse.urlencode(params)
        data = data.encode('utf-8')
        req = urllib.request.Request(uniprot_url, data)

        # retrieve UniProt response
        try:
            with urllib.request.urlopen(req) as f:
                response = f.read()
        except HTTPError:
            try:
                failed_queries[query] += 1
            except KeyError:
                failed_queries[query] = 1
            if failed_queries[query] > args.retries:
                del failed_queries[query]
            else:
                uniprot_rest_queries.append(query)

        uniprot_batch_response = response.decode('utf-8')

        uniprot_batch_response = uniprot_batch_response.split('\n')

        for line in uniprot_batch_response[1:]:  # the first line includes the titles, last line is an empty str
            if line == '':  # add check incase last line is not an empty str 
                continue
            uniprot_accession = line.split('\t')[1]
            genbank_accession = line.split('\t')[0]
            db_id = protein_db_dict[genbank_accession]
            uniprot_gbk_dict[uniprot_accession] = {'gbk_acc': genbank_accession, 'db_id': db_id}

    logger.info(
        f"Retrieved {len(genbank_accessions)} gbk accessions from the local db\n"
        f"{len(list(uniprot_gbk_dict.keys()))} were assoicated with records in UniProt"
    )

    return uniprot_gbk_dict


def get_chunks_list(lst, chunk_length):
    """Separate the long list into separate chunks.
    :param lst: list to be separated into smaller lists (or chunks)
    :param chunk_length: int, the length of the lists the longer list is to be split up into
    Return a list of nested lists.
    """
    chunks = []
    for i in range(0, len(lst), chunk_length):
        chunks.append(lst[i:i + chunk_length])
    return chunks


def get_uniprot_data(uniprot_gbk_dict, cache_dir, args):
    """Batch query UniProt to retrieve protein data. Save data to cache directory.
    
    Bioservices requests batch queries no larger than 200.
    :param uniprot_gbk_dict: dict, keyed by GenBank accession and valued by UniProt accession
    :param cache_dir: path to directory to write out cache
    :param args: cmd-line args parser
    
    Return
    Dict of data retrieved from UniProt and to be added to the db 
        {uniprot_acccession: {gbk_acccession: str, uniprot_name: str, pdb: set, ec: set}}
    Set of all retrieved EC numbers
    """
    logger = logging.getLogger(__name__)
    
    # break up list into nested list of shorter lists for batch querying
    bioservices_queries = get_chunks_list(
        list(uniprot_gbk_dict.keys()),
        args.bioservices_batch_size,
    )

    # data can be inserted straight away for the following tables
    substrate_binding_inserts = set()  # (protein_id, position, note, evidence)
    glycosylation_inserts = set()  # (protein_id, note, evidence)
    temperature_inserts = set()  # (protein_id, lower_opt, upper_opt, lower_therm, upper_therm, lower_lose, upper_lose, note, evidence)
    ph_inserts = set()  # (protein_id, lower_ph, upper_ph, note, evidence)
    citation_inserts = set()  # (protein_id, citation)
    transmembrane_inserts = set()  # (protein_id, transmembrane bool,)
    uniprot_protein_data = set()  # {protein_id: (uniprot_acc, uniprot_name)}

    # data has to be added in stages across the related tables
    active_sites_inserts = set()  # (protein_id, position, type_id, activity_id, note, evidence)
    active_site_types_inserts = set()  # (site_type)
    associated_activities_inserts = set()  # (associated_activity)

    metal_binding_inserts = set()  # (protein_id, ion_id, ion_number, note, evidence)
    metals_inserts = set()  # (ion)

    cofactors_inserts = set()  # (protein_id, molecule_id, note, evidence)
    cofactor_molecules_inserts = []   # (molecule,)

    protein_pdb_inserts = set()  # (pdb_id, protein_id)
    pdbs_inserts = set()  # (pdb_accession,)

    protein_ec_inserts = set()  # (protein_id,)
    ec_inserts = set()  # (ec_number,)

    protein_table_updates = set()  # (protein_id, uniprot_id, protein_names)

    uniprot_record_ids = list(uniprot_gbk_dict.keys())
    parsed_uniprot_ids = set()  # IDs of records that have been parsed

    for query in tqdm(bioservices_queries, "Batch retrieving protein data from UniProt"):
        uniprot_df = UniProt().get_df(entries=query)
       
        # filter for the columns of interest
        uniprot_df = uniprot_df[[
            'Entry',  ## UniProt record ID
            'Protein names',
            'Active site',
            'Binding site',
            'Metal binidng',
            'Cofactor',
            'Temperature dependence',
            'pH dependence',
            'PubMed ID',
            'EC number',
            'Cross-reference (PDB)',
        ]]

        index = 0
        for index in tqdm(range(len(uniprot_df['Entry'])), desc="Parsing UniProt response"):
            row = uniprot_df.iloc[index]

            uniprot_acc = row['Entry']

            if uniprot_acc not in uniprot_record_ids:
                continue
            
            uniprot_name = row['Protein names']

            protein_db_id = uniprot_gbk_dict[uniprot_acc]['db_id']

            # checked if parsed before incase bioservices returned duplicate proteins
            if uniprot_acc in parsed_uniprot_ids:
                logger.warning(
                    f'Multiple entries for UniProt:{uniprot_acc}, '
                    f'GenBank:{uniprot_gbk_dict[uniprot_acc]} retrieved from UniProt,\n'
                    'Add ALL data to the database'
                )

            # data for adding to the Proteins table
            protein_table_updates.add( (protein_db_id, uniprot_acc, uniprot_name) )



    return


if __name__ == "__main__":
    main()
