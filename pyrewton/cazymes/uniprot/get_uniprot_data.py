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


if __name__ == "__main__":
    main()
