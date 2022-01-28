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


import json
import logging
import urllib.parse
import urllib.request

import numpy as np

from typing import List, Optional
from urllib.error import HTTPError

from bioservices import UniProt
from saintBioutils.utilities.logger import config_logger
from tqdm import tqdm

from pyrewton.cazymes.uniprot import parse_uniprot
from pyrewton.sql.sql_orm import get_cazome_db_connection
from pyrewton.sql.sql_interface.add_data import bulk_insert, add_uniprot_data
from pyrewton.sql.sql_interface.load_data import load_uniprot_data
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
    protein_db_dict = load_uniprot_data.get_protein_db_ids(connection)

    # {uniprot_id: {'gbk_acc': str, 'db_id': int}}
    uniprot_id_dict = get_uniprot_ids(protein_db_dict, args)

    uniprot_cache_path = "uniprot_acc_cache.json"
    with open(uniprot_cache_path, 'w', 'r') as fh:
        json.dump(uniprot_id_dict, fh)

    (
        substrate_binding_inserts,
        glycosylation_inserts,
        temperature_inserts,
        ph_inserts,
        citation_inserts,
        transmembrane_inserts,
        active_sites_inserts,
        active_site_types_inserts,
        associated_activities_inserts,
        metal_binding_inserts,
        metals_inserts,
        cofactors_inserts,
        cofactor_molecules_inserts,
        protein_pdb_inserts,
        pdbs_inserts,
        protein_ec_inserts ,
        ec_inserts,
        protein_table_updates,
    ) = get_uniprot_data(uniprot_id_dict, args)

    add_simple_uniprot_data(
        glycosylation_inserts,
        temperature_inserts,
        ph_inserts,
        citation_inserts,
        transmembrane_inserts,
        connection,
    )

    add_uniprot_data.insert_substrate_data(substrate_binding_inserts, connection)

    add_uniprot_data.insert_active_site_data(
        active_sites_inserts,
        active_site_types_inserts,
        associated_activities_inserts,
        connection,
    )
    
    add_uniprot_data.insert_metal_binding_data(metal_binding_inserts, metals_inserts, connection)

    add_uniprot_data.insert_cofactor_data(
        cofactors_inserts,
        cofactor_molecules_inserts,
        connection,
    )

    add_uniprot_data.insert_pdb_data(
        protein_pdb_inserts,
        pdbs_inserts,
        connection,
    )

    add_uniprot_data.insert_ec_data(
        protein_ec_inserts ,
        ec_inserts,
        connection,
    )

    add_uniprot_data.insert_protein_names(protein_table_updates, connection)

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

    uniprot_rest_queries = get_chunks_list(genbank_accessions, args.batch_size)

    uniprot_gbk_dict = {}  # {uniprot_accession: gbk_accession}
    failed_queries = {}  # {query: tries}

    for query_chunk in tqdm(
        uniprot_rest_queries,
        desc='Batch retrieving UniProt accessions',
    ):
        if type(query_chunk) != str:
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

        for line in tqdm(uniprot_batch_response[1:], "Parsing retrieved UniProt query"):  # the first line includes the titles, last line is an empty str
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


def get_uniprot_data(uniprot_gbk_dict, args):
    """Batch query UniProt to retrieve protein data. Save data to cache directory.
    
    Bioservices requests batch queries no larger than 200.

    :param uniprot_gbk_dict: dict, keyed by GenBank accession and valued by UniProt accession
    :param args: cmd-line args parser
    
    Return
    Dict of data retrieved from UniProt and to be added to the db 
        {uniprot_acccession: {gbk_acccession: str, uniprot_name: str, pdb: set, ec: set}}
    Set of all retrieved EC numbers
    """
    logger = logging.getLogger(__name__)

    uniprot_record_ids = list(uniprot_gbk_dict.keys())
    
    # break up list into nested list of shorter lists for batch querying
    bioservices_queries = get_chunks_list(
        uniprot_record_ids,
        args.batch_size,
    )

    # add PDB column to columns to be retrieved
    UniProt()._valid_columns.append('database(PDB)')

    # data can be inserted straight away for the following tables
    substrate_binding_inserts = set()  # (protein_id, position, note, evidence)
    glycosylation_inserts = set()  # (protein_id, note, evidence)
    temperature_inserts = set()  # (protein_id, lower_opt, upper_opt, lower_therm, upper_therm, lower_lose, upper_lose, note, evidence)
    ph_inserts = set()  # (protein_id, lower_ph, upper_ph, note, evidence)
    citation_inserts = set()  # (protein_id, citation)
    transmembrane_inserts = set()  # (protein_id, transmembrane bool,)

    # data has to be added in stages across the related tables
    active_sites_inserts = set()  # (protein_id, position, type_id, activity_id, note, evidence)
    active_site_types_inserts = set()  # (site_type)
    associated_activities_inserts = set()  # (associated_activity)

    metal_binding_inserts = set()  # (protein_id, ion_id, ion_number, note, evidence)
    metals_inserts = set()  # (ion)

    cofactors_inserts = set()  # (protein_id, molecule_id, note, evidence)
    cofactor_molecules_inserts = set()  # (molecule,)

    protein_pdb_inserts = set()  # (pdb_id, protein_id)
    pdbs_inserts = set()  # (pdb_accession,)

    protein_ec_inserts = set()  # (protein_id,)
    ec_inserts = set()  # (ec_number,)

    protein_table_updates = set()  # (protein_id, uniprot_id, protein_names)

    parsed_uniprot_ids = set()  # IDs of records that have been parsed

    last_renaming_number = len(uniprot_record_ids) - len(parsed_uniprot_ids)

    while len(parsed_uniprot_ids) < len(uniprot_record_ids):

        uniprot_ids_to_parse = [uniprot_id for uniprot_id in uniprot_record_ids if uniprot_id not in parsed_uniprot_ids]

        bioservices_queries = get_chunks_list(
            uniprot_ids_to_parse,
            args.batch_size,
        )

        (
            latest_parsed_uniprot_ids,
            substrate_binding_inserts,
            glycosylation_inserts,
            temperature_inserts,
            ph_inserts,
            citation_inserts,
            transmembrane_inserts,
            active_sites_inserts,
            active_site_types_inserts,
            associated_activities_inserts,
            metal_binding_inserts,
            metals_inserts,
            cofactors_inserts,
            cofactor_molecules_inserts,
            protein_pdb_inserts,
            pdbs_inserts,
            protein_ec_inserts ,
            ec_inserts,
            protein_table_updates,
        ) = query_uniprot(
            bioservices_queries,
            uniprot_record_ids,
            parsed_uniprot_ids,
            uniprot_gbk_dict,
            substrate_binding_inserts,
            glycosylation_inserts,
            temperature_inserts,
            ph_inserts,
            citation_inserts,
            transmembrane_inserts,
            active_sites_inserts,
            active_site_types_inserts,
            associated_activities_inserts,
            metal_binding_inserts,
            metals_inserts,
            cofactors_inserts,
            cofactor_molecules_inserts,
            protein_pdb_inserts,
            pdbs_inserts,
            protein_ec_inserts ,
            ec_inserts,
            protein_table_updates,
        )

        parsed_uniprot_ids = parsed_uniprot_ids.union(latest_parsed_uniprot_ids)
        logger.warning(
            f"UniProt data not retrieved for all UniProt IDs\n"
            f"The {len(uniprot_record_ids)} proteins in the db have a UniProt ID\n"
            f"Data for only {len(parsed_uniprot_ids)} proteins was retrieved"
        )

        new_remaining_number = len(uniprot_record_ids) - len(parsed_uniprot_ids)
        if new_remaining_number == last_renaming_number:
            break
        else:
            last_renaming_number = new_remaining_number
            uniprot_ids_to_parse = [uniprot_id for uniprot_id in uniprot_record_ids if uniprot_id not in parsed_uniprot_ids]
            # # # # # for uniprot_id in uniprot_ids_to_parse:
            # # # # #     logger.warning(f"Data not retrieved for UniProt ID {uniprot_id}")

    # if len(parsed_uniprot_ids) < len(uniprot_record_ids):

    #     # for uniprot_id in uniprot_record_ids:
    #     #     if uniprot_id not in parsed_uniprot_ids:
    #     #         logger.warning(f"Data not retrieved for UniProt ID {uniprot_id}")

    return (
        substrate_binding_inserts,
        glycosylation_inserts,
        temperature_inserts,
        ph_inserts,
        citation_inserts,
        transmembrane_inserts,
        active_sites_inserts,
        active_site_types_inserts,
        associated_activities_inserts,
        metal_binding_inserts,
        metals_inserts,
        cofactors_inserts,
        cofactor_molecules_inserts,
        protein_pdb_inserts,
        pdbs_inserts,
        protein_ec_inserts ,
        ec_inserts,
        protein_table_updates,
    )


def query_uniprot(
    bioservices_queries,
    uniprot_record_ids,
    parsed_uniprot_ids,
    uniprot_gbk_dict,
    substrate_binding_inserts,
    glycosylation_inserts,
    temperature_inserts,
    ph_inserts,
    citation_inserts,
    transmembrane_inserts,
    active_sites_inserts,
    active_site_types_inserts,
    associated_activities_inserts,
    metal_binding_inserts,
    metals_inserts,
    cofactors_inserts,
    cofactor_molecules_inserts,
    protein_pdb_inserts,
    pdbs_inserts,
    protein_ec_inserts ,
    ec_inserts,
    protein_table_updates,
):
    logger = logging.getLogger(__name__)

    for query in tqdm(bioservices_queries, "Batch retrieving protein data from UniProt"):
        uniprot_df = UniProt().get_df(entries=query)

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

            # data for adding to the SubstrateBindingSites table
            substrate_binding_inserts = substrate_binding_inserts.union(
                parse_uniprot.get_sites_data(
                    row,
                    'Binding site',
                    'BINDING',
                    parse_uniprot.get_substrate_binding_site_data,
                    protein_db_id,
                )
            )
            # data for adding to the Glycosylations table
            glycosylation_inserts = glycosylation_inserts.union(
                parse_uniprot.get_sites_data(
                    row,
                    'Glycosylation',
                    'CARBOHYD',
                    parse_uniprot.get_glycosylation_data,
                    protein_db_id,
                )
            )
            # data for adding to the Temperatures table
            temperature_inserts = temperature_inserts.union(
                parse_uniprot.get_temp_data(row, protein_db_id)
            )
            # data for adding to the OptimalPHs table
            ph_inserts = ph_inserts.union(
                parse_uniprot.get_optimum_ph(row, protein_db_id)
            )
            # data for adding to the Citations table
            citation_inserts = citation_inserts.union(
                parse_uniprot.get_citations(row, protein_db_id)
            )
            
            # data for adding to the Transmembranes table
            try:
                np.isnan(row['Transmembrane'])
                transmembrane_inserts.add( (protein_db_id, False,) )
            except TypeError:
                # raised if value is returened for 'Transmembrane'
                # becuase transmembrane region is annotated
                transmembrane_inserts.add( (protein_db_id, True,) )
            
            # data to be added to the ActiveSites, SiteTypes and AssociatedActivities tables
            new_active_sites, new_site_types, new_activities = parse_uniprot.get_active_site_data(
                row,
                protein_db_id,
            )

            active_sites_inserts = active_sites_inserts.union(new_active_sites)
            active_site_types_inserts = active_site_types_inserts.union(new_site_types)
            associated_activities_inserts = associated_activities_inserts.union(new_activities)

            # data to be added to Metals and MetalBindingSites tables
            new_metal_sites, new_metals = parse_uniprot.get_metal_binding_sites(
                row,
                protein_db_id,
            )
            metal_binding_inserts = metal_binding_inserts.union(new_metal_sites)
            metals_inserts = metals_inserts.union(new_metals)

            # data to be add to Cofactors
            new_cofactors, new_molecules = parse_uniprot.get_cofactor_data(
                row,
                protein_db_id,
            )
            cofactors_inserts = cofactors_inserts.union(new_cofactors)
            cofactor_molecules_inserts = cofactor_molecules_inserts.union(new_molecules)

            new_protein_pdb_inserts, new_pdb_insert = parse_uniprot.get_pdb_ecs(
                row,
                protein_db_id,
                'Cross-reference (PDB)',
            )
            protein_pdb_inserts = protein_pdb_inserts.union(new_protein_pdb_inserts)
            pdbs_inserts = pdbs_inserts.union(new_pdb_insert)

            new_protein_ec_inserts, new_ec_inserts = parse_uniprot.get_pdb_ecs(
                row,
                protein_db_id,
                'EC number',
            )
            protein_ec_inserts = protein_ec_inserts.union(new_protein_ec_inserts)
            ec_inserts = ec_inserts.union(new_ec_inserts)

            parsed_uniprot_ids.add(uniprot_acc)

    return (
        parsed_uniprot_ids,
        substrate_binding_inserts,
        glycosylation_inserts,
        temperature_inserts,
        ph_inserts,
        citation_inserts,
        transmembrane_inserts,
        active_sites_inserts,
        active_site_types_inserts,
        associated_activities_inserts,
        metal_binding_inserts,
        metals_inserts,
        cofactors_inserts,
        cofactor_molecules_inserts,
        protein_pdb_inserts,
        pdbs_inserts,
        protein_ec_inserts ,
        ec_inserts,
        protein_table_updates,
    )


def add_simple_uniprot_data(
        glycosylation_inserts,
        temperature_inserts,
        ph_inserts,
        citation_inserts,
        transmembrane_inserts,
        connection,
):
    """Bulk insert data into tables, which require no additional parsing"""
    if len(glycosylation_inserts) != 0:
        try:
            bulk_insert.insert_data(
                connection,
                'Glycosylations',
                ['protein_id', 'position', 'note', 'evidence'],
                list(glycosylation_inserts),
            )
        except Exception:
            print("GLY ERROR:", )

    if len(temperature_inserts) != 0:
        bulk_insert.insert_data(
            connection,
            'Temperatures',
            [
                'protein_id',
                'lower_optimal',
                'upper_optimal',
                'lower_thermostable',
                'upper_thermostable',
                'lower_lose_activity',
                'upper_lose_activity',
                'note',
                'evidence',
            ],
            list(temperature_inserts),
        )

    if len(ph_inserts) != 0:
        bulk_insert.insert_data(
            connection,
            'OptimalPHs',
            ['protein_id', 'lower_pH', 'upper_pH', 'note', 'evidence'],
            list(ph_inserts),
        )

    if len(citation_inserts) != 0:
        bulk_insert.insert_data(
            connection,
            'Citations',
            ['protein_id', 'citation'],
            list(citation_inserts),
        )

    if len(transmembrane_inserts) != 0:
        bulk_insert.insert_data(
            connection,
            'Transmembranes',
            ['protein_id', 'uniprot_transmembrane'],
            list(transmembrane_inserts),
        )

    return


if __name__ == "__main__":
    main()