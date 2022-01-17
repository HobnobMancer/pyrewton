#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) University of St Andrews 2022
# (c) University of Strathclyde 2022
# (c) James Hutton Institute 2022
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
"""Extract the sequences of proteins matching predefined criteria from a db, and write to FASTA"""


import logging

from typing import List, Optional

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from saintBioutils.utilities.file_io import make_output_directory
from saintBioutils.utilities.logger import config_logger
from tqdm import tqdm

from pyrewton.sql.sql_orm import get_cazome_db_connection
from pyrewton.sql.sql_interface.load_data.load_genbank_data import get_protein_records
from pyrewton.utilities.parsers.cmd_parser_extract_db_seq import build_parser


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

    if args.output is None:
        output_path = "cazome_protein_seqs.fasta"
    else:
        if str(args.output.parent) != '.':
            make_output_directory(args.output.parent, args.force, args.nodelete)
        output_path = args.output

    connection = get_cazome_db_connection(args.database, args)

    protein_accessions = get_user_protein_accessions(args)

    protein_records = get_protein_records(connection, protein_accessions)
    
    seq_records = []
    
    for protein_record in tqdm(protein_records, "Compiling SeqRecords"):
        if protein_record.uniprot_id is None:
            uniprot_id = "Unknown UniProt ID"
        else:
            uniprot_id = protein_record.uniprot_id
        
        new_record = SeqRecord(
            Seq((protein_record.sequence).strip()),
            id=protein_record.genbank_accession,
            description=uniprot_id,
        )
        seq_records.append(new_record)

    SeqIO.write(seq_records, output_path, "fasta")


def get_user_protein_accessions(args):
    """Retrieve protein accessions provided by the user in a txt file.
    
    :param args: cmd-line args parser
    
    Return list of unique protein accessions.
    """
    logger = logging.getLogger(__name__)

    try:
        with open(args.proteins, "r") as fh:
            lines = fh.read().splitlines()
    except FileNotFoundError:
        logger.error(
            f"Could not find file of protein accessions at {args.proteins}\nTerminating program"
        )

    protein_accessions = [line.strip() for line in lines]

    return list(set(protein_accessions))


if __name__ == "__main__":
    main()
