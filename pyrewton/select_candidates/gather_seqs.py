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
"""Extract protein sequences from two FASTA files and write to a single FASTA"""


import logging

from typing import List, Optional

from Bio import SeqIO
from saintBioutils.utilities.file_io import make_output_directory
from saintBioutils.utilities.logger import config_logger


def main(args: Optional[List[str]] = None, logger: Optional[logging.Logger] = None):
    """Set up parser and logger, and coordinate prediction of CAZymes."""
    if logger is None:
        config_logger(args)
    logger = logging.getLogger(__package__)

    if args.output is None:
        output_path = "protein_seqs.fasta"
    else:
        if str(args.output.parent) != '.':
            make_output_directory(args.output.parent, args.force, args.nodelete)
        output_path = args.output
    
    seq_record_ids = set()
    seq_records = []
    i = 0
    for record in SeqIO.parse(args.fasta_1, "fasta"):
        record.description=""
        if record.id not in seq_record_ids:
            seq_records.append(record)
            seq_record_ids.add(record.id)
        i += 1
    logger.warning(f"Loaded {i} sequences from: {args.fasta_1}")
    i = 0
    for record in SeqIO.parse(args.fasta_2, "fasta"):
        record.description=""
        if record.id not in seq_record_ids:
            seq_records.append(record)
            seq_record_ids.add(record.id)
        i += 1
    logger.warning(f"Loaded {i} sequences from: {args.fasta_2}")

    logger.warning(f"Writing {len(seq_records)} to {output_path}")
    SeqIO.write(seq_records, output_path, "fasta")


if __name__ == "__main__":
    main()
