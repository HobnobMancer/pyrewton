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
"""Cluster cazyme seqs using MMseqs"""


import logging
import subprocess

from typing import List, Optional

from Bio import SeqIO
from saintBioutils.utilities.file_io import make_output_directory
from saintBioutils.utilities.logger import config_logger


def main(args: Optional[List[str]] = None, logger: Optional[logging.Logger] = None):
    """Set up parser and logger, and coordinate prediction of CAZymes."""
    if logger is None:
        config_logger(args)
    logger = logging.getLogger(__package__)

    if str(args.mmseqs_db.parent) != '.':
        make_output_directory(args.mmseqs_db.parent, args.force, args.nodelete)

    if str(args.mmseqs_out.parent) != '.':
        make_output_directory(args.mmseqs_out.parent, args.force, args.nodelete)    
    
    if str(args.out_tsv.parent) != '.':
        make_output_directory(args.out_tsv.parent, args.force, args.nodelete) 

    create_db_args = [
        "mmseqs",
        "createdb",
        args.fasta,
        args.mmseqs_db,
    ]

    txt = ' '.join(create_db_args)
    print(f"Running command: {txt}")

    theproc = subprocess.Popen([
        "mmseqs",
        "createdb",
        args.fasta,
        args.mmseqs_db,
    ])

    cluster_args = [
        "mmseqs",
        "cluster",
        args.mmseqs_db,
        args.mmseqs_out,
        args.mmseqs_db.parent,
        "--min-seq-id",
        args.pident,
        "-c",
        args.cov,
    ]

    txt = ' '.join(cluster_args)
    print(f"Running command: {txt}")

    theproc = subprocess.Popen([
        "mmseqs",
        "cluster",
        args.mmseqs_db,
        args.mmseqs_out,
        args.mmseqs_db.parent,
        "--min-seq-id",
        args.pident,
        "-c",
        args.cov,
    ])

    createtsv_args = [
        "mmseqs",
        "createtsv",
        args.mmseqs_db,
        args.mmseqs_db,
        args.mmseqs_out,
        args.out_tsv,
    ]

    txt = ' '.join(createtsv_args)
    print(f"Running command: {txt}")

    theproc = subprocess.Popen([
        "mmseqs",
        "createtsv",
        args.mmseqs_db,
        args.mmseqs_db,
        args.mmseqs_out,
        args.out_tsv,
    ])


if __name__ == "__main__":
    main()
