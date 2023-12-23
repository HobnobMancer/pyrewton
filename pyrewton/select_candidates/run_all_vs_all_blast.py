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
"""Explore sequence diversity by running all versus all BLAST pairwise alignments"""


import logging
import subprocess


from typing import List, Optional

from Bio.Blast.Applications import NcbiblastpCommandline
from saintBioutils.utilities.file_io import make_output_directory
from saintBioutils.utilities.logger import config_logger


def main(args: Optional[List[str]] = None, logger: Optional[logging.Logger] = None):
    """Set up parser and logger, and coordinate prediction of CAZymes."""
    if logger is None:
        config_logger(args)
    logger = logging.getLogger(__package__)

    if str(args.output.parent) != '.':
        make_output_directory(args.output.parent, args.force, args.nodelete)

    if str(args.db_path.parent) != '.':
        make_output_directory(args.db_path.parent, args.force, args.nodelete)

    if args.method == 'BLAST':
        logger.warning("Running BLAST")

        all_v_all_blastp = NcbiblastpCommandline(
            query=args.input,
            subject=args.input,
            out=args.output,
            evalue=args.evalue,
            outfmt="6 qseqid sseqid pident qcovs qlen slen length bitscore evalue",
        )
        stdout, stderr = all_v_all_blastp()

        # check if alignment was successful
        if len(stderr) != 0:
            print(stderr)

        logger.warning(f"Written alignemnt output to:\n{args.output}")

    elif args.method == 'DIAMOND':
        logger.warning("Running DIAMOND")

        create_db = [
            "diamond",
            "makedb",
            "--in",
            args.input,
            "--db",
            args.db_path,
        ]

        txt = ' '.join(create_db)
        print(f"Running command: {txt}")

        theproc = subprocess.Popen([
            "diamond",
            "makedb",
            "--in",
            args.input,
            "--db",
            args.db_path,
        ])

        run_diamond = [
            "diamond",
            "blastp",
            "--db",
            args.db_path,
            "--query",
            args.input,
            "--out",
            args.output,
            "--outfmt",
            "6 qseqid sseqid qlen slen length pident evalue bitscore",
            "--evalue",
            args.evalue,
            "--max-target-seqs",
            0,
        ]

        txt = ' '.join(run_diamond)
        print(f"Running command: {txt}")

        theproc = subprocess.Popen([
            "diamond",
            "blastp",
            "--db",
            args.db_path,
            "--query",
            args.input,
            "--out",
            args.output,
            "--outfmt",
            "6 qseqid sseqid qlen slen length pident evalue bitscore",
            "--evalue",
            args.evalue,
            "--max-target-seqs",
            0,
        ])


if __name__ == "__main__":
    main()
