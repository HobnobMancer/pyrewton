# /usr/bin/env python
# -*- coding: utf-8 -*-
#!/usr/bin/env bash

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
"""extract_cds.py
Script to match aligned single-copy orthologue protein sequences from
MAFFT alignments with the corresponding CDS sequences, prior to
backtranslation/threading with t-coffee.
"""
import os

from pathlib import Path

from Bio import SeqIO
from tqdm import tqdm

# PROTDIR contains the MAFFT alignments
PROTDIR = Path("oomycetes_orthologs_1/Single_Copy_Orthologue_Sequences")

# CDSDIR contains the CDS/proteome for each genome
CDSDIR = Path("ncbi-genomes-2021-12-14/oomycetes/genomes")

# OUTDIR is the location to which CDS sequences corresponding to the
# protein alignments will be written
OUTDIR = Path("sco_cds")

# add args parser

# The protein sequences have IDs in the form NC_004547.2_1; this is not
# enough to identify the corresponding genome sequence/CDS file. We take
# a brute force approach and load all CDS into memory
print("Loading CDS sequences into memory.")
seqdict = {}
for cdsfpath in tqdm(CDSDIR.iterdir(), total=len(list(CDSDIR.iterdir()))):
    for record in SeqIO.parse(cdsfpath, "genbank"):
       	for (index, feature) in enumerate(record.features):
               if feature.type == "CDS":
                   protein_id = feature.qualifiers["protein_id"][0]
                   cds = feature.location.extract(record).seq
                   seqdict[protein_id] = cds
print(f"Loaded {len(seqdict)} CDS sequences.")

# We iterate over the SCO protein files, and extract the corresponding
# CDS sequences for each, then write these to OUTDIR
print(f"Writing extracted CDS files to {OUTDIR}")
os.makedirs(OUTDIR, exist_ok=True)  # Create output directory, if needed
for protfpath in tqdm(PROTDIR.iterdir()):
    with (OUTDIR / protfpath.stem).with_suffix(".fasta").open("w") as ofh:
        cds = [seqdict[_.id] for _ in SeqIO.parse(protfpath, "fasta")]
        fixed_cds = []
        for _ in cds:
            if str(_.seq[-3:]) not in ("TAA", "TAG", "TGA"):
                fixed_cds.append(_)
            else:
                fixed_cds.append(_[:-3])
        SeqIO.write(
            fixed_cds,
            ofh,
            "fasta",
        )

