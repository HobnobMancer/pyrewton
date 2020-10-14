#!/usr/bin/env python
# -*- coding: utf-8 -*-
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
"""Module invoking CAZyme prediction tools: dbCAN, CUPP and eCAMI"""

import subprocess


def invoke_prediction_tools(query, args, logger):
    """Pass paramaters to CAZyme prediction tools and invoke tool run."""
    fasta_path = query.fasta
    outdir = query.prediction_dir

    # variables to store successful/failed process return codes
    dbcan_returncode = 0
    cupp_returncode = 0
    ecami_returncode = 0

    # create list of args to invoke run_dbCAN
    dbcan_args = [
        "tools/dbcan/run_dbcan",
        str(fasta_path),
        "protein",
        "--output_dir",
        str(outdir),
    ]

    # create log of dbCAN run
    with open(f"{outdir}/dbcan.log") as fh:
        process = subprocess.run(dbcan_args, stdout=fh, text=True)

    # check if successul
    if process.returncode != 0:  # return code is 0 for successful run
        logger.warning(
            (
                f"dbCAN ran into error for {outdir[-1]}\n."
                "dbCAN error:\n"
                f"{process.stderr}"
            )
        )
        dbcan_returncode = 1

    # create list of args to invoke cupp
    cupp_args = [
        "tools/cupp/prediction.py",
        "-query",
        str(fasta_path),
        "-output_path",
        f"{outdir}/cupp_output.fasta",
    ]

    # create log of CUPP run
    with open(f"{outdir}/cupp.log") as fh:
        process = subprocess.run(cupp_args, stdout=fh, text=True)

    # check if successful
    if process.returncode != 0:
        logger.warning(
            (
                f"CUPP ran into error for {outdir[-1]}\n."
                "CUPP error:\n"
                f"{process.stderr}"
            )
        )
        cupp_returncode = 1

    # create list of args for ecami
    ecami_args = [
        "tools/ecami/prediction.py",
        "-input",
        str(fasta_path),
        "-kmer_db",
        "CAZyme",
        "-output",
        f"{outdir}/ecami_output.txt",
    ]

    # create log of eCAMI run
    with open(f"{outdir}/ecami.log", "w") as fh:
        process = subprocess.run(ecami_args, stdout=fh, text=True)

    # check if successful
    if process.returncode != 0:
        logger.warning(
            (
                f"eCAMI ran into error for {outdir[-1]}\n."
                "eCAMI error:\n"
                f"{process.stderr}"
            )
        )
        ecami_returncode = 1

    return dbcan_returncode, cupp_returncode, ecami_returncode
