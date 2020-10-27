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

import os
import subprocess
import sys

from pathlib import Path


def invoke_prediction_tools(query, args, logger):
    """Pass paramaters to CAZyme prediction tools and invoke tool run."""
    # create complete path to fasta file to negate changing cwd
    current_path = os.getcwd()
    fasta_path = query.fasta
    input_path = current_path / fasta_path
    # create complete path to output and directory, to negate changing cwd
    outdir = query.prediction_dir
    output_dir = current_path / outdir

    # variables to store successful/failed process return codes
    dbcan_returncode = 0
    cupp_returncode = 0
    ecami_returncode = 0

    # create list of args to invoke run_dbCAN
    dbcan_args = [
        "run_dbcan.py",
        str(input_path),
        "protein",
        "--out_dir",
        str(output_dir),
    ]

    # change cwd to dbCAN directory to be able to access database files
    os.chdir('tools/dbcan/')

    # create log of dbCAN run
    with open(f"{output_dir}/dbcan.log", "w+") as fh:
        process = subprocess.run(dbcan_args, stdout=fh, text=True)

    # check if successul
    if process.returncode != 0:  # return code is 0 for successful run
        logger.warning(
            (
                f"dbCAN ran into error for {output_dir[-1]}\n."
                "dbCAN error:\n"
                f"{process.stderr}"
            )
        )
        dbcan_returncode = 1

    print("done dbCAN!")
    sys.exit(1)

    # move to cupp directory so can access CUPP
    os.chdir('../')  # moves up to pyrewton/cazymes/prediction/tools
    os.chdir('cupp')

    # create list of args to invoke cupp
    cupp_args = [
        "prediction.py",
        "-query",
        str(fasta_path),
        "-output_path",
        f"{output_dir}/cupp_output.fasta",
    ]

    # create log of CUPP run
    with open(f"{output_dir}/cupp.log", "w+") as fh:
        process = subprocess.run(cupp_args, stdout=fh, text=True)

    # check if successful
    if process.returncode != 0:
        logger.warning(
            (
                f"CUPP ran into error for {output_dir[-1]}\n."
                "CUPP error:\n"
                f"{process.stderr}"
            )
        )
        cupp_returncode = 1

    # move to ecami directory so can access eCAMI
    os.chdir('../')  # moves up to pyrewton/cazymes/prediction/tools
    os.chdir('ecami')

    # create list of args for ecami
    ecami_args = [
        "prediction.py",
        "-input",
        str(fasta_path),
        "-kmer_db",
        "CAZyme",
        "-output",
        f"{output_dir}/ecami_output.txt",
    ]

    # create log of eCAMI run
    with open(f"{output_dir}/ecami.log", "w+") as fh:
        process = subprocess.run(ecami_args, stdout=fh, text=True)

    # check if successful
    if process.returncode != 0:
        logger.warning(
            (
                f"eCAMI ran into error for {output_dir[-1]}\n."
                "eCAMI error:\n"
                f"{process.stderr}"
            )
        )
        ecami_returncode = 1

    # move back to 'predictions/' directory
    os.chdir('../')  # moves to 'tools/'
    os.chdir('../')  # moves to 'predictions/'

    return dbcan_returncode, cupp_returncode, ecami_returncode
