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


def invoke_prediction_tools(query, logger):
    """Pass paramaters to CAZyme prediction tools and invoke tool run.

    :param query: Query class instance, input fasta file for tools
    :param logger: logger object

    Return full/absolute path to directory containing output files for all predictions for input FASTA file.
    """
    # create complete path to fasta file to negate changing cwd
    current_path = os.getcwd()
    fasta_path = query.fasta
    input_path = current_path / fasta_path

    # create complete path to output and directory, to negate changing cwd
    out_dir = query.prediction_dir
    out_dir = current_path / out_dir

    # change cwd to dbCAN directory to be able to access database files
    os.chdir('tools/dbcan/')
    print("dbCAN is predicting CAZymes")
    invoke_dbcan(input_path, out_dir, logger)

    # move to cupp directory so can access CUPP
    os.chdir('..')  # moves up to pyrewton/cazymes/prediction/tools
    os.chdir('cupp')
    print("CUPP is predicting CAZymes")
    invoke_cupp(input_path, out_dir, logger)

    # move to ecami directory so can access eCAMI
    os.chdir('..')  # moves up to pyrewton/cazymes/prediction/tools
    os.chdir('ecami')
    print("eCAMI is predicting CAZymes")
    invoke_ecami(input_path, out_dir, logger)

    # move back to 'predictions/' directory
    os.chdir('..')  # moves to 'tools/'
    os.chdir('..')  # moves to 'predictions/'

    return out_dir


def invoke_dbcan(input_path, out_dir, logger):
    """Invoke the prediction tool (run-)dbCAN.

    :param input_path: path to input FASTA file
    :param out_dir: path to output directory for input FASTA file query
    :param logger: logger object

    Return nothing
    """
    # create list of args to invoke run_dbCAN
    dbcan_args = [
        "run_dbcan.py",
        str(input_path),
        "protein",
        "--out_dir",
        str(out_dir),
    ]
    print("dbCAN here!!!")
    # # create log of dbCAN run
    # with open(f"{out_dir}/dbcan.log", "w+") as fh:
    #     process = subprocess.run(dbcan_args, stdout=fh, text=True)

    # # check if successul
    # if process.returncode != 0:  # return code is 0 for successful run
    #     logger.warning(
    #         (
    #             f"dbCAN ran into error for {out_dir}\n."
    #             "dbCAN error:\n"
    #             f"{process.stderr}"
    #         )
    #     )

    return


def invoke_cupp(input_path, out_dir, logger):
    """Invoke the prediction tool CUPP.

    :param input_path: path to input FASTA file
    :param out_dir: path to output directory for input FASTA file query
    :param logger: logger object

    Return nothing
    """
    # create list of args to invoke cupp
    cupp_args = [
        "python3",
        "CUPPprediction.py",
        "-query",
        str(input_path),
        "-output_path",
        f"{out_dir}/cupp_output.fasta",
    ]

    # create log of CUPP run
    with open(f"{out_dir}/cupp.log", "w+") as fh:
        process = subprocess.run(cupp_args, stdout=fh, text=True)

    # check if successful
    if process.returncode != 0:
        logger.warning(
            (
                f"CUPP ran into error for {out_dir}\n."
                "CUPP error:\n"
                f"{process.stderr}"
            )
        )

    return


def invoke_ecami(input_path, out_dir, logger):
    """Invoke the prediction tool eCAMI.

    :param input_path: path to input FASTA file
    :param out_dir: path to output directory for input FASTA file query
    :param logger: logger object

    Return nothing
    """
    # create list of args for ecami
    ecami_args = [
        "python3",
        "prediction.py",
        "-input",
        str(input_path),
        "-kmer_db",
        "CAZyme",
        "-output",
        f"{out_dir}/ecami_output.txt",
    ]

    # create log of eCAMI run
    with open(f"{out_dir}/ecami.log", "w+") as fh:
        process = subprocess.run(ecami_args, stdout=fh, text=True)

    # check if successful
    if process.returncode != 0:
        logger.warning(
            (
                f"eCAMI ran into error for {out_dir}\n."
                "eCAMI error:\n"
                f"{process.stderr}"
            )
        )

    return
