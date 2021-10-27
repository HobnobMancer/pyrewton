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
"""Funcs for handling file I/O when creating tests sets for the CAZyme prediction tool evaluation"""


import logging
import json
import shutil
import sys
import yaml

from pathlib import Path


def retrieve_assemblies_dict(yaml_path):
    """Retrieve dict of genomic assemblies to create test sets from.

    :param yaml_path: Path to yaml file, keyed by NCBI:txid and valued by genomic assemblies

    Return dict, keyed by NCBI:txid and valued by list genomic assemblies"""
    logger = logging.getLogger(__name__)
    # open the YAML file
    try:
        with open(yaml_path, "r") as fh:
            assemblies_dict = yaml.full_load(fh)
    except FileNotFoundError:
        logger.error(
            "Did not find the configuration file. Check the path is correct.\n"
            "Terminating programme"
        )
        sys.exit(1)
    return assemblies_dict


def get_cazy_dict(cazy_path):
    logger = logging.getLogger(__name__)
    try:
        with open(cazy_path, "r") as fh:
            cazy_dict = json.load(fh)
    except FileNotFoundError:
        logger.error(
            "Did not find the configuration file. Check the path is correct.\n"
            "Terminating programme"
        )
        sys.exit(1)
    return cazy_dict


def prepare_output_dir(output_dir):
    """Delete and make temporary output directory"""
    output_dir = Path(output_dir)
    try:
        if output_dir.exists:
            shutil.rmtree(output_dir)
    except (FileExistsError, FileNotFoundError) as e:
        pass
    try:
        output_dir.mkdir(exist_ok=True)
    except FileExistsError:
        print("output_dir_exists")
    return
