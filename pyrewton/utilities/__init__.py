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
"""Module for building loggers and cmd-line args parsers for pyrewton modules."""

import logging
import os

from pathlib import Path


def build_logger(output, file_name):
    """Build loggers with pre-defined parameters for writing out errors and failed scrapes.

    :param output: Path to output dir to write log file to or None
    :param file_name: str, name of output log file

    Return logger object.
    """
    logger = logging.getLogger(file_name[:-4])

    if output is None:
        output = os.getcwd()
        path_ = Path(f"{output}/{file_name}")
    else:
        path_ = output / f"{file_name}"

    # Set format of loglines
    log_formatter = logging.Formatter(file_name + ": {} - {}".format("%(asctime)s", "%(message)s"))

    # Setup file handler to log to a file
    file_log_handler = logging.FileHandler(path_)
    file_log_handler.setLevel(logging.WARNING)
    file_log_handler.setFormatter(log_formatter)
    logger.addHandler(file_log_handler)

    return logger
