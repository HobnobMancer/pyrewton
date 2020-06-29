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
"""Build loggers for pyrewton modules."""

import logging


def build_logger(script_name, args) -> logging.Logger:
    """Return a logger for this script.

    Enables logger for script, sets parameters and creates new file to store log.

    :param script_name: str, name of script
    :param args: parser argument

    Return logger object.
    """
    logger = logging.getLogger(script_name)

    # Set format of loglines
    log_formatter = logging.Formatter(
        script_name + ": {} - {}".format("%(asctime)s", "%(message)s")
    )

    # Setup console handler to log to terminal
    console_log_handler = logging.StreamHandler()
    if args.verbose is True:
        console_log_handler.setLevel(logging.INFO)
    else:
        console_log_handler.setLevel(logging.WARNING)
    console_log_handler.setFormatter(log_formatter)
    logger.addHandler(console_log_handler)

    # Setup file handler to log to a file
    if args.log is not None:
        file_log_handler = logging.FileHandler(args.log)
        if args.verbose is True:
            file_log_handler.setLevel(logging.INFO)
        else:
            file_log_handler.setLevel(logging.WARNING)
        file_log_handler.setFormatter(log_formatter)
        logger.addHandler(file_log_handler)

    return logger
