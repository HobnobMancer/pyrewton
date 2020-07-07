#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author:
# Emma E. M. Hobbs

# Contact
# eemh1@st-andrews.ac.uk

# Emma E. M. Hobbs,
# Biomolecular Sciences Building,
# University of St Andrews,
# North Haugh Campus,
# St Andrews,
# KY16 9ST
# Scotland,
# UK

# The MIT License

"""Tests build of loggers for pyrewton.

These tests are inteded to be run from the root repository using:
pytest -v
"""

from argparse import Namespace
from pathlib import Path

import pytest

from pyrewton.loggers import build_logger


@pytest.fixture
def logger_output():
    output = Path("tests/")
    return output


@pytest.fixture
def logger_args(logger_output):
    argsdict = {"args": Namespace(verbose=False, log=logger_output)}
    return argsdict


@pytest.mark.skip(reason="Target files do not exist in repository")
def test_build_logger(null_logger, logger_args):
    """Tests building of logger"""
    build_logger(null_logger, logger_args["args"])
