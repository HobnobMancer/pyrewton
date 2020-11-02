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

import pytest

from pyrewton.utilities import build_logger


@pytest.fixture
def logger_output(test_dir):
    logger_output = (
        test_dir / "test_targets" / "bld_lggr_test_targets" / "test_bld_logger.log"
    )
    return logger_output


# create fixture to test when args.verbose is false
@pytest.fixture
def logger_args_false(logger_output):
    argsdict = {"args": Namespace(verbose=False, log=logger_output)}
    return argsdict


# create fixture to test when args.verbose is true
@pytest.fixture
def logger_args_true(logger_output):
    argsdict = {"args": Namespace(verbose=True, log=logger_output)}
    return argsdict


@pytest.mark.run(order=4)
def test_build_logger_v_false(null_logger, logger_args_false):
    """Tests building of logger"""
    build_logger("test_logger", logger_args_false["args"])


@pytest.mark.run(order=5)
def test_build_logger_v_true(null_logger, logger_args_true):
    """Tests building of logger"""
    build_logger("test_logger", logger_args_true["args"])
