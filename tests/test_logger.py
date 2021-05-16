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

from pyrewton.utilities import build_logger, config_logger


@pytest.fixture
def logger_output(test_dir):
    logger_output = (
        test_dir / "test_targets" / "bld_lggr_test_targets"
    )
    return logger_output


# create fixture to test when args.verbose is false
@pytest.fixture
def logger_args_false(logger_output):
    argsdict = {"args": Namespace(verbose=False, log=(logger_output / "test_logger.log"))}
    return argsdict


# create fixture to test when args.verbose is true
@pytest.fixture
def logger_args_true(logger_output):
    argsdict = {"args": Namespace(verbose=True, log=(logger_output / "test_logger.log"))}
    return argsdict


def test_build_logger_v_false(logger_output):
    """Tests building of logger"""
    build_logger(logger_output, "test_logger.log")


def test_build_logger_v_true(logger_output):
    """Tests building of logger"""
    build_logger(None, "test_logger")


def test_config_logger_v_false(logger_args_false):
    """Tests building of logger"""
    config_logger(logger_args_false["args"])


def test_config_logger_v_true(logger_args_true):
    """Tests building of logger"""
    config_logger(logger_args_true["args"])
