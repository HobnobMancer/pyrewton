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
"""Configuration file for pytest files.

Contains fixtures used by multiple test files.
"""

import logging

from pathlib import Path

import pytest


@pytest.fixture
def test_dir():
    return Path("tests/")


@pytest.fixture
def null_logger():
    logger = logging.getLogger("Test_ac_number_retrieval logger")
    logger.addHandler(logging.NullHandler())
    return logger
