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

import pytest

from pyrewton.loggers import build_logger


@pytest.mark.run(order=4)
def test_build_logger(null_logger, logger_args):
    """Tests building of logger"""
    build_logger("test_logger", logger_args["args"])
