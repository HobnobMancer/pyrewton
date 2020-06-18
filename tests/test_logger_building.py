#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import logging
import unittest

from argparse import Namespace
from pathlib import Path

import pytest

from pyrewton import loggers


class Test_housekeeping_functions(unittest.TestCase):

    """Class defining tests of get_ncbi_genomes.py housekeeping functions.

    These include creating the parser, the logger and output dir.
    """

    # Establish inputs for tests and expected outputs

    def setUp(self):
        """"Retrieve inputs and targets for tests."""

        # Define test directories
        self.test_dir = Path("tests")
        self.output_dir = self.test_dir / "test_targets" / "bld_lggr_test_targets"
        self.log_output = self.output_dir / "test_bld_logger.log"

        # Null logger instance
        self.logger = logging.getLogger("Test_logger_build")
        self.logger.addHandler(logging.NullHandler())

        # Define test inputs
        self.test_logger = "test_logger"

        # Define Namespace and disable genbank download
        self.argsdict = {"args": Namespace(verbose=False, log=self.log_output)}

    # Define function to test

    @pytest.mark.run(order=3)
    def test_build_logger(self):
        """Tests building of logger"""
        loggers.logger_pyrewton_main.build_logger(
            self.test_logger, self.argsdict["args"]
        )
