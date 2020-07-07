#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import logging
import unittest

from argparse import Namespace
from pathlib import Path

import pandas as pd
import pytest

from pyrewton.file_io import make_output_directory, write_out_dataframe


class Test_housekeeping_functions(unittest.TestCase):

    """Class defining tests of output_dir_handling_main.py."""

    # Establish inputs for tests and expected outputs

    def setUp(self):
        """"Retrieve inputs and targets for tests."""

        # Define test directories
        self.test_dir = Path("tests")
        self.output_dir = self.test_dir / "test_targets" / "out_dir_hndng_test_targets"
        self.df_output = self.output_dir / "test_writing_df.csv"

        self.args = Namespace(output=self.output_dir, nodelete=False, force=True)

        # Null logger instance
        self.logger = logging.getLogger("Test_logger_parser_output")
        self.logger.addHandler(logging.NullHandler())

        # Create dataframe (df)
        self.df_data = [["A", "B", "C"]]
        self.df = pd.DataFrame(self.df_data, columns=["C1", "C2", "C3"])

    # Define function to test

    @pytest.mark.run(order=9)
    def test_output_dir_creation(self):
        """Tests function for creating output dir"""
        make_output_directory(self.args, self.logger)

    @pytest.mark.run(order=10)
    def test_writing_df(self):
        """Tests function for writing out created dataframe"""
        write_out_dataframe(self.df, self.logger, self.df_output, True, False)