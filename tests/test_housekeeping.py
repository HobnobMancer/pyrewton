#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import logging
import unittest

from argparse import Namespace
from pathlib import Path

import pandas as pd
import pytest

from Bio import Entrez

from pyrewton import directory_handling, loggers, parsers

Entrez.email = "my_email@my.domain"


class Test_housekeeping_functions(unittest.TestCase):

    """Class defining tests of get_ncbi_genomes.py housekeeping functions.

    These include creating the parser, the logger and output dir.
    """

    # Establish inputs for tests and expected outputs

    def setUp(self):
        """"Retrieve inputs and targets for tests."""

        # Define test directories
        self.test_dir = Path("tests")
        self.output_dir = self.test_dir / "test_targets" / "gt_ncbi_gnms_test_targets"
        self.df_output = self.output_dir / "gt_ncbi_gnms_test_df.csv"

        # Null logger instance
        self.logger = logging.getLogger("Test_logger_parser_output")
        self.logger.addHandler(logging.NullHandler())

        # Define test inputs
        self.test_logger = "test_logger"

        # Define Namespace and disable genbank download
        self.argsdict = {"args": Namespace(verbose=False, log=None)}

        # Create dataframe (df)
        self.df_data = [["A", "B", "C"]]
        self.df = pd.DataFrame(self.df_data, columns=["C1", "C2", "C3"])

    # Define function to test

    @pytest.mark.run(order=1)
    def test_build_parser(self):
        """Tests building of parser"""
        parsers.parser_get_ncbi_genomes.build_parser()
        parsers.parser_get_cazyme_annotations.build_parser()

    @pytest.mark.run(order=2)
    def test_build_logger(self):
        """Tests building of logger"""
        loggers.logger_pyrewton_main.build_logger(
            self.test_logger, self.argsdict["args"]
        )

    @pytest.mark.run(order=3)
    def test_output_dir_creation(self):
        """Tests function for creating output dir"""
        directory_handling.output_dir_handling_main.make_output_directory(
            self.output_dir, self.logger, True, True
        )

    @pytest.mark.run(order=4)
    def test_writing_df(self):
        """Tests function for writing out created dataframe"""
        directory_handling.output_dir_handling_main.write_out_dataframe(
            self.df, self.logger, self.df_output, True, False
        )
