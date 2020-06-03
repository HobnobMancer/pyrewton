#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import logging
import unittest

from argparse import Namespace
from pathlib import Path

import pandas as pd
import pytest

from Bio import Entrez

from Section1_Extracting_Genomes import Extract_genomes_NCBI

Entrez.email = "proteng.ext_gnm_ncbi@my.domain"


class Test_housekeeping_functions(unittest.TestCase):

    """Class defining tests of Extract_genomes_NCBI.py housekeeping functions.

    These include creating the parser, the logger and output dir.
    """

    # Establish inputs for tests and expected outputs

    def setUp(self):
        """"Retrieve inputs and targets for tests."""

        # Define test directories
        self.test_dir = Path("tests")
        self.output_dir = (
            self.test_dir / "test_targets" / "test_ext_gnm_ncbi" / "EgN_target_dir"
        )
        self.df_output = self.output_dir / "test_df.csv"

        # Null logger instance
        self.logger = logging.getLogger("Test_name_and_ID_Retrieval logger")
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
        Extract_genomes_NCBI.build_parser()

    @pytest.mark.run(order=2)
    def test_build_logger(self):
        "Tests building of logger" ""
        Extract_genomes_NCBI.build_logger(self.test_logger, self.argsdict["args"])

    @pytest.mark.run(order=3)
    def test_output_dir_creation(self):
        """Tests function for creating output dir"""
        Extract_genomes_NCBI.make_output_directory(
            self.output_dir, self.logger, True, True
        )

    @pytest.mark.run(order=4)
    def test_writing_df(self):
        """Tests function for writing out created dataframe"""
        Extract_genomes_NCBI.write_out_dataframe(
            self.df, self.logger, self.df_output, True, False
        )
