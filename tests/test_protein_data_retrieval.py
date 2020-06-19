#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import logging
import pytest
import unittest

import pandas as pd

from argparse import Namespace
from pathlib import Path

from pyrewton.genbank.get_cazyme_annotations.get_cazyme_annotations import (
    get_df_foundation_data,
    get_uniprotkb_data,
)


class Test_parsing_input_file(unittest.TestCase):
    """Class defining tests of get_cazyme_annotations.py to build pandas dataframe."""

    # Establish inputs for tests and expected outputs

    def setUp(self):
        """"Set up null logger instance and input file path."""

        # Define test directories
        self.test_dir = Path("tests")
        self.input_dir = self.test_dir / "test_inputs" / "gt_czym_anntns_test_inputs"

        # Null logger instance
        self.logger = logging.getLogger("test_input_handling_logger")
        self.logger.addHandler(logging.NullHandler())

        # establish input paths
        self.input_df_path = self.input_dir / "test_input_df.csv"
        self.genbank_file_path = self.input_dir / "test_genbank_files"

        # Define Namespace for input dataframe and GenBank file paths
        self.argsdict = {"args": Namespace(genbank=self.genbank_file_path)}

        # Retrieve series from input dataframe
        self.input_df = pd.read_csv(
            self.input_df_path,
            header=0,
            names=["Genus", "Species", "NCBI Taxonomy ID", "NCBI Accession Numbers"],
        )
        self.df_series_0 = self.input_df.iloc[0]
        self.df_series_1 = self.input_df.iloc[1]

    # Define tests
    # Tests pass by being able to access and read input objects without incurring errors

    @pytest.mark.run(order=15)
    def test_genbank_protein_data_df(self):
        """Tests the function which builds dataframe from GenBank files protein data."""
        get_df_foundation_data(self.df_series_0, self.argsdict["args"], self.logger)

    @pytest.mark.run(order=16)
    def test_uniprot_result_handling(self):
        """Test function which retrieves, reads and processes input dataframe(df)."""
        get_uniprotkb_data(self.df_series_1, self.logger)
