#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import logging
import pytest
import unittest

from argparse import Namespace
from pathlib import Path

from pyrewton.directory_handling.input_dir_get_cazyme_annotations import (
    check_input,
    get_input_df,
    get_genbank_file,
)


class Test_parsing_input_file(unittest.TestCase):

    """Class defining tests of get_cazyme_annotations.py ability to parse
    provided input"""

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

        self.argsdict = {
            "args": Namespace(
                df_input=self.input_df_path, genbank=self.genbank_file_path
            )
        }

    # Define tests

    @pytest.mark.run(order=6)
    def test_input_validity_checking(self):
        """Tests the function which checks both inputs are reachable/available."""
        check_input(self.argsdict["args"], self.logger)

    @pytest.mark.run(order=7)
    def test_input_dataframe_retrieval(self):
        """Test function which retrieves, reads and processes input dataframe(df)."""
        get_input_df(self.input_df_path, self.logger)

    @pytest.mark.run(order=8)
    def test_genbank_reading(self):
        """Tests ability to read and process GenBank file."""
        get_genbank_file("GCA_011074995.1", self.argsdict["args"], self.logger)
