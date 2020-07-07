#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import logging
import unittest

import pytest

from pyrewton import parsers


class Test_housekeeping_functions(unittest.TestCase):

    """Class defining tests of get_ncbi_genomes.py housekeeping functions.

    These include creating the parser, the logger and output dir.
    """

    # Establish inputs for tests and expected outputs

    def setUp(self):
        """"Retrieve inputs and targets for tests."""

        # Null logger instance
        self.logger = logging.getLogger("Test_parser_output")
        self.logger.addHandler(logging.NullHandler())

        # Define test inputs
        self.test_logger = "test_logger"

    # Define function to test

    @pytest.mark.run(order=1)
    def test_build_parser_gt_ncb_gnms(self):
        """Tests building of parser"""
        parsers.parser_get_ncbi_genomes.build_parser()

    @pytest.mark.run(order=2)
    def test_build_parser_gt_unprt_prtns(self):
        """Tests building of parser"""
        parsers.parser_get_uniprot_proteins.build_parser()

    @pytest.mark.run(order=3)
    def test_build_parser_gt_czym_anno(self):
        """Tests building of parser"""
        parsers.parser_get_cazyme_annotations.build_parser()
