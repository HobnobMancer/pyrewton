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

"""Tests for building of parsers for pyrewton.

These tests are inteded to be run from the root repository using:
pytest -v
"""


import pytest

from pyrewton.utilities.parsers import (
    cmd_parser_calc_stats,
    cmd_parser_get_evaluation_dataset_from_dict,
    cmd_parser_get_evaluation_dataset,
    cmd_parser_get_genbank_annotations,
    cmd_parser_get_ncbi_genomes,
    cmd_parser_get_uniprot_proteins,
    cmd_parser_predict_cazymes,
)


def test_parser_gt_ncb_gnms():
    """Tests building of parser"""
    cmd_parser_get_ncbi_genomes.build_parser()


def test_parser_gt_ncb_gnms_argc():
    """Tests building of parser for get_ncbi_genomes when argv is not None"""
    cmd_parser_get_ncbi_genomes.build_parser(["dummy_email"])


def test_parser_gt_gnbnk_anno():
    """Tests building of parser"""
    cmd_parser_get_genbank_annotations.build_parser()


def test_parser_gt_gnbnk_anno_argv():
    """Tests building of parser for get_genbank_annotations when argv is not None"""
    cmd_parser_get_genbank_annotations.build_parser(["tests", "tests"])


def test_parser_gt_unprt_prtns():
    """Tests building of parser"""
    cmd_parser_get_uniprot_proteins.build_parser()


def test_parser_gt_unprt_prtns_argv():
    """Tests building of parser for get_uniprot_proteins when argv is not None"""
    cmd_parser_get_uniprot_proteins.build_parser(["tests/"])


def test_parser_prdct_czyms():
    """Test building parser for predict_cazymes.py"""
    cmd_parser_predict_cazymes.build_parser()


def test_parser_prdct_czyms_argv():
    """Test building parser for predict_cazymes.py when argv is not none"""
    cmd_parser_predict_cazymes.build_parser(["tests/"])


def test_parser_calc_stat():
    """Test building parser for calculate_stats.py"""
    cmd_parser_calc_stats.build_parser()


def test_parser_calc_stat_argv():
    """Test building parser for calculate_stats.py when argv is not none"""
    cmd_parser_calc_stats.build_parser(["tests/", "cazy_dict"])


def test_parser_get_eval_dict():
    """Test building parser for building evaluation test sets from a dict"""
    cmd_parser_get_evaluation_dataset_from_dict.build_parser()


def test_parser_get_eval_dict_argv():
    """Test building parser for building evaluation test sets from a dict, when argv is not none"""
    cmd_parser_get_evaluation_dataset_from_dict.build_parser(["tests/", "yaml", "cazy", "output"])


def test_parser_eval_testsets():
    """Test building parser for building evaluation test sets from a db"""
    cmd_parser_calc_stats.build_parser()


def test_parser_eval_testsets_argv():
    """Test building parser for building evaluation test sets from a db, when argv is not none"""
    cmd_parser_calc_stats.build_parser(["tests/", "test_dict"])
