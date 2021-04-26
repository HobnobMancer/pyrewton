#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author:
# Emma E. M. Hobbs
#
# Contact
# eemh1@st-andrews.ac.uk
#
# Emma E. M. Hobbs,
# Biomolecular Sciences Building,
# University of St Andrews,
# North Haugh Campus,
# St Andrews,
# KY16 9ST
# Scotland,
# UK
#
# The MIT License
"""Module to parse the output from dbCAN, CUPP and eCAMI. Create standarised output."""


import numpy as np


class CazymeDomain:
    """Single CAZyme domain in a protein, predicted by a CAZyme prediction tool.

    Each unique CAZy domain per protein is identifiable by a unique CAZy
    family-subfamily combination.

    Every CAZyme domain has a source CAZyme prediction tool that predicted the CAZyme
    domain, a parent CAZyme protein (represented by the protein accession), and CAZy
    family and subfamily combination. If no CAZy subfamily is predicted, the
    subfamily will be listed as a null value.

    Hotpep, CUPP and eCAMI predict EC numbers for each CAZyme domain.
    HMMER and CUPP predict amino acid domain ranges.
    Multiple EC numbers and domain ranges can be predicted for a single CAZyme domain,
    therefore, these attritbutes are stored as lists.
    """

    def __init__(
        self,
        prediction_tool,
        protein_accession,
        cazy_family,
        cazy_subfamily=None,
        ec_numbers=None,
        domain_range=None,
    ):
        """Initiate instance

        :attr prediction_tool: str, CAZyme prediciton tool which predicted the domain
        :attr protein_accession: str
        :attr cazy_family: str
        :attr cazy_subfamily: str
        :attr ec_numbers: list (list of str, each str contains a unique EC number)
        :attr domain_range: list (list of str, each str contains a unique domain range)
        """
        self.prediction_tool = prediction_tool
        self.protein_accession = protein_accession
        self.cazy_family = cazy_family

        # not all CAZyme domans are catalogued under a CAZy subfamily
        if cazy_subfamily is None:
            self.cazy_subfamily = np.nan
        else:
            self.cazy_subfamily = cazy_subfamily

        # EC numbers are not predicted for all CAZyme domains
        if ec_numbers is None:
            self.ec_numbers = []  # enables adding in EC#s included in another line of the file
        else:
            self.ec_numbers = ec_numbers

        # Not all prediction tools predict CAZyme domains
        if domain_range is None:
            self.domain_range = []  # enables adding domain range listed in another line of the file
        else:
            self.domain_range = domain_range

    def __str__(self):
        return(
            f"-CazymeDomain in {self.protein_accession}, "
            f"fam={self.cazy_family}, subfam={self.cazy_subfamily}-"
        )

    def __repr__(self):
        return(
            f"<CazymeDomain parent={self.protein_accession} fam={self.cazy_family}, "
            f"subfam={self.cazy_subfamily}>"
        )


class CazymeProteinPrediction:
    """Single protein and CAZyme/non-CAZyme prediction by a CAZyme prediction tool"""

    def __init__(
        self,
        prediction_tool,
        protein_accession,
        cazyme_classification,
        cazyme_domains=None,
    ):
        """Initate class instance.

        :attr prediction_tool: str, name of CAZyme prediction tool
        :attr protein_accession: str
        :attr cazyme_classification: int, 1=CAZyme, 0=non-CAZyme
        :attr cazyme_domains: list of CazymeDomain instances, domain predicted to be in the CAZyme
        """
        self.prediction_tool = prediction_tool
        self.protein_accession = protein_accession
        self.cazyme_classification = cazyme_classification  # CAZyme=1, non-CAZyme=0

        # non-CAZymes will have no cazyme_domains
        if cazyme_domains is None:
            self.cazyme_domains = []
        else:
            self.cazyme_domains = cazyme_domains

    def __str__(self):
        if self.cazyme_classification == 0:
            return f"-CazymeProteinPrediction, protein={self.protein_accession}, non-CAZyme-"
        else:
            return(
                f"-CazymeProteinPrediction, protein={self.protein_accession}, "
                f"CAZyme domains={len(self.cazyme_domains)}-"
            )

    def __repr__(self):
        return(
            f"<CazymeProteinPrediction, protein={self.protein_accession}, "
            f"cazyme_classification{self.cazyme_classification}>"
        )
