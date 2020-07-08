#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Author : Emma E. M. Hobbs
#
# Contact:
# eemh1@st-andrews.ac.uk
#
# Emma Hobbs,
# School of Biology,
# University of St Andrews,
# Biomedical Sciences Research Complex,
# St Andrews,
# Fife,
# KY16 9ST
# Scotland,
# UK
#
# MIT License

import setuptools

from pathlib import Path


# get long description from README.md
with Path("README.md").open("r") as long_description_handle:
    long_description = long_description_handle.read()


setuptools.setup(
    name="pyrewton",
    version="0.1.1",
    # Metadata
    author="Emma E. M. Hobbs",
    author_email="eemh1@st-andrews.ac.uk",
    description="".join(
        [
            (
                "pyrewton provides the scripts for the EASTBio protein "
                "engineering PhD project, supporting the identification "
                "of engineer candidates for catalysis in biofuel production"
            )
        ]
    ),
    long_description=long_description,
    long_description_content_type="text/markdown",
    license="MIT",
    keywords="genome bioinforamtics protein engineering",
    platforms="Posix, MacOS X",
    url="https://github.com/HobnobMancer/PhD_Project_Scripts",  # Github repository
    entry_points={
        "console_scripts": [
            "get_ncbi_genomes.py = pyrewton.genbank.get_ncbi_genomes.get_ncbi_genomes:main",
            "get_genbank_annotations.py = pyrewton.genbank.get_genbank_annotations."
            "get_genbank_annotations:main",
            "get_uniprot_proteins.py = pyrewton.annotations.get_uniprot_proteins:main",
            "search_uniprot_proteins.py = pyrewton.annotations.search_uniprot_proteins:main",
            "search_genbank_annotations.py = pyrewton.annotations.search_genbank_annotations:main",
        ]
    },
    # Ensure all additional requirements are installed
    install_requires=["biopython>=1.76", "pandas>=1.0.3"],
    # Include conda microenvironment
    # and template input file for Extract_genomes_NCBI.py
    package_data={
        "Conda microenvironment": ["environment.yml"],
        "Get NCBI genomes input file": ["get_ncbi_genomes_template_input_file.txt"],
    },
    include_package_data=True,
    classifiers=[
        "Environment :: Console",
        "Intended Audience :: Science/Research",
        "Licence :: OSI Approved :: MIT License",
        "Natural Language :: English",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: POSIX :: Linux",
        "Programming Language :: Python :: 3.7",
        "Topic :: Scientific/Engineering :: Bioinformatics",
    ],
)
