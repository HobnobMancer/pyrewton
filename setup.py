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

import os
import setuptools
import subprocess

from setuptools import Command

from pathlib import Path


class InstallCPTs(Command):
    """A custom command to install CAZyme prediction tools (CPTs): dbCAN, eCAMI and CUPP."""

    description = "Install dbCAN, eCAMI and CUPP"
    user_options = [('pyrewton-dir=', 'p', 'path to dir containing pyrewton setup.py file')]

    def initialize_options(self):
        self.pyrewton_dir = None

    def finalize_options(self):
        if self.pyrewton_dir is None:
            self.pyrewton_dir = os.path.dirname(os.path.abspath(__file__))  # get abspath to dir of setup.py
        elif os.path.isdir(self.pyrewton_dir):
            self.pyrewton_dir = os.path.dirname(os.path.abspath(__file__))  # get abspath to dir of setup.py
        elif self.pyrewton_dir == ".":
            self.pyrewton_dir = os.path.dirname(os.path.abspath(__file__))  # get abspath to dir of setup.py

    def run(self):
        """Run command"""
        installation_path = self.pyrewton_dir + "/installation.sh"
        tools_dir = self.pyrewton_dir + "/pyrewton/cazymes/prediction/tools"
        subprocess.check_call([installation_path, tools_dir])


# get long description from README.md
if __file__ == "setup.py":
    with Path("README.md").open("r") as long_description_handle:
        long_description = long_description_handle.read()

else:
    path = __file__
    path = path.replace("setup.py", "")
    path = Path(path) / "README.md"
    with path.open("r") as long_description_handle:
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
            "get_uniprot_proteins.py = pyrewton.cazymes.uniprot.get_uniprot_proteins:main",
            "predict_cazymes.py = pyrewton.cazymes.prediction.predict_cazymes:main",
        ]
    },
    # Ensure all additional requirements are installed
    install_requires=[
        "biopython>=1.76",
        "bioservices>=1.7.9",
        "numpy>=1.19.4",
        "pandas>=1.0.3",
        "pyyaml>=5.3.1",
        "run-dbcan==2.0.11",
        "scipy>=1.5.4",
        "tqdm>=4.53.0",
    ],
    # Include conda microenvironment
    # and template input file for Extract_genomes_NCBI.py
    package_data={
        "Conda microenvironment": ["environment.yml"],
        "get_ncbi_genomes input file": ["get_ncbi_genomes_template_input_file.txt"],
        "get_uniprot_proteins config file": ["uniprot_config.yaml"],
        "get_uniprot_proteins ec numbers config file": ["uniprot_cazy_ec_retrieval.yaml"],
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
    cmdclass={'cpt': InstallCPTs},
)
