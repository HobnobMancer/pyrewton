# EastBIO PhD Project Scripts, packaged into pyrewton

[![DOI](https://zenodo.org/badge/243783792.svg)](https://zenodo.org/badge/latestdoi/243783792)
[![Funding](https://img.shields.io/badge/Funding-EASTBio-blue)](http://www.eastscotbiodtp.ac.uk/)
[![PhD licence](https://img.shields.io/badge/Licence-MIT-green)](https://github.com/HobnobMancer/PhD_Project_Scripts/blob/master/LICENSE)
[![CircleCI](https://circleci.com/gh/HobnobMancer/PhD_Project_Scripts.svg?style=shield)](https://circleci.com/gh/HobnobMancer/PhD_Project_Scripts)
[![codecov](https://codecov.io/gh/HobnobMancer/PhD_Project_Scripts/branch/master/graph/badge.svg)](https://codecov.io/gh/HobnobMancer/PhD_Project_Scripts)
[![Documentation Status](https://readthedocs.org/projects/phd-project-scripts/badge/?version=latest)](https://phd-project-scripts.readthedocs.io/en/latest/?badge=latest)
[![Python](https://img.shields.io/badge/Python-v3.7.---orange)](https://www.python.org/about/)
[![Research](https://img.shields.io/badge/Bioinformatics-Protein%20Engineering-ff69b4)](http://www.eastscotbiodtp.ac.uk/eastbio-student-cohort-2019)
[![Version](https://img.shields.io/badge/Version-v0.1.1-9cf)]]

_Please find more detailed documentation at for operation and troubleshooting at [Read the Docs](https://phd-project-scripts.readthedocs.io/en/latest/)_

## Contents

1. [Overview](#Overview)
2. [Installation](#Installation)
    - [Requirements](#Requirements)
3. [Directories](#Directories)
4. [Modules](#Modules)
    - [genbank](#genbank)
    - [cazymes](#cazymes)
        - [uniprot](#uniprot)
        - [prediction](#prediction)
        
## Overview

`pyrewton` is a Python3 script package for the automated identification of CAZyomes (call carbohydrates encoded within the genome of an given species). The package is run at the command line and free to use under the MIT license.

All modules, submodules and associated Python scripts are located within the `pyrewton` directory in this repository. Specifically, `pyrewton` supports:

- Downloading of all genomic assemblies (as GenBank files .gbff) from the [NCBI Assembly database](https://www.ncbi.nlm.nih.gov/assembly)
associated with each species passed to the programme
- Retrieval of all annotated protein sequences from GenBank (.gbff) files
- Retrieve proteins entries from [UniProtKB](https://www.uniprot.org/), using a JSON file to configure the queries

Features currently in development:
- Use the tools [dbCAN](https://github.com/linnabrown/run_dbcan), [CUPP](https://www.bioengineering.dtu.dk/english/researchny/research-sections/section-for-protein-chemistry-and-enzyme-technology/enzyme-technology/cupp), and [eCAMI](https://github.com/zhanglabNKU/eCAMI) to predict the which query protein sequences are CAZymes and predict their CAZy family
- Evaluate the accuracy of the CAZyme prediction tools to distinguish between CAZyme and non-CAZyme protein sequences
- Evaluate the accuracy of the CAZyme prediction tools to the correct CAZy family
- Produce a report of the CAZyme prediction tool evaluation

Inital plans and devleopment plans are stored within the [Wiki](https://github.com/HobnobMancer/PhD_Project_Scripts/wiki).

## Installation

The easiest way to install `pyrewton` is to use `pip3`:  
`pip3 install -e <path to directory containing pyrewton setup.py>`.  
Pass the path to the **directory** containing the setup.py file, **not** the path to the setup.py file.
 
Using this method of installation will install all required Python package dependencies.

## Requirements

Python version 3.7+
Miniconda3 managed microenvironment, incorporated code checkers are included in list form in 'requirements.txt'.
Miniconda3 environment file is also available in the GitHub repository: 'environment.yml'.
For all required Python libraries please read 'requirements.txt'.

## Directories

Below is a directory plan of this repositorty, followed by a brief overview of each directories role , to facilitate navigation through the repository.

### assets

Directory containing all files needed for the GitHub page, created for easy access to accompnaying Jupyter notebooks.

### docs

Directory containing files to build documentation hosted at ReadTheDocs.

### notebooks

Directory containing all jupyter notebooks, and html copies used for easier in-browser viewing via the GitHub pages.

### tests

Directory containing all `pytest` files for testing `pyrewton`, including subdirectories for test inputs and targets. Each module/submodule has its own specific test input and target subdirectory.

### pyrewton

Directory containing all `pyrewton` program modules (including all submodules and Python scripts).

## Modules

_Please find more detailed documentation at for operation and troubleshooting at [Read the Docs](https://phd-project-scripts.readthedocs.io/en/latest/)_

This is a brief overview of the functionalities of each module within `pyrewton`. For more detailed documentation on the operation of each module and indiviudal Python scripts please see the documentation at [Read the Docs](https://phd-project-scripts.readthedocs.io/en/latest/) or the module specific READMEs, found in their respecitve module repository.

### utilities

Contains all Python scripts for building command-line parsers and loggers.

### file_io

Contains functions for handling directories and files in `pyrewton`, including retrieving program inputs and creating output directories.

### genbank

Directory containing all submodules that are involved in retrieving handling GenBank files. This includes retrieval of GenBank files from GenBank, and retrieval of protein sequences from GenBank files.

### cazymes

Contains all submodules associated with identifying CAZymes, including:
    - `uniprot` which retrieves proteins from UniProt which meet query criteria set out in a configuration JSON file
    - `prediction` which predicates which query protein sequences are CAZymes and non-CAZymes, evaluate the accuracy of the tools and produce a report of this evaluation

## Repository renamed 2020-10-05
**Note:** This repository was renamed from 'PhD_Project_Scripts' to 'pyrewton'.