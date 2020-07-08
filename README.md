# EastBIO PhD Project Scripts, packaged into pyrewton

[![DOI](https://zenodo.org/badge/243783792.svg)](https://zenodo.org/badge/latestdoi/243783792)
[![Funding](https://img.shields.io/badge/Funding-EASTBio-blue)](http://www.eastscotbiodtp.ac.uk/)
[![PhD licence](https://img.shields.io/badge/Licence-MIT-green)](https://github.com/HobnobMancer/PhD_Project_Scripts/blob/master/LICENSE)
[![CircleCI](https://circleci.com/gh/HobnobMancer/PhD_Project_Scripts.svg?style=shield)](https://circleci.com/gh/HobnobMancer/PhD_Project_Scripts)
[![codecov](https://codecov.io/gh/HobnobMancer/PhD_Project_Scripts/branch/master/graph/badge.svg)](https://codecov.io/gh/HobnobMancer/PhD_Project_Scripts)
[![Python](https://img.shields.io/badge/Python-v3.7.---orange)](https://www.python.org/about/)
[![Research](https://img.shields.io/badge/Bioinformatics-Protein%20Engineering-ff69b4)](http://www.eastscotbiodtp.ac.uk/eastbio-student-cohort-2019)

Version v0.1.1 2020/06/04

This repository contains all scripts associated with the EastBio PhD project ‘Identifying Engineering Candidates for Advanced Biocatalysis in Biofuel Production'.

## Contents

1. [Overview](#Overview)
2. [Installation](#Installation)
    - [Requirements](#Requirements)
3. [Directories, modules and files](#Directories)

## Overview

Scripts are packaged into `pyrewton`, which is a Python3 script package run at the command line, and free to use under the MIT license. All modules, submodules and associated Python scripts are located within the `pyrewton` directory in this repository. Specifically, `pyrewton` supports:

- Downloading of all genomic assemblies (as GenBank files .gbff) from the [NCBI Assembly database](https://www.ncbi.nlm.nih.gov/assembly)
associated with each species passed to the programme
- Retrieving cazyme annotations (and associated data from UniProtKB) from GenBank files
- Retrieve all cazyme entries in UniProt for a given species

Inital plans and devleopment plans are stored within the [Wiki](https://github.com/HobnobMancer/PhD_Project_Scripts/wiki).

## Installation

The easiest way to install `pyrewton` is to use `pip`:  
`pip3 install -e <path to directory containing pyrewton setup.py>`.  
Pass a path to the **directory** containing the setup.py file, **not** the path to the setup.py file.  
Using this method of installation will install all required Python package dependencies.

## Requirements

Python version 3.7+
Miniconda3 managed microenvironment, incorporated code checkers are included in list form in 'requirements.txt'.
Miniconda3 environment file is also available in the GitHub repository: 'environment.yml'.

## Directories

Below is a directory plan of the pyrewton module structure, followed by a brief overview of each directories role in the repository, to facilitate navigation through the repository.

### assets

Directory containing all files needed for the GitHub page.

### docs

Directory containing files to build documentation hosted at ReadTheDocs.

### notebooks

Directory containing all jupyter notebooks, and html copies used for easier in-browser viewing via the GitHub pages.

### tests

Directory containing all `pytest` files for testing `pyrewton` during development, including subdirectories for test inputs and targets, with each module/submodule possessing its own specific test input and target subdirectory.

### pyrewton

Directory containing all `pyrewton` program modules (including all submodules and Python scripts).

### pyrewton module: parsers

Directory containing all Python scripts for building command-line parsers.

### pyrewton module: loggers

Directory containing Python scripts for building loggers.

### pyrewton module: file_io

Directory contains functions for handling directories and files in `pyrewton` Python scripts, including retrieving program inputs and creating output directories.

### pyrewton module: genbank

Directory containing all submodules that are involved in handling GenBank files.

**Genbank submodule: get_ncbi_genomes** - Takes a list of species and downloads all directly linked GenBank (.gbff) files in the NCBI Assembly database.

**Genbank submodule: get_genbank_annotations** - Retrieve all protein annotations from GenBank files.

### pyrewton module: annotations

Directory containing all submodules that are involved in retrieving genomic annotations.

**get_uniprot_proteins.py** - retrieves all protein entries in UniProt in the species passed to the script.  
**search_uniprot_proteins.py** - search proteins retrieved from UniProt to identify those with link to the CAZy database, and those with potential cazyme functionality inferred from their EC number and/or GO (Gene Ontology) annotated function.  
**search_genbank_annotations.py** - search dataframe of cazymes retrieved from UniProt against dataframe of proteins retrieved from GenBank files to identify cazymes retrieved from the GenBank files. Then search GenBank proteins for those with inferred cazyme functionality and write these out to another dataframe.

