# PhD_Project_Scripts

## Overview

PhD_Project_Scripts contains all scripts associated with the EASTBio Proteng project.

Completed scripts are packaged into `proteng`.

`proteng` is a Python3 script package that supports the pull down of genomic assemblies from the [NCBI Assembly database](https://www.ncbi.nlm.nih.gov/assembly).

[In GitHub](https://github.com/HobnobMancer/PhD_Project_Scripts) each section of the project is given its own folder and subfolders within.

## Installation

The easiest way to install `proteng` is to use `pip`:
`pip3 install proteng`.
This will also install all required Python packages dependencies.

### Requirements

Scripts are written using:
Python = 3.7
Miniconda3 managed microenvironment, incorporated code checkers are included in list form in 'requirements.txt'
Miniconda3 environment file is also available in the GitHub repository: 'environment.yml'.

## Scripts

A brief overview of each script.

### Extract_genomes_NCBI

Extract_genomes_NCBI takes a plain text file containing species scientific names or NCBI taxonomy as input. The script will find the corresponding taxonomy ID or scientific name, as appropriate, as well as retrieve all directly linked accession numbers, generating a dataframe of 'Genus', 'Species', 'NCBI Taxonomy ID' and 'NCBI Accession number'.
