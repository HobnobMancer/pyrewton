# EastBIO PhD Project Supplementary Material

## Identifying, Characterising & Engineering Plant Cell Wall Degrading Enzymes for Enhanced Biocatalysts in Biofuel Production

[![DOI](https://zenodo.org/badge/243783792.svg)](https://zenodo.org/badge/latestdoi/243783792)
[![Funding](https://img.shields.io/badge/Funding-EASTBio-blue)](http://www.eastscotbiodtp.ac.uk/)
[![PhD licence](https://img.shields.io/badge/Licence-MIT-green)](https://github.com/HobnobMancer/PhD_Project_Scripts/blob/master/LICENSE)
[![CircleCI](https://circleci.com/gh/HobnobMancer/PhD_Project_Scripts.svg?style=shield)](https://circleci.com/gh/HobnobMancer/PhD_Project_Scripts)
[![codecov](https://codecov.io/gh/HobnobMancer/PhD_Project_Scripts/branch/master/graph/badge.svg)](https://codecov.io/gh/HobnobMancer/PhD_Project_Scripts)
[![Documentation Status](https://readthedocs.org/projects/phd-project-scripts/badge/?version=latest)](https://phd-project-scripts.readthedocs.io/en/latest/?badge=latest)
[![Python](https://img.shields.io/badge/Python-v3.7.---orange)](https://www.python.org/about/)
[![Research](https://img.shields.io/badge/Bioinformatics-Protein%20Engineering-ff69b4)](http://www.eastscotbiodtp.ac.uk/eastbio-student-cohort-2019)

These pages are published from the repository containing files used during the EastBIO PhD Project, for the identification of engineering candidates to enhance biocatalysis for biofuel production.

The `Jupyter` notebook environments have been used to facilitate explanation of the function and operation of the command-line python programmes used. These notebooks include explanatory text and snippets of code extract to illustrate code function, and thus not all code snippets are runnable in the `Jupyter` environment.

Please find more detailed information for operation and troubleshooting pyrewton (the Python 3 program written to perform the compuational work during the PhD project) at [Read the Docs](https://phd-project-scripts.readthedocs.io/en/latest/).

## Table of Contents

1. [Files and Directories](#linkfiles)
2. [Quick Start](#linkquick)
3. [Notebooks](#linkuse)
4. [Help](#linkhelp)

### Files and Directories<a id="linkfile"><a/>
  
The repository is separated into documentation and program modules.

Inital plans and devleopment plans are stored within the [Wiki](https://github.com/HobnobMancer/PhD_Project_Scripts/wiki).

- `assets/`: directory containing all files needed for the GitHub page.
- `docs/`: directory containing all files used for build documentation.
- `notebooks/`: directory containing all jupyter notebooks.
- `tests/`: directory containing all unit-tests, and required test inputs and targets.
- `supplementary/`: directory containing supplemenatry data, accompanying the PhD thesis.
- `pyrewton/`: directory containing all `pyrewton` program modules (including all submodules and Python scripts).
  - `parsers/`: directory containing all Python scripts for building parsers used in `pyrewton`.
  - `loggers/`: directory containing all Python scripts for building loggers in `pyrewton`.
  - `file_io/`: directory containing all Python scripts for handling directories and files in `pyrewton`.
  - `genbank/`: directory containing all submodules that are involved in retrieving NCBI GenBank files, and retrieving data from GenBank files.
    - `genbank/get_genomes_ncbi/`: directory containing the submodule `get_genomes_ncbi`, which takes a list of species (in a plain text file) and downloads all directly linked GenBank (.gbff) files from the NCBI Assembly database.
    - `genbank/get_genbank_annotations/`: directory containing the submodule `get_genbank_annoations`, which retrieves all protein annotations from GenBank (.gbff) files.
  - `annotations/`: directory containing Python scripts that make the up `annotations` module of `pyrewton`, which are involved in identifying cazyme annotations.
    - `get_uniprot_proteins.py`: script that retrieves all protein entries from UniProt for each species containing in a dataframe passed to the script.
    - `search_uniprot_proteins.py`: script that searches the proteins retrieved from UniProt for those with indicated cazyme functionality
    - `search_genbank_annotations.py`: script that searches retrieved GenBank annotations for those with indicated cazyme functionality.


### Quick Start<a id="linkquick"><a/>

The quickest way to install the programme `pyrewton` is to use pip: `pip3 install -e <path to DIRECTORY containing setup.py file>`. This will install all required Python packages and dependencies.<br/>
<font color="red"><b>Note:</b></font> [Conda](https://docs.conda.io/en/latest/) will need to be installed on your system.

### Notebooks <a id="linkuse"><a/>

Jupyter notebooks are accessible here, but also via the terminal if the programme `pyrewton` is installed on your system.

**Access the notebooks in browser:**

- [01_downloading_genbank_files](https://hobnobmancer.github.io/PhD_Project_Scripts/notebooks/01_downloading_genbank_files.html)<br/>
Notebook for invoking `get_ncbi_genomes.py`, which explains the function and operation of the Python script for downloading GenBank files from the NCBI Assembly database.
- [02_retrieving_genbank_annotation_summary](https://hobnobmancer.github.io/PhD_Project_Scripts/notebooks/02_retrieving_genbank_annotations.html)<br/>
Notebook for `get_genbank_annotations` for the retrieval of all protein annotations from the downloaded genomic assemblies.

**Access the notebooks via the terminal:**<br/>
From the top level of the repository directory on your system, start the `Jupyter` notebook server by issuing the command:
`Jupyter notebook`
This will opena  new browser window or tab containing the Jupyter homepage, with a listing of all files and directories as laid out in the repository.<br/>
Navigate to the 'Notebook/' directory, containing all Jupyter notebooks for the repository. To open a notebook simple click on the title of the notebook.

For more information please see the Jupyter notebook [Quick-start guide](https://jupyter-notebook-beginner-guide.readthedocs.io/en/latest/?fbclid=IwAR1yIwkYCDjcw5FJZ7CfKES3l72HubqGYGcFrVrUKwWZoYh4NHy3VVu0AgQ) and [tutorial](https://www.tutorialspoint.com/jupyter/jupyter_quick_guide.htm).

Each notebook is named to match the number of the section of the programme to which it is associated, such that the first notebook (01_Extracting_GenBank_files) refers to Section1_Extracting_genomes within the repository, as made explicitly clear in the notebook itself.
Again note that these notebooks contain exerts of code to facilitate the explanation of the code architecture and function, as well as details provided on how the code was implemented during the project. Therefore, the code in many of the code cells is not runnable and replication or exploration of the data should be performed using the original scripts provided within the repository.

If you are new or have little experience of using command-line programmes and/or the Python coding language, please read the appropriate notebooks  before using the associated scripts and posting queries/issues on the GitHub pages. The notebooks have been written so that those of all experience levels should be able to understand how and why the script operates.

### Help<a id="linkhelp"><a/>

Many of the common errors expected to arise during the operation of the scripts provided in this repository are named in the repository README, including the probable causes of these issues.

Please raise any issues with any of the programmes at the GitHub issues pages for this repository, by following [the link](https://github.com/HobnobMancer/PhD_Project_Scripts/issues).
