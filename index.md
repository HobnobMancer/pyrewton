# Pyrewton
## Independent, comprehensive benchmarking of CAZyme classifiers and Creation of a CAZome SQL database


[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3876218.svg)](https://doi.org/10.5281/zenodo.3876218)
[![Funding](https://img.shields.io/badge/Funding-EASTBio-blue)](http://www.eastscotbiodtp.ac.uk/)
[![PhD licence](https://img.shields.io/badge/Licence-MIT-green)](https://github.com/HobnobMancer/PhD_Project_Scripts/blob/master/LICENSE)
[![CircleCI](https://circleci.com/gh/HobnobMancer/pyrewton.svg?style=shield)](https://circleci.com/gh/HobnobMancer/PhD_Project_Scripts)
[![codecov](https://codecov.io/gh/HobnobMancer/pyrewton/branch/master/graph/badge.svg)](https://codecov.io/gh/HobnobMancer/pyrewton)
[![Documentation Status](https://readthedocs.org/projects/pyrewton/badge/?version=latest)](https://pyrewton.readthedocs.io/en/latest/?badge=latest)
[![Python](https://img.shields.io/badge/Python-v3.7.---orange)](https://www.python.org/about/)
[![Research](https://img.shields.io/badge/Bioinformatics-Protein%20Engineering-ff69b4)](http://www.eastscotbiodtp.ac.uk/eastbio-student-cohort-2019)

This page is published as supplementary material to the main [`pyrewton` repo](), and the independet evaluation of the CAZyme classifiers [dbCAN](), [CUPP](), and [eCAMI]() presented in [Hobbs et al., 2021]().

This page covers includes the specific method used to perform the presented evaluation of the CAZyme classifiers listed above, including additional supplemenatry material to faciltiate the reporduction of the analysis.

Information for the general operation of `pyrewton` and support troubleshooting can be found at [Read the Docs](https://phd-project-scripts.readthedocs.io/en/latest/).

## Table of Contents

1. [Files and Directories](#linkfiles)
2. [Quick Start](#linkquick)
3. [Notebooks](#linkuse)
4. [Help](#linkhelp)

### Files and Directories<a id="linkfile"><a/>
  
The `pyrewton` repository is separated into documentation and program modules.

- `assets/`: directory containing all files needed for the GitHub page.
- `docs/`: directory containing all files used for build documentation.
- `notebooks/`: directory containing all jupyter notebooks. The notebooks present the operation of some scripts, providing an explanation of the specific computational processes of some scripts.
- `tests/`: directory containing all unit-tests, and required test inputs and targets.
- `supplementary/`: directory containing supplemenatry data, accompanying the presentation of the results in _Hobbs et al., 2021_.
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

  
### Installation

1. Create a new virtual environment  
  _To install Conda please see the [Conda documentation](https://docs.conda.io/en/latest/)._
  ```bash
  conda create -n pyrewton
  ```
2. Clone the `pyrewton` repository
  ```bash
  git clone https://github.com/HobnobMancer/pyrewton.git
  cd pyrewton  # navigate into the repo root
  ```
3. Install `pyrewton`
  ```bash
  pip3 install -e .
  ```
4. Install the CAZyme classifiers into the correct directories
  ```bash
  python3 . cpt -p .
  ```
 
Following this method ensure all requirments are installed, as well as installing the CAZyme classifers into the correct directorys.
  
### Requirements

POISx or Mac OS, or linux emulator
Python version 3.7+
Miniconda3 or Anaconda managed microenvironment, incorporated code checkers are included in list form in 'requirements.txt'.
Miniconda3 environment file is also available in the GitHub repository: 'environment.yml'.
For all required Python libraries please read 'requirements.txt'.
  
## Independent evaluation of CAZyme classifiers
 
At the time of publishing `pyrewton` no independent evaluation of the widely used CAZyme classifiers dbCAN, CUPP and eCAMI had been performed. Additionally, previous evaluations had not included an evaluation of the performance of the tools to differentiate between CAZymes and non-CAZymes predict the CAZy class annotation and had not evaluated the multi-label CAZy family annotation performance of the classifiers.
  
### Creation of the test sets
  
The Python script `` was used to generate the test sets, and was invoked using the following command:
```bash

```
  
The yaml file containing the genomic accessions and taxonomy IDs of the genomes selected for the inclusion in the study is lcoated in [`supplementary/test_set_data`]().
  
In total 70 test sets were created for the evaluation, 39 from Bacteria genomes and 31 from Eukaryote genomes.
  
### Invokving the classifiers

The Python script `` was used to invoke each CAZyme classifier (dbCAN, CUPP and eCAMI) for each test set created in the last step. The command listed below was used and should be run from ``. This is because dbCAN, CUPP and eCAMI use hard coded path to find their respective data, which requires `pyrewton` to navigate to the correct directors for the CAZyme classifers to function properly and this cannot be achieved if the script `` is not invoked in the correct directory.

```bash
  
```

### Calculating statistics

To statistical evaluate the performance of the CAZyme classifiers, the Python script `` was invoked using the following command:
```bash
```

### Presenting the findings

The R notebook ``, located in `` was used to generate visualsiations of the statistical output from ``.
 
To reuse this R notebook, the hard coded file import paths on lines .... need to be changed to match the location of the new `` ouput.

## Creation of a local CAZome database
 
  
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
