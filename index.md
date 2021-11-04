# Pyrewton
## Independent, comprehensive benchmarking of CAZyme classifiers & Creation of a CAZome SQL database


[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3876218.svg)](https://doi.org/10.5281/zenodo.3876218)
[![Funding](https://img.shields.io/badge/Funding-EASTBio-blue)](http://www.eastscotbiodtp.ac.uk/)
[![PhD licence](https://img.shields.io/badge/Licence-MIT-green)](https://github.com/HobnobMancer/PhD_Project_Scripts/blob/master/LICENSE)
[![CircleCI](https://circleci.com/gh/HobnobMancer/pyrewton.svg?style=shield)](https://circleci.com/gh/HobnobMancer/PhD_Project_Scripts)
[![codecov](https://codecov.io/gh/HobnobMancer/pyrewton/branch/master/graph/badge.svg)](https://codecov.io/gh/HobnobMancer/pyrewton)
[![Documentation Status](https://readthedocs.org/projects/pyrewton/badge/?version=latest)](https://pyrewton.readthedocs.io/en/latest/?badge=latest)
[![Python](https://img.shields.io/badge/Python-v3.7.---orange)](https://www.python.org/about/)
[![Research](https://img.shields.io/badge/Bioinformatics-Protein%20Engineering-ff69b4)](http://www.eastscotbiodtp.ac.uk/eastbio-student-cohort-2019)

This page is published as supplementary material to the main [`pyrewton` repo](https://github.com/HobnobMancer/pyrewton), and the independet evaluation of the CAZyme classifiers [dbCAN](zhange _et al._, 2018), [CUPP](), and [eCAMI]() presented in [Hobbs et al., 2021]().

> Hobbs, E. E. M., Gloster, T. M., Chapman, S., Pritchard, L. (2021): Microbiology Society Annual Conference 2021. figshare. Poster. https://doi.org/10.6084/m9.figshare.14370836.v
> Zhang, H., Yohe, T., Huang, L., Entwistle, S., Wu, P., Yang, Z., Busk, P.K., Xu, Y., Yin, Y. (2018) 'dbCAN2: a meta server for automated carbohydrate-active enzyme annotation', _Nucleic Acids Res._, 46(W1), pp. W95-W101. doi: 10.1093/nar/gky418
> Barrett, K., Lange, L. (2019) 'Peptide-based functional annotation of carbohydrate-active enzymes by conserved unique peptide patterns (CUPP)', _Biotechnology for Biofuels_, 12(102).  https://doi.org/10.1186/s13068-019-1436-5
> Xu, J., Zhang, H., Zheng, J., Dovedo, P., Yin, Y. (2020) 'eCAMI: simultaneous classification and motif identification for enzyme annotation', _Bioinformatics_, 36(7), pp. 2068-2075 https://doi.org/10.1093/bioinformatics/btz908

This page covers includes the specific method used to perform the presented evaluation of the CAZyme classifiers listed above, including additional supplemenatry material to faciltiate the reporduction of the analysis.

Information for the general operation of `pyrewton` and support troubleshooting can be found at [Read the Docs](https://phd-project-scripts.readthedocs.io/en/latest/).

## Table of Contents

- [Program overview and installation](#program-overview-and-installation)
  - [Files and directories](#files-and-directories)
  - [Installation](#installation)
  - [Requirements](#requirements)
- [Independent evaluation of CAZyme classifiers](#independent-evaluation-of-cazyme-classifiers)
  - [Creation of the test sets](#creation-of-test-sets)
  - [Invoking the classifiers](#invoking-the-classifiers)
  - [Calculating the statistics](#calculating-the-statistics)
  - [Presenting the findings](#presenting-the-findings)
- [Creating a local CAZome database](#creating-a-local-cazome-database)
- [Notebooks](#notebooks)
- [Help](#help)

## Program overview and installation

### Files and directories
  
The `pyrewton` repository is separated into documentation and program modules.

- `assets`: Directory containing all files needed for the GitHub page.
- `docs`: Directory containing all files used for building the [documentation](https://pyrewton.readthedocs.io/en/latest/).
- `notebooks`: Directory containing all jupyter notebooks. The notebooks present the operation of some scripts, providing an explanation of the specific computational processes of some scripts.
- `tests`: Directory containing all unit-tests, and required test inputs and targets.
- `supplementary`: Directory containing supplemenatry data, accompanying the presentation of the results in _Hobbs et al., 2021_.
- `pyrewton`: Directory containing all `pyrewton` program modules (including all submodules and Python scripts).
  - `cazymes`: Module for predicting CAZy family annotations and evaluating CAZyme prediction tools
  - `genbank`: Module for retrieving and parsing data from [NCBI GenBank](https://www.ncbi.nlm.nih.gov/genbank/)
  - `utilities`: Module containing program utility functions and methods
    - `file_ui`: Submodule for handling file inputs and outputs
    - `parsers`: Submodule for building cmd-line args parsers
  
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
 
Following this method ensures all requirments are installed, as well as installing the CAZyme classifers into the correct directorys.
  
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

## Creating a local CAZome database
 
  
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

### Help

More indepth documentation is hosted at [Read the Docs](https://pyrewton.readthedocs.io/en/latest/). This includes listing all optional flags for configuring the operation of `pyrewton` and more details on the operation of each script/entry point in `pyrewton`.
 
Please raise any issues with any of the programmes at the GitHub issues pages for this repository, by following [the link](https://github.com/HobnobMancer/PhD_Project_Scripts/issues).
