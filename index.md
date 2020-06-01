# EastBIO PhD Project Supplementary Material
## Identifying, Characterising & Engineering Plant Cell Wall Degrading Enzymes for Enhanced Biocatalysts in Biofuel Production


[![Funding](https://img.shields.io/badge/Funding-EASTBio-blue)](http://www.eastscotbiodtp.ac.uk/)
[![PhD licence](https://img.shields.io/badge/Licence-MIT-green)](https://opensource.org/licenses/MIT)
[![CircleCI](https://img.shields.io/badge/CircleCI-Passing-brightgreen)](https://circleci.com/product/)
[![Python](https://img.shields.io/badge/Python-3.7.---orange)](https://www.python.org/about/)
[![Research](https://img.shields.io/badge/Bioinformatics-Protein%20Engineering-ff69b4)](http://www.eastscotbiodtp.ac.uk/eastbio-student-cohort-2019)

This repository contains files used during the EastBIO PhD Project, for the identification of engineering candidates to enhance biocatalysis for biofuel production.

The `Jupyter` notebook environments has been used to facilitate explanation the function and operation of the command-line python programmes used. These notebooks include explantory text and snippets of code extract to illustrate code function, and thus not all code snippets are runnable in the `Jupyter` environment.

## Table of Contents

1. [Files and Directories](#Files_and_Directories)
2. [Quick Start](#Quick_Start)
3. [Notebooks](#Notebooks)
4. [Help](#Help)

### Files and Directories

- `Section1_Extracting_Genomes/`: directory containing Python script and input file for downloading GenBank files for species of interest.
- `Section2_CAZy_genbank_summary/`: directory containing Python script for summarising the already annotated CAZy content within GenBank files.

### Quick Start

The quickest way to install the programme `ProtEng` is to use pip: `pip3 install Proteng -e`. This will install all required Python packages and dependicies.\
**Note:** [Conda](https://docs.conda.io/en/latest/) will need to be installed on your system.

### Notebooks

Jupyter notebooks are accessible via these pages, but also via the terminal if the programme `Proteng` is installed on your system.

From the top level of the repository directory on your system, start the `Jupyter` notebook server by issuing the command:
`jupyter notebook`
This will opena new browser window or tab containing the Jupyter homepage, with a listing of all files and directories as laid out in the repository.\
Navigate to the 'Notebook/' directory, containing all Jupyter notebooks for the repository. To open a notebook simple click on the title of the notebook.

For more information please see the Jupyter notebook [Quick-start guide](https://jupyter-notebook-beginner-guide.readthedocs.io/en/latest/?fbclid=IwAR1yIwkYCDjcw5FJZ7CfKES3l72HubqGYGcFrVrUKwWZoYh4NHy3VVu0AgQ) and [tutorial](https://www.tutorialspoint.com/jupyter/jupyter_quick_guide.htm).

Each notebook is named to match the number of the section of the programme to which it is associated, such that the first notebook (01_Extracting_GenBank_files) refers to Section1_Extracting_genomes within the repository, as made explictiy clear in the notebook itself.
Again note that these notebooks contain exerts of code to facilitate the explanation of the code architecture and function, as well as details provided on how the code was implemented during the project. Therefore, the code in many of the code cells is not runnable and replication or exploration of the data should be performed using the original scripts provided within the repository.

If you are new or have little experience of using command-line programmes and/or the Python coding language, please read the approrpirate notesbooks before using the associated scripts and posting queries/issues on the GitHub pages. The notebooks have been written so that those of all experience levels should be able to understand how and why the script operates.

### Help

Many of the common errors expected to arise during the operation of the scripts provided in this repository are named in the repository README, including the probable causes of these issues.

Please raise any issues with any of the programmes at the GitHub issues pages for this repositroy, by following [the link](https://github.com/HobnobMancer/PhD_Project_Scripts/issues).
