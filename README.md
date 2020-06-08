# EastBIO PhD Project Scripts, packaged into `pyrewton`

[![DOI](https://zenodo.org/badge/243783792.svg)](https://zenodo.org/badge/latestdoi/243783792)
[![Funding](https://img.shields.io/badge/Funding-EASTBio-blue)](http://www.eastscotbiodtp.ac.uk/)
[![PhD licence](https://img.shields.io/badge/Licence-MIT-green)](https://opensource.org/licenses/MIT)
[![CircleCI](https://img.shields.io/badge/CircleCI-Passing-brightgreen)](https://circleci.com/product/)
[![codecov](https://codecov.io/gh/HobnobMancer/PhD_Project_Scripts/branch/master/graph/badge.svg)](https://codecov.io/gh/HobnobMancer/PhD_Project_Scripts)
[![Python](https://img.shields.io/badge/Python-v3.7.---orange)](https://www.python.org/about/)
[![Research](https://img.shields.io/badge/Bioinformatics-Protein%20Engineering-ff69b4)](http://www.eastscotbiodtp.ac.uk/eastbio-student-cohort-2019)

## Contents

1. [Overview](#Overview)
2. [Installation](#Installation)
    - [Requirements](#Requirements)
3. [Directories, modules and files](#Directories,-modules-and-files)
    - [pyrewton modules](#pyrewton-modules)
4. [Troubleshooting and common errors](#Troubleshooting-and-common-errors)
    - [get_ncbi_genomes](#get_ncbi_genomes)
5. [Testing](#Testing)
    - [Extract_genomes_NCBI](#Extract_genomes_NCBI)

## Overview

The repository PhD_Project_Scripts contains all scripts associated with the EASTBio PhD project 'Identifying Engineering Candidates for Advanced Biocatalysis in Biofuel Production'.

Completed scripts are packaged into `pyrewton`, which is a programme free to use under the MIT license. `pyrewton` is a Python3 script package that can be run at the command line. `pyrewton` is built up of multiple Python modules, which either perform a 'housekeeping' task such as logger building, or perform a main operation of `pyrewton`, such as coordinating the downloading of GenBank files. All modules, submodules and associated Python scripts are located within the `pyrewton` directory. Specifically, it supports the downloading of all genomic assemblies (as GenBank files .gbff) from the [NCBI Assembly database](https://www.ncbi.nlm.nih.gov/assembly) associated with each species passed to the program, and summarising the cazyme annotations within those GenBank files.

## Installation

The easiest way to install `pyrewton` is to use `pip`:  
`pip3 install -e <path to directory containing pyrewton setup.py>`.  
This will install all required Python packages dependencies.

## Requirements

Scripts are written using:  
Python = 3.7+  
Miniconda3 managed microenvironment, incorporated code checkers are included in list form in 'requirements.txt'  
Miniconda3 environment file is also available in the GitHub repository: 'environment.yml'.

## Directories, modules and files

Below is a directory plan of the pyrewton module structure, followed by a brief overview of each directories role in the repository, to facilitate navigation through the repository.

```bash
.
├── codecov.yml
├── Computational_Overall_Flow_Chart.pptx
├── config.yml
├── conftest.py
├── coverage.xml
├── environment.yml
├── index.md
├── LICENSE
├── README.md
├── requirements.txt
├── setup.py
├── assets
│   ├── _documents required for GitHub pages_
├── tests
│   ├── test_accession_number_retrieval.py
│   ├── test_housekeeping.py
│   ├── test_input_reading.py
│   ├── test_inputs
│   │   └── gt_ncbi_gnms_test_inputs
│   │       ├── gt_ncbi_gnms_reading_test_input.txt
│   │       └── gt_ncbi_gnms_test_inputs.json
│   ├── test_name_and_id_retrieval.py
│   └── test_targets
│       └── gt_ncbi_gnms_test_targets
│           ├── gt_ncbi_gnms_test_df.csv
│           └── gt_ncbi_gnms_test_targets.json
└── pyrewton
    ├── __init__.py
    ├── directory_handling
    │   ├── __init__.py
    │   ├── output_dir_handling_main.py
    │   └── __pycache__
    │       ├── __init__.cpython-37.pyc
    │       └── output_dir_handling_main.cpython-37.pyc
    ├── loggers
    │   ├── __init__.py
    │   ├── logger_pyrewton_main.py
    │   └── __pycache__
    │       ├── __init__.cpython-37.pyc
    │       └── logger_pyrewton_main.cpython-37.pyc
    ├── parsers
    │   ├── __init__.py
    │   ├── parser_get_ncbi_genomes.py
    │   └── __pycache__
    │       ├── __init__.cpython-37.pyc
    │       └── parser_get_ncbi_genomes.cpython-37.pyc
    └── genbank
        └── get_ncbi_genomes
            ├── get_ncbi_genomes.py
            ├── get_ncbi_genomes_template_input_file.txt
            ├── __init__.py
            ├── __pycache__
            │   ├── get_ncbi_genomes.cpython-37.pyc
            │   └── __init__.cpython-37.pyc
            ├── selected_species_list.txt
            └── working_species_list.txt
```

### assets

Directory containing all files needed for the GitHub page.

### notebooks

Directory containing all jupyter notebooks, and html copies used for easier in-browser viewing via the GitHub pages.

### tests

Directory containing all `pytest` files for testing `pyrewton` during development, including subdirectories for test inputs and targets, with each module/submodule possessing its own specific test input and target subdirectory.

### pyrewton

Directory containing all `pyrewton` program modules (including all submodules and Python scripts).

#### pyrewton module: parsers

Directory containing all Python scripts for building command-line parsers.

#### pyrewton module: loggers

Directory containing Python scripts for building loggers.

#### pyrewton module: directory_handling

Directory containing all Python scripts for handling directories in `pyrewton` Python scripts, including retrieving program inputs and creating output directories.

#### pyrewton module: genbank

Directory containing all submodules that are involved in retrieving NCBI GenBank files and data from said files:

##### Genbank submodule: get_ncbi_genomes

Directory containing the submodule that takes a list of species and downloads all directly linked GenBank (.gbff) files in the NCBI Assembly database.

### pyrewton modules

Below describes the operation and functionality of each pyrewton module/submodule, including a listing of the command-line options.

#### get_ncbi_genomes

`get_ncbi_genomes` takes a plain text file containing species scientific names or NCBI taxonomy as input. The script will find the corresponding taxonomy ID or scientific name, as appropriate, as well as retrieve all directly linked accession numbers, generating a dataframe of 'Genus', 'Species', 'NCBI Taxonomy ID' and 'NCBI Accession number'.

**Invoking the script from the command line:**  
**Compulsory arguments:**  
`-u, --user`  
&emsp;&emsp;Although indicated as optional, Entez requires an email address must be provided. If not provided the program will log this as an error and terminate.

**Optional arguments:**  
`-d, --dataframe`  
&emsp;&emsp;Specify output path for dataframe, which will be saved as a .csv file (inclusion of file extensions is optional). If not provided dataframe will be written out to STDOUT.

`-f, --force`  
&emsp;&emsp;Enable writing in specified output directory if output directory already exists.

`-g, --genbank`  
&emsp;&emsp;Enable or disable downloading of GenBank files.

`-h, --help`  
&emsp;&emsp;Display help messages and exit

`-i, --input`  
&emsp;&emsp;Specify input filename (with extension) input file.
If not given input will be taken from STDIN.

`-l, --log`  
&emsp;&emsp;Specify name of log file (With extension).
If only filename is given, log file will be written out
to the current working directory, otherwise provide path
including filename.
If not option is given no log file will be written out,
however, logs will still be printed to the terminal.

`-n, --nodelete`  
&emsp;&emsp;Enable not deleting files in existing output directory. If not enabled, output directory exists and writing in output directory is 'forced' then files in output directory will not be deleted, and new files will be written to the output directory.

`-o, --output`  
&emsp;&emsp;Specify filename (with extension) of output file.
If not option is given output will be written to STDOUT.

`-r, --retries`  
&emsp;&emsp;Specifiy maximum number of retries before cancelling call to NCBI
if a network error is encountered. The default is a maximum of 10 retries. When
maximum is reached, a value of 'NA' is returned.

`-t, --timeout`  
&emsp;&emsp;Specify timeout limit of URL connection when downloading GenBank files. Default is 10 seconds.

`-v, --verbose`  
&emsp;&emsp;Enable verbose logging - changes logger level from WARNING to INFO.

## Troubleshooting and common errors

This section covers common errors expected to arise when invoking\
each script, and the probably causes of each error.

### get_ncbi_genomes

The majority of issues will arise due to errors in the input file.
Always ensure the input file does not contain any errors or blank
lines. Additionally, always ensure the correct path to the input
file is provided.

**IOError**  
This error will occur if there is a network issue when using Entrez
to call to NCBI. The script will automatically retry the call the set
maximum number of times.

If the maximum number of retries is met before connecting to NCBI
without encountering a network error, 'NA' is returned and stored
in the dataframe.

**FileNotFoundError**  
This error will occur is the incorrect path is provided as the
input argument at the command line, or no input argument is
provided and STDIN contains no data. Ensure the path includes
the file name, with extension.

If this error occurs the program will terminate.

**IndexError during scientific name retrieval**  
This occurs when Entrez fails to retrieve a scientific name for the given taxonomy ID.

This is potentially caused by a typo in the taxonomy id provided
in the input file.

If this error occurs the string 'NA' will be returned.

**IndexError during taxonomy ID retrieval**  
This occurs when Entrez fails to retrieve a taxonomy ID for
the given scientific name. Returns 'NA'.

This is potentially caused by a typo in the species name in the
input file, or a typo in a taxonomy ID 'NCBI:txid' prefix,
causing the program to misinterpret the ID as a species name
and use it to try and retrieve a scientific name.

If this error occurs the string 'NA' will be returned.
If no taxonomy ID is available for the retrieval of accession numbers,
the retrieval of accession numbers is cancelled, and a value of 'NA' is
returned.

**IndexError during assembly ID retrieval**  
This occurs when Entrez fails to retrieve assembly IDs from
NCBI.

This may be because there are no directly linked assemblies
for the given taxonomy ID. Check the NCBI Taxonomy database
to ensure there are 'directly' linked assemblies and not
only 'subtree' assemblies.

If this error occurs the program with exit the retrieve of
the assembly IDs and not retrieve the NCBI accession numbers,
and return the string 'NA'.
This allows for troubleshooting using on the specie(s)
for which it is required, to reduce demand on NCBI.

**RuntimeError during posting of assembly IDs to NCBI**  
This error occurs when Entrez fails to post the retrieved
assembly IDs, causing it to fail to retrieve the document
summary.

This is potentially caused by incorrect formatting of the
the assembly IDs or the request is too large for Entrez/
NCBI. If this is the case, repeat the procedure in batches.

If this error occurs the program with exit the posting of
the assembly IDs and not retrieve the NCBI accession numbers,
and return the string 'NA'.
This allows for troubleshooting using on the specie(s)
for which it is required, to reduce demand on NCBI.

**RuntimeError or IndexError during retrieval of accession numbers**  
This error occurs when Entrez fails to retrieve the document
summary from NCBI.

This is potentially caused by incorrect formatting of the
the assembly IDs or the request is too large for Entrez/
NCBI. If this is the case, repeat the procedure in batches.

If this error occurs the program with exit the retrieval of
the NCBI accession numbers, and return the string 'NA'.
This allows for troubleshooting using on the specie(s)
for which it is required, to reduce demand on NCBI.

## Testing

This section describes the tests that are included with the
program to ensure each of the scripts is operating as
expected.

Inputs and targets for tests are stored in the respective `test_inputs`
and `test_targets` directories. Each script is provided its own subdirectory
containing its associated tests input and target files. These subdirectories
are named using the contracted version of the scripts name, such that submodule
`get_ncbi_genomes` respective test subdirectories are called `gt_ncbi_gnms_test_targets`.

### get_ncbi_genomes

**test_housekeeping.py**  
This tests the building of the parser, logger and output directory via their respective functions, which must all operate correctly in order for the programme to be functional.

**test_input_reading.py**  
This tests the ability of the `parse_inputfile()` function's
ability to open and read a given input file, defined using
the `pathlib` `Path` module.

**test_name_and_id_retrieval.py**  
This tests the functions `get_genus_species_name()` and `get_tax_id()`, using provided inputs. The result of the Entrez calls to the NCBI database are compared against
expected results.

**test_accession_number_retrieval.py**  
This test test the function `get_accession_numbers()` to ensure that Entrez can call to the NCBI database. Owing to the frequent updating of the NCBI Assembly database it is not possible to compare the results from the NCBI call against expected results.
