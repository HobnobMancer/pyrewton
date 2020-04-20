# PhD_Project_Scripts

## Overview

PhD_Project_Scripts contains all scripts associated with the EASTBio Proteng project.

Completed scripts are packaged into `proteng`.

`proteng` is a Python3 script package that supports the pull down of genomic assemblies from the [NCBI Assembly database](https://www.ncbi.nlm.nih.gov/assembly).

[In GitHub](https://github.com/HobnobMancer/PhD_Project_Scripts) each section of the project is given its own directory and subdirectories within.

## Installation

The easiest way to install `proteng` is to use `pip`:\
`pip3 install proteng`.
This will install all required Python packages dependencies.

### Requirements

Scripts are written using:\
Python = 3.7\
Miniconda3 managed microenvironment, incorporated code checkers are included in list form in 'requirements.txt'\
Miniconda3 environment file is also available in the GitHub repository: 'environment.yml'.

## Scripts

A brief overview of each script.

### Extract_genomes_NCBI

Extract_genomes_NCBI takes a plain text file containing species scientific names or NCBI taxonomy as input. The script will find the corresponding taxonomy ID or scientific name, as appropriate, as well as retrieve all directly linked accession numbers, generating a dataframe of 'Genus', 'Species', 'NCBI Taxonomy ID' and 'NCBI Accession number'.

**Invoking the script from the command line:**\
**Positional arguments:**\
User email address: Email address of user must be provided

**Optional arguments:**\
`-h, --help`\
&emsp;&emsp;Display help messages and exit

`-i, --input`\
&emsp;&emsp;Specify input filename (with extension) input file.
If only the filename is supplied Extract_genomes_NCBI.py
will only look in the curent working directory, otherwise
provide path to input file.
If no option is given the default input file
"Extract_genomes_NCBI_input_file.txt", located in the
same directory as Extract_genomes_NCBI.py, will be used.

`-l, --log`\
&emsp;&emsp;Specify name of log file (With extension).
If only filename is given, log file will be written out
to the current working directory, otherwise provide path
including filename.
If not option is given no log file will be written out,
however, logs will still be printed to the terminal.

-`o, --output`\
&emsp;&emsp;Specify filename (with extension) of output file.
If only the filename is given, Extract_genomes_NCBI.py
the output file will be written to the current working
directory, otherwise provide path to desired directory,
including filename.
If not option is given an output file will be written to
the current directory with the standard name:
"Extract_genomes_NCBI_Date_Time", where _Date_ and _Time_
are the date and time script was invoked, respectively.
