# PhD_Project_Scripts

## Contents

1. [Overview](#Overview)
2. [Installation](#Installation)
    - [Requirements](#Requirements)
3. [Scripts](#Scripts)
    - [Extract_genomes_NCBI](#Extract_genomes_NCBI)
4. [Troubleshooting and common errors](#Troubleshooting-and-common-errors)
    - [Extract_genomes_NCBI](#Extract_genomes_NCBI)
5. [Testing](#Testing)
    - [Extract_genomes_NCBI](#Extract_genomes_NCBI)

## Overview

PhD_Project_Scripts contains all scripts associated with the EASTBio Proteng project.

Completed scripts are packaged into `Proteng`.

`Proteng` is a Python3 script package that supports the pull down of genomic assemblies from the [NCBI Assembly database](https://www.ncbi.nlm.nih.gov/assembly).

[In GitHub](https://github.com/HobnobMancer/PhD_Project_Scripts) each section of the project is given its own directory and subdirectories within.

## Installation

The easiest way to install `proteng` is to use `pip`:\
`pip3 install Proteng`.\
This will install all required Python packages dependencies.

## Requirements

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

## Troubleshooting and common errors

This section covers common errors expected to arise when invoking\
each script, and the probably causes of each error.

### Extract_genomes_NCBI

The majoirty of issues will arise due to errors in the input file.
Always ensure the input file does not contain any errors or blank
lines. Additionally, always ensure the correct path to the input
file is provided.

**IOError**\
This error will occur if there is a network issue when using Entrez
to call to NCBI. The script will automatically retry the call.

If a network error occurs during the retrieval of a scientific name
or taxonomy ID the program will terminate, to avoid errors occuring
during downstream processing that will cause the program to
automatically terminate.

If a network error occurs during the retrievel of the accession
numbers, the program will exit the process of the retrieval of
the accession numbers and start processing the retrieval of the
taxonomy ID in the dataframe. This is so a second attempt to
retrieve the accession numbers can be performed by reinvoking
the script, only for the species for which the retrieval
previously failed, thus reducing demand on the NCBI database.

**FileNotFoundError**\
This error will occur is the incorrect path is provided as the
input argument at the command line, or no input argument is
provided and STDIN contains no data. Ensure the path includes
the file name, with extension.

If this error occurs the program will terminate.

**IndexError during scientific name retrieval**\
This occurs when Entrez fails to retrieve a sceintific name for
the given taxonomy ID.

This is potentially caused by a typo in the taxonomy id provided
in the input file.

If this error occurs the program terminates to avoid errors
occuring during downstream processing that will cause the
program to automatically terminate.

**IndexError during taxonomy ID retrieval**\
This occurs when Entrez fails to retrieve a taxonomy ID for
the given scientific name.

This is potentially caused by a typo in the species name in the
input file, or a typo in a taxonomy ID 'NCBI:txid' prefix,
causing the program to misinturpret the ID as a species name
and use it to try and retrieve a scientific name.

If this error occurs the program terminates to avoid errors
occuring during downstream processing that will cause the
program to automatically terminate.

**IndexError during assembly ID retrieval**\
This occurs when Entrez fails to retrieve assembly IDs from
NCBI.

This may be becuase there are no directly linked assemblies
for the given taxonomy ID. Check the NCBI Taxonomy database
to ensure there are 'directly' linked assemblies and not
only 'subtree' assemblies.

If this error occurs the program with exit the retrieve of
the assembly IDs and not retrieve the NCBI accession numbers.
This allows for troubleshooting using on the specie(s)
for which it is required, to reduce demand on NCBI.

**RuntimeError during posting of assembly IDs to NCBI**\
This error occurs when Entrez fails to post the retrieved
assembly IDs, causing it to fail to retrieve the document
summary.

This is potentially caused by incorrect formating of the
the assembly IDs or the request is too large for Entrez/
NCBI. If this is the case, repeat the procedure in batches.

If this error occurs the program with exit the posting of
the assembly IDs and not retrieve the NCBI accession numbers.
This allows for troubleshooting using on the specie(s)
for which it is required, to reduce demand on NCBI.

**RuntimeError or IndexError during retrieval of accession numbers**\
This error occurs when Entrez fails to retrieve the document
summary from NCBI.

This is potentially caused by incorrect formating of the
the assembly IDs or the request is too large for Entrez/
NCBI. If this is the case, repeat the procedure in batches.

If this error occurs the program with exit the retrieval of
the NCBI accession numbers.
This allows for troubleshooting using on the specie(s)
for which it is required, to reduce demand on NCBI.

## Testing

This section describes the tests that are included with the
program to ensure each of the scripts is operating as
expected.

Inputs and targets for tests are stored in the respective `test_inputs`
and `test_targets` directories. Each script is provided its own subdirectory
containing its associated tests input and target files. These subdirectories
are named using the contracted version of the scripts name, such that
`Extract_genomes_NCBI.py` respective test subdirectories are called `test_ext_gnm_ncbi`.

### Extract_genomes_NCBI

**test_input_reading.py**\
This tests the ability of the `parse_inputfile()` function's
ability to open and read a given input file, defined using
the `pathlib` `Path` module.

**test_name_and_id_retrieval.py**\
This test tests the functions `get_genus_species_name()`, `get_g_s_name_retry`,
`get_tax_id()` and `get_t_ix_retry`, using provided inputs. The result of the
Entrez calls to the NCBI database are compared against
expected results.

**test_accession_number_retrieval.py**\
This test test the `get_accession_numbers()`, `get_a_id_retry`, `post_a_ids_retry` and `get_a_n_retry` functions, to ensure that Entrez can call to the NCBI database. Owing to the frequent updating of the NCBI Assembly database it is not possible to
compare the results from the NCBI call against expected results.
