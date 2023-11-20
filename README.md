# pyrewton

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3876218.svg)](https://doi.org/10.5281/zenodo.3876218)
[![Funding](https://img.shields.io/badge/Funding-EASTBio-blue)](http://www.eastscotbiodtp.ac.uk/)
[![PhD licence](https://img.shields.io/badge/Licence-MIT-green)](https://github.com/HobnobMancer/PhD_Project_Scripts/blob/master/LICENSE)
[![CircleCI](https://circleci.com/gh/HobnobMancer/pyrewton.svg?style=shield)](https://circleci.com/gh/HobnobMancer/PhD_Project_Scripts)
[![codecov](https://codecov.io/gh/HobnobMancer/pyrewton/branch/master/graph/badge.svg)](https://codecov.io/gh/HobnobMancer/pyrewton)
[![Documentation Status](https://readthedocs.org/projects/pyrewton/badge/?version=latest)](https://pyrewton.readthedocs.io/en/latest/?badge=latest)
[![Python](https://img.shields.io/badge/Python-v3.7.---orange)](https://www.python.org/about/)
[![Research](https://img.shields.io/badge/Bioinformatics-Protein%20Engineering-ff69b4)](http://www.eastscotbiodtp.ac.uk/eastbio-student-cohort-2019)

_Please find more detailed documentation for operation and troubleshooting at [Read the Docs](https://phd-project-scripts.readthedocs.io/en/latest/)_

## Contents

1. [Overview](#Overview)
2. [Installation](#Installation)
    - [Requirements](#Requirements)
3. [Development](#Development)
3. [Directories](#Directories)
4. [Modules](#Modules)
    - [genbank](#genbank)
    - [cazymes](#cazymes)
        - [uniprot](#uniprot)
        - [evaluate_tools](#evaluate_tools)
5. [Evaluations](#Evaluations)
        
## Overview

Pyrewton is a Python3 script package for the automated identification of CAZyomes (all carbohydrates encoded within the genome of a given species). The package is run at the command line and free to use under the MIT license.

Pyrewton supports:
- Downloading of all genomic assemblies (as GenBank files .gbff) from the [NCBI Assembly database](https://www.ncbi.nlm.nih.gov/assembly)
associated with each species passed to the programme
- Retrieval of all annotated protein sequences from GenBank (.gbff) files
- Retrieve proteins entries from [UniProtKB](https://www.uniprot.org/), using a JSON file to configure the queries

Features currently in development:
- Use the 3rd-party tools [dbCAN](https://github.com/linnabrown/run_dbcan), [CUPP](https://www.bioengineering.dtu.dk/english/researchny/research-sections/section-for-protein-chemistry-and-enzyme-technology/enzyme-technology/cupp), and [eCAMI](https://github.com/zhanglabNKU/eCAMI) to predict the which query protein sequences are CAZymes and predict their CAZy family
- Evaluate the accuracy of the CAZyme prediction tools to distinguish between CAZyme and non-CAZyme protein sequences
- Evaluate the accuracy of the CAZyme prediction tools to the correct CAZy family
- Produce a report of the CAZyme prediction tool evaluation

Development plans are stored within the [Wiki](https://github.com/HobnobMancer/pyrewton/wiki).
<p>&nbsp;</p>

If you're coming from the [Microbiology Society Conference poster](https://doi.org/10.6084/m9.figshare.14370836.v1), and want to see how to repeat the analysis, please navigate to `pyrewton/cazymes/precition`. This directory contains an additional README specifically for the independent evaluation of the CAZyme predictiont tools dbCAN, CUPP and eCAMI.

In the above study, Hobbs _et al._, 2021, the following version of the CAZyme classifiers were evaluated:
- dbCAN: v2.0.11
- CUPP: v1.0.14
- eCAMI: no version given, April 2020 release

## Installation

The easiest method is to use pip to install `pyrewton` and all requirements.

1. Create a virtual environment with dependencies, then activate the environment - _where venv_name is an chosen name for the virtual environment_
`conda create -n <venv_name> python=3.8 diamond hmmer prodigal -c conda-forge -c bioconda`   
`conda activate <venv_name>`

2. Clone the repository
`git clone https://github.com/HobnobMancer/pyrewton.git`

3. Install pyrewton
`pip3 install -e <path to directory containing setup.py file>`   
Do not forget to use the **-e** option when install using pip3, otherwise each time pyrewton is invoked a ModuleNotFound error will be raised. Pass the path to the **directory** containign the setup.py file not the path to the setup.py file; if you are currently in the root directory of the repoistory where the file is located, simply use '.' to indicate the current working directory.

4. Install the third party CAZyme prediction tools
The easiest way to do this, and ensure they are installed into the correct directories is to use:  
`python3 <path to pyrewton setup.py> cpt -p <path to directory to install tools>`

Each tool is installed to a subdirectory in the specified directory, such that if the specified directory was `tools/` the resulting structure would be:
```bash
tools/
    - dbcan/
    - ecami/
    - cupp/
```

For alternative methods of installation see the full documentation at [Read the Docs](https://phd-project-scripts.readthedocs.io/en/latest/).

<p>&nbsp;</p>

## Requirements

POISx or Mac OS, or linux emulator   
Python version 3.7+   
Miniconda3 or Anaconda managed microenvironment, incorporated code checkers are included in list form in 'requirements.txt'.   
Miniconda3 environment file is also available in the GitHub repository: 'environment.yml'.   
For all required Python libraries please read 'requirements.txt'.   

<p>&nbsp;</p>


## Development

This section of the README lists the areas that are currently being worked upon and expanded:
- update unit tests to cover the newly added scripts, and match the last repo restructuring
- continue indepth evaluation of the CAZyme prediction tools, visualising the data in an Rnotebook
- writing documentation to be hosted at ReadTheDocs for detailed instructions on invoking `pyrewton`

<p>&nbsp;</p>


# Documentation

Full documentation can be found at [Read the Docs](https://pyrewton.readthedocs.io/en/latest/?badge=latest).

Here is a summary of using the `pyrewton` pipeline to automate annotating and exploring CAZymes and CAZomes of species of interest.

## Download genomes.

Use the `download_genomes` subcommand to download the genomic assemblies of all genomes available for a set of candidate species. 

* Download GenBank assemblies to increase the probability of finding CAZy annotated CAZymes
* If a GenBank assembly is not available, the RefSeq assembly is retrieved
* Candidate species are specified using a plain text file
    * List a unique species per line
    * Identified species by their scientific name or NCBI Taxonomy ID
    * If a species is specified, assemblies for all its stains are downloaded
* Genomes downloaded in .GBFF format
* Leave genomes compressed if working with `pyrewton` - to help manage disk space
* A CSV file lising the scientific name, NCBI Taxonomy ID and the assembly version accession of all downloaded genomes is created to log the process

**Note:** All NCBI taxonomy IDs need the prefix 'NCBI:txid'. For example "NCBI:txid318829" not "318829"

**To include comments:** To include a comment (i.e. a piece of text not to be read by `pyrewton`), start the line with a hashtag '#'. For example:
```bash
# Fungal species
Aspergillus fumigatus
Aspergillus nidulans
Aspergillus niger
Aspergillus sydowii
# Oomycete species
Albugo candida
Hyaloperonospora arabidopsidis
Phytophthora cinnamomi
```

An example/template input file is included in the GitHub repo at `templates/get_ncbi_genomes_template_input_file`.

An example output CSV file listing the genomic assemblies is included in the GitHub repo at `templates/example_genome_file.csv`.

Example command:
```bash
pyrewton download_genomes \
    $1 \
    -d data/genomes/genome_dataframe.csv \
    -i data/species/species_list \
    -o data/genomes/ 
```

Flags: 
```bash
positional arguments:
  user email address    Email address of user, this must be provided for Entrez

options:
  -h, --help            show this help message and exit
  -d DATAFRAME, --dataframe DATAFRAME
                        Location of file for species table to be written to, include .csv extention (default: <_io.TextIOWrapper name='<stdout>'
                        mode='w' encoding='utf-8'>)
  -f, --force           Force file over writting (default: False)
  -g, --genbank         Disable pulldown of GenBank (.gbff) files (default: True)
  -i input file name, --input_file input file name
                        Input filename (default: <_io.TextIOWrapper name='<stdin>' mode='r' encoding='utf-8'>)
  -l log file name, --log log file name
                        Defines log file name and/or path (default: None)
  -n, --nodelete        enable/disable deletion of exisiting files (default: False)
  -o output directory name, --output output directory name
                        Path to output directory. STDOUT. Directory will be created if needed. (default: <_io.TextIOWrapper name='<stdout>' mode='w'     
                        encoding='utf-8'>)
  -r maximum number of retries, --retries maximum number of retries
                        Defines the maximum number of retries if network errors are encountered (default: 10)
  -t TIMEOUT, --timeout TIMEOUT
                        Timeout for URL connections, in seconds (default: 10)
  -v, --verbose         Set logger level to 'INFO' (default: False)
```

## Extract protein sequences

The `extract_protein_seqs` subcommand extracts protein sequences from all CDS features in local, compressed GenBank Flat File Format (.gbff) genomic assemblies (.gbff.gz).

* Extracts the protein sequence, protein version accession, locus tag, locus / location, and functional annotation
* Creates a multiple sequence FASTA file of extracted protein sequences per genomic assembly
* Each output FASTA file is named with the corresponding NCBI genomic version accesion
* Creates a CSV file summarising the extraction, i.e. the number of proteins retrieved from each genome

An example input CSV file listing the genomic assemblies is included in the GitHub repo at `templates/example_genome_file.csv`.

Example command:
```bash
pyrewton extract_protein_seqs \
    data/genomes/genome_dataframe.csv \
    data/genomes/ \
    data/proteins/proteomes
```

Flags:
```bash
positional arguments:
  input dataframe name  Path to input dataframe
  genome_directory      Path to directory containing compressed genomic assemblies in GenBank Flat File Format (.gbff.gz)

options:
  -h, --help            show this help message and exit
  -d OUTPUT_DF, --output_df OUTPUT_DF
                        path to output directory to write FASTA files to (default: <_io.TextIOWrapper name='<stdout>' mode='w' encoding='utf-8'>)        
  -f, --force           Force file over writting (default: False)
  -l LOG, --log LOG     Defines log file name and/or path (default: None)
  -n, --nodelete        enable/disable deletion of exisiting files (default: False)
  -o OUTPUT, --output OUTPUT
                        Output directory. Path to directory to which FASTA files are written (default: <_io.TextIOWrapper name='<stdout>' mode='w'       
                        encoding='utf-8'>)
  -v, --verbose         Set logger level to 'INFO' (default: False)
```

## Automate running CAZyme classifiers to predict CAZymes

The `pyrewton` subcommand `predict_cazymes` coordinates `pyrewton` to automate running the specified CAZymes classifies (from CUPP, dbCAN and eCAMI) for each multi-sequence FASTA file in an input directory.

* All tools (CAZyme classifiers) must be installed to an individual subdirectory for each tool, and which are all located in the same parent directory.
* Name as many tools to run as desired, in a space separated list.
* Provide the path to the parent directory containing all the tools
* `pyrewton` presumes the tool subdirectories are written in lower case. If not (e.g. CUPP was downloaded to a directory called `Cupp/`), call dbCAN as 'Cupp' when using the `predict_cazymes` subcommand.
* If no output directory is specified the output is written to the current working directory
* A new output directory is created for each multi-sequence FASTA file parsed by the CAZyme classifiers. Inside this subdirectory is the output from each classifier. Each subdirectory is named with the corresponding NCBI genomic version accesion.

Example command (presuming a directory called dbcan/ containing the dbcan databases is availabe in a directory called `tools/`):
```bash
pyrewton predict_cazymes \
    data/proteins \
    dbcan,
    tools \
    -o data/dbcan_output
```

flags:
```bash
positional arguments:
  input directory       Path to directory containing FASTA files for prediction tools
  {cupp,ecami,dbcan}    CAZyme classifiers to run. Pick as many as wanted. Case insensitive
  tool_dir              Path to parent directory where CAZyme classifiers are installed. E.g. point to the parent directory containing the directory called 'dbcan/'

options:
  -h, --help            show this help message and exit
  -f, --force           Force file over writting (default: False)
  -l log file name, --log log file name
                        Defines log file name and/or path (default: None)
  -n, --nodelete        enable/disable deletion of exisiting files (default: False)
  -o output directory, --output output directory
                        Directory to which all outputs are written (default: <_io.TextIOWrapper name='<stdout>' mode='w' encoding='utf-8'>)
  -v, --verbose         Set logger level to 'INFO' (default: False)
  ```

## Create a comprehensive CAZome database

Compile predicted CAZyme classifications from CAZyme classifiers and canconical CAZyme classifications from CAZy into a single, shareable, local SQLite3 database using the `compile_cazome_db` subcommand.

To specify which CAZyme classifers are used, a config YAML file is provided to `pyrewton`. This YAML file is keyed by the names of the tools, and under each tool name are tool additional keys:
* `version:` the version number of the tool
* `cazy_release:` the approximate data of the CAZy database release the tool 
was trained against. If unknown leave blank
* `dir:` path to the directory containing the output directories
For example:
```yaml
dbCAN_2:
  version: 2.0.11
  cazy_release: March 2020
  dir: data/dbcan_2
dbCAN_3:
  version: 3.0.7
  cazy_release: April 2022
  dir: data/dbcan_3
```

The names of the tools will appear as is (as stated in the YAML file) in the resulting local CAZome database. If using multiple version of dbCAN, for clarity in the dataset, label these with different names. If not pyrewton will do this mannually using the version number provided.

`pyrewton` defines each classifier-CAZy family pair as a new domain, thus a single protein can be associated with multiple domains in the CAZome database.

### Flags
```bash
positional arguments:
  config_file           Path to directory containing FASTA files of protein seqs extract from genomic assemblies

options:
  -h, --help            show this help message and exit
  --new_db NEW_DB       Path to build a new CAZome database. Path to create database file. CANNOT be used at same time as --db (default: None)
  --db DB               Path to existing CAZome database to add data to. CANNOT be used at same time as --db (default: None)
  --genome_csv GENOME_CSV
                        Path to CSV file listing taxonomies and genomic accessions. Created using pyrewton download_genomes subcommand (default: None)
  --protein_dir PROTEIN_DIR
                        Path to directory containing FASTA files of protein seqs extract from genomic assemblies (default: None)
  --cazy CAZY           Path to local CAZyme db of CAZy annotations of proteins created using cazy_webscraper Will add CAZy annotations to the output. (default: None)
  --cazy_date CAZY_DATE
                        Date CAZy data was pulled down (default: None)
  -f, --force           Force file over writting (default: False)
  -l log file name, --log log file name
                        Defines log file name and/or path (default: None)
  -n, --nodelete        enable/disable deletion of exisiting files (default: False)
  --sql_echo            Set SQLite echo to True, adds verbose SQL messaging (default: False)
  -v, --verbose         Set logger level to 'INFO' (default: False)
```

### Adding CAZyme annotations to an existing local CAZome database

The `pyrewton` subcommand `compile_cazome_db` can also be used to add new CAZyme classifications from CAZy and new dbCAN, CUPP and eCAMI analyses to an existing CAZome database created using `pyrewton`. 

To do this, use the `--db` flag to provide a path to an existing CAZome database. Conversly, use the `--new_db` flag to provide a path to create a new CAZome database.

`pyrewton` will parse all specified CAZy, dbCAN, CUPP and eCAMI data, adding new proteins, genomes, taxonomies, CAZy families, classifiers and CAZyme domains to the local CAZome database without adding any duplicates to any tables.

# Repo structure

## Directories

Below is a directory plan of this repository, followed by a brief overview of each directories role , to facilitate navigation through the repository.

### **assets**

Directory containing all files needed for the GitHub page, created for easy access to accompanying Jupyter notebooks.

### **docs**

Directory containing files to build documentation hosted at ReadTheDocs.

### **notebooks**

Directory containing all Jupyter notebooks, and html copies used for easier in-browser viewing via the GitHub pages.

### **tests**

Directory containing all `pytest` files for testing `pyrewton`, including subdirectories for test inputs and targets. Each module/submodule has its own specific test input and target subdirectory.

### **pyrewton**

Directory containing all `pyrewton` program modules (including all submodules and Python scripts).
<p>&nbsp;</p>

## Modules

_Please find more detailed documentation at for operation and troubleshooting at [Read the Docs](https://phd-project-scripts.readthedocs.io/en/latest/)_

This is an overview of the functionalities of each module within `pyrewton`, as well as basics of operation. For more detailed documentation on the operation of each module and indiviudal Python scripts please see the documentation at [Read the Docs](https://phd-project-scripts.readthedocs.io/en/latest/)

### **utilities**

Contains all functions that are called from other Python scripts for building command-line parsers and loggers. Includes the submodule **file_io**, which contains functions that are called from other Python scripts for handling directories and files in `pyrewton`, including retrieving program inputs and creating output directories.

### **genbank**

Directory containing all submodules that are involved in retrieving handling GenBank files. This includes retrieval of GenBank files from GenBank, and retrieval of protein sequences from GenBank files.

#### **get_ncbi_genomes**
This submodule is for the retrieval of genomic assemblies (as .gbff files) for each species listed in a plain text file (.txt). Each line of the plain text file must contain a single, unique species, for an example see 'get_ncbi_genomes_template_input_file.txt' within the directory. The species can be specified by taxonomy ID (using the 'NCBI:txid' prefix) or scientific name. `get_ncbi_genomes` will retrieve the scientific name or taxonomy depending on which is given, and will write out a dataframe containing the scientific name, NCBI taxonomy ID and all accession numbers of all genomic assemblies retrieved from NCBI.

**Note:** What is meant by all genomic assemblies is the latest version of all genomic assemblies, taking preference for GenBank files over reference assemblies. If not GenBank (identified by the 'GCA' prefix) assembly if available then the latest version of the reference assembly (identified by the 'GCF' prefix) will be retrieved.

When invoking `get_ncbi_genomes` a user email must be provided. This is a requirements of Entrez, the search and retrieval system of NCBI, which is accessed during the retrieval of taxonomy information and genomic assemblies. 

An example of the basic operation is:
`python3 get_ncbi_genomes <user_email> <-i path_to_input_.txt> <-o directory_to_store_assemblies> <-d species_dataframe_output_path`

All command options can be viewed by using `python3 get_ncbi_genomes -h` or `python3 get_ncbi_genomes --help`, and at [ReadtheDocs](https://phd-project-scripts.readthedocs.io/en/latest/genbank.html#get-ncbi-genomes).

#### **get_genbank_proteins**
This submodule is for the retrieval of proteins sequences from GenBank (.gbff) genomic assemblies. The protein sequences are identified as 'CDS' annotated features. The locus tag, gene start/end, gene ID, annotated function and protein sequence are retrieved and written out to a dataframe, with a unique protein on each line and including a columns for the host species scientific name, NCBI taxonomy and accession number of the host genomic assembly. The protein sequences are also written out to FASTA files, with a single FASTA file containing all the protein sequences from only one genomic assembly. Therefore, each genomic assembly input results in one FASTA file output.

When invoking `get_genbank_proteins` the path to the input dataframe (which is the output from `get_ncbi_genomes`) and the directory containing the genomic assemblies must be parsed, and in this order.

An example of basic operation is:
`python3 get_genbank_proteins <path_to_input_df> <path_to_assemly_dir>`

All command options can be viewed by using `python3 get_ncbi_genomes -h` or `python3 get_genbank_proteins --help`, and at [ReadtheDocs](https://phd-project-scripts.readthedocs.io/en/latest/annotations.html).

### **cazymes**
This module is involved in the identification and prediction of CAZymes.

#### **uniprot**
This submodule retrieved protein from UniProtKB. The query criteria are configured by a YAML file. The configuration file incorporates two keys: 'tax_ids' and 'queries'. Under 'tax_ids' list the NCBI taxonomy ID of the species the search is to be restricted to. Under 'queries' list the queries to be performed (for each taxonomy ID if given) using the UniProt query [syntax](https://www.uniprot.org/help/text-search) and query [fields](https://www.uniprot.org/help/query-fields). If only the taxonomy ID is given then all proteins for that taxonomy ID will be given. For an example configuration file see 'uniprot_config.yaml' in the `uniprot` directory.

The submodule writes out a dataframe containing
- NCBI taxonomy ID of the host species
- scientific name of the host species
- UniProt entry ID
- UniProt entry name
- Protein names
- EC numbe
- Protein length in amino acids
- Protein mass (Da)
- Domains
- Protein families in external database such as its [CAZy](www.cazy.org/) family
- Gene Ontology ID, molecular function and biological process annotations
- Protein sequence

The submodule also writes out all protein sequences to FASTA files, with each query to UniProt producing a single FASTA file containing all the resulting protein sequences retrieved from the query. The FASTA file names follow the format: 'uniprot_{uniprot_query}_{time_stamp}'. If only the taxonomy ID is given then the taxonomy ID will be written in the FASTA file name in the place of 'uniprot_query'.

When invoking the submodule `uniprot`, invoke the script `get_uniprot_proteins.py` and the path to the configuration must be provided.

An example of basic operation is:
`python3 get_uniprot_proteins <path_to_config_file.yaml>`

All command options can be viewed by using `python3 get_uniprot_proteins -h` or `python3 get_uniprot_proteins --help`. 'Read the Docs' documentation coming soon!

#### **evaluate_tools**

_Detailed documentation for hosting at ReadTheDocs is still in development._

This submodule is for the prediction if a query protein sequence is a CAZyme or non-CAZyme, and the prediction of the CAZy family if the protein is predicated to be a CAZyme.  
This directory contains an additional README specifically for the independent evaluation of the CAZyme predictiont tools dbCAN, CUPP and eCAMI.
<p>&nbsp;</p>

## Evaluations

For a summary of the inital evaluation of the performance of dbCAN, CUPP and eCAMI (March 2021), see the [poster](https://hobnobmancer.github.io/pyrewton/pyrewton_2021_03_31_dbcancuppecami_evaluation.pdf) we presented at the Microbiology Society Annual Conference 2021.

The Rnotebook (as a Rmarkdown file and HTML file), containing the full evaluation of the dbCAN, CUPP and eCAMI and which was used to generate the figures for Microbiology Society Conference poster are located in `pyrewton/cazymes/evaluate/tools/report`. We recommend using the Rnotebook as a template for any evaluations you conduct yourself.

## Repository renamed 2020-10-05
**Note:** This repository was renamed from 'PhD_project_scripts' to 'pyrewton'.
