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

This page is published as supplementary material to the main [`pyrewton` repo](https://github.com/HobnobMancer/pyrewton), and the independent evaluation of the CAZyme classifiers [dbCAN](zhange _et al._, 2018), [CUPP](), and [eCAMI]() presented in [Hobbs et al., 2021]().

> Hobbs, E. E. M., Gloster, T. M., Chapman, S., Pritchard, L. (2021): Microbiology Society Annual Conference 2021. figshare. Poster. https://doi.org/10.6084/m9.figshare.14370836.v

> Zhang, H., Yohe, T., Huang, L., Entwistle, S., Wu, P., Yang, Z., Busk, P.K., Xu, Y., Yin, Y. (2018) 'dbCAN2: a meta server for automated carbohydrate-active enzyme annotation', _Nucleic Acids Res._, 46(W1), pp. W95-W101. doi: 10.1093/nar/gky418

> Barrett, K., Lange, L. (2019) 'Peptide-based functional annotation of carbohydrate-active enzymes by conserved unique peptide patterns (CUPP)', _Biotechnology for Biofuels_, 12(102).  https://doi.org/10.1186/s13068-019-1436-5

> Xu, J., Zhang, H., Zheng, J., Dovedo, P., Yin, Y. (2020) 'eCAMI: simultaneous classification and motif identification for enzyme annotation', _Bioinformatics_, 36(7), pp. 2068-2075 https://doi.org/10.1093/bioinformatics/btz908

In Hobbs _et al._, 2021, the following version of the CAZyme classifiers were evaluated:
- dbCAN: v2.0.11
- CUPP: v1.0.14
- eCAMI: no version given, April 2020 release

This page covers the specific method used to perform the presented evaluation of the CAZyme classifiers dbCAN, CUPP and eCAMI, including additional supplemenatry material to faciltiate the reporduction of the analysis, as presented in [Hobbs et al., 2021]().

Information for the general operation of `pyrewton` and support troubleshooting can be found at [Read the Docs](https://phd-project-scripts.readthedocs.io/en/latest/).

The supplementary tables and figures for the evaluation presented in Hobbs _et al_., can be [found here](https://hobnobmancer.github.io/pyrewton/pyrewton/cazymes/evaluate_tools/report/cazyme_prediction_tool_evaluation-2021-03-30.html).

## Citation

When using `pyrewton` please cite:

> Hobbs, E. E. M., Gloster, T. M., Chapman, S., Pritchard, L. (2021): Microbiology Society Annual Conference 2021. figshare. Poster. https://doi.org/10.6084/m9.figshare.14370836.v

## Table of Contents

<!-- TOC -->
- [pyrewton](#pyrewton)
- [Citation](#citation)
- [Program overview and installation](#program-overview-and-installation)
  - [Files and directories](#files-and-directories)
  - [Installation](#installation)
  - [Requirements](#requirements)
- [Independent evaluation of CAZyme classifiers](#independent-evaluation-of-cazyme-classifiers)
  - [Creation of the test sets](#creation-of-the-test-sets)
  - [Invoking the classifiers](#invoking-the-classifiers)
  - [Calculating the statistics](#calculating-statistics)
  - [Presenting the findings](#presenting-the-findings)
- [Creating a local CAZome database](#creating-a-local-cazome-database)
- [Notebooks](#notebooks)
- [Help](#help)
<!-- /TOC -->

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

1. Create a new virtual environment and install some dependencies.
  _To install Conda please see the [Conda documentation](https://docs.conda.io/en/latest/)._
  ```bash
  conda create -n <venve name> python=3.8 diamond hmmer prodigal -c conda-forge -c bioconda
  conda activate <venve name>
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
4. Install the CAZyme classifiers into the correct directories. This **must** be run from the root of the repository to ensure the CAZyme classifiers are written to the correct directories.
  ```bash
  python3 setup.py cpt -p .
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

The supplementary tables and figures for the evaluation presented in Hobbs _et al_., can be [found here](https://hobnobmancer.github.io/pyrewton/pyrewton/cazymes/evaluate_tools/report/cazyme_prediction_tool_evaluation-2021-03-30.html).

### Retrieving ground truth annotations

CAZy is used as the reference annotation database during the evaluation. CAZy family annotations of proteins retrieved from CAZy can be provided to `pyrewton` as an SQLite3 database created using [`cazy_webscraper`](https://github.com/HobnobMancer/cazy_webscraper). The instructions for using `cazy_webscraper` can be found in the `cazy_webscraper` documentation.

Alternatively, a JSON file, keyed by GenBank protein accessions and valued by list of CAZy family annotations can be provided. This can also be generated using `cazy_webscraper`.

It is recommended to use an SQLite3 database becuase `cazy_webscraper` incorprates a log of every scrape of CAZy, UniProt and GenBank which is performed to add data to the local CAZyme database, thus facilitate reproduction of the database and the subsequent evaluation of the CAZyme prediction tools.

### Creation of the test sets
  
The Python scripts `create_test_sets*.py` are used to generate the test sets. They are invoked using the same command structure:
```bash
create_test_sets*.py \
  <user email> \
  <yaml file path> \
  <path to cazy json or db> \
  <path to output dir>
```
The script `create_test_sets_from_db.py` is used when retrieving CAZy family annotations from a local CAZyme database created using [`cazy_webscraper`](https://github.com/HobnobMancer/cazy_webscraper).
The script `create_test_sets_from_dict.py` is used when retrieving CAZy family annotations from a JSON file keyed by GenBank accessions and valued by list of CAZy family annotations.
The input aurgements for both scripts are the same:
- User email address: Required by NCBI Entrez for downloading genomes to produce the test sets from
- YAML file path: Path to a YAMl file keyed by NCBI taxonomy IDs and valued by list of GenBank assemblies, this defines the genomes to be downloaded and used to create the test sets
- Path to CAZy JSON or local CAZyme database db file
- Path to an output directory to write the output to - _this output directory and its parents do not need to already exist_

Specifically, in the YAML file, each selected GenBank genomic assembly selected to be used to create a test set 
is represented as a list. The first element is the **genomic accession number**, the second element 
is the **assembly name**. For example:
```yaml
txid498019:
    - [GCA_003013715.2, ASM301371v2]
    - [GCA_008275145.1, ASM827514v1]
txid573826:
    - [GCA_000026945.1, ASM2694v1]
```

`create_test_sets_*.py` create an output directory, at the location specified by the user. Inside this output directory, four directories and a plain text file are produced:
- `alignment_scores`: contains `.csv` files of the BLAST all-versus-all scores of the selected CAZymes query against all non-CAZymes
- `extracted_proteins_seqs`: FASTA files containing all extracted protein sequences from the downloaded genomic assemblies, one FASTA  file is created per genome and contains all extracted protein sequences from that one genome
- `genomes`: contains the downloaded genomic assemlbies in `gbff.gz` format
- `test_sets`: contains the tests sets (FASTA) files for be used as input for each CAZyme classifier. One test set is created per genome.
- `cazome_coverage_<time stamp>.txt`: Contains the following headers and data:
  - Genomic_accession: The accession of the genomic assembly
  - Total_proteins: Total number of proteins extracted from the genomic assembly
  - Total_CAZymes: Total number of CAZy annotated CAZymes extracted from the genomic assembly
  - Genome_CAZome_percentage: Percentage of the proteome contained within the CAZome
  - CAZome_coverage_percengate: Percentage of the identified CAZome included in the test set compiled from the genomic assembly
  - CAZyme_sample_size: Number of CAZymes included in the test set compiled from the genomic assembly

The default test set size is 100 CAZymes and 100 non-CAZymes, which was used in Hobbs _et al._, 2021. To change the sample size, call the `--sample_size` flag and define the number of **CAzymes** to be included in the test set, an equal number of non-CAZymes wills be added to the test set. For example, `--sample_size 200` will produce test sets of 200 CAZymes and 200 non-CAZymes.

A full list of optional arguments are listed at [Read the Docs]().

The yaml file containing the genomic accessions and taxonomy IDs of the genomes selected for the inclusion in the study is lcoated in [`supplementary/mar_2021_eval/selected_genomes.yaml`](), as well as the JSON file of CAZy family annotations [`supplementary/mar_2021_eval/cazy_dict_2021-03-25.json`]() - _although we strongly recommend using the local CAZyme database appraoch for future evaluations_.
  
In total 70 test sets were created for the evaluation, 39 from Bacteria genomes and 31 from Eukaryote genomes. A complete list of the proteins (identified by their GenBank accession number) used in the test sets for the evaluation presented in Hobbs _et al._,_ is located in [`supplementary/mar_2021_eval/test_set_composition.json`]()
  
### Invokving the classifiers

The Python script `predict_cazymes.py` is used to invoke each CAZyme classifier (dbCAN, CUPP and eCAMI) for each test set created in the last step.  

`predict_cazymes.py` **must** be run from the `pyrewton/cazymes/evaluation/` directory. This is because dbCAN, CUPP and eCAMI use hard coded paths to find their respective datasets. This requires `pyrewton` to navigate to the correct directories for the CAZyme classifers to function properly and this cannot be achieved if the script `predict_cazymes.py` is not invoked in the correct directory.

`predict_cazymes.py` takes 2 positional arguments:
1. The path to the directory containing the test sets created using `create_test_sets_*.py`
2. The path to the output directory (this directory and it's parents do not need to already exist. `pyrewton` will construct these)

The command for `predict_cazymes.py` takes the following structure:
```bash
python3 predict_cazymes.py \
  <path to dir containing test sets> \
  <path to output diir>
```

Additional optional flags are laid out in the [documentation]().

The output printed to the terminal when invoking each of the CAZyme classifiers is still printed to the terminal when invoking the classifiers via `predict_cazymes.py`. Additionally, the terminal output is logged and written to a `.log` file for each CAZyme classifier and stored in the respecitive test set output directory.

### Calculating statistics

To statistical evaluate the performance of the CAZyme classifiers, the Python scripts `calculate_stats_from_*.py` are used.

Similar to when creating the test sets, known CAZy family annotations (the ground truths) can be provided to `calculate_stats_from_*.py` as either:
- a JSON file keyed by GenBank accession, and valued by set of CAZy family annotations (this can be compiled using [`cazy_webscraper`](https://github.com/HobnobMancer/cazy_webscraper))
- a path to a local CAZyme database created using [`cazy_webscraper`](https://github.com/HobnobMancer/cazy_webscraper)

The Python script `calculate_stats_from_dict.py` is used when using a JSON file ofCAZy family annotations.   
`calculate_stats_from_db.py` is used when using a local CAZyme SQLite3 database.

Both Python scripts have the same command structure:
```bash
python3 calculate_stats_from_*.py \
  <path to dir containing predict_cazymes.py output> \
  <path to dir containing create_test_sets_*.py output, containing the dirs alignment_scores and test_sets> \
  <path to JSON or SQLite3 database file of CAZy family annotations from CAZy>
```

`python3 calculate_stats_from_*.py` produces several output files, which can be grouped to 3 levels of evaluating/benchmarking CAZyme classifier performance:

**Binary CAZyme/non-CAZyme differentiation:** 
1. The directory `binary_classifications` contains a `.csv` file per input test set, with contains a unique protien on each row, listing the Protein_accession, Genomic_accession and binary CAZyme (1) and non-CAZyme (0) classifictaion per CAZyme prediciton tool and the ground truth retrieved from CAZy.
2. `binary_classification_evaluation_<date>.csv`, written in long form with the columns: Statistic_parameters, Genomic_assembly, Prediciton_tool, Statistic_value. Statistical parameters calcualted are specificity, sensitivity, precition, F1-score and accuracy.
3. `binary_bootstrap_accuracy_evaluation_<date>.csv`, written in long form with the columns: Genomic_accession, Prediction_tool, Bootstrap_number, accuracy. This contains the accuracy score of each bootstrap sample.
4. `binary_false_positive_classifications_<date>.csv`, with the columns: Genomic_accession, Protein_accession, then one column per prediction tool (listing a 0 if the tool did not annotate the protien as a CAZyme and a 1 if it did), BLAST_score_ratio (from creating the test sets), CAZyme_subject (the CAZyme if had the greatest sequence identitiy with in the genome when creating the test set).
5. `binary_false_negative_classifications_<date>.csv`, with the columns: Genomic_accession, Protein_accession, then one column per prediction tool (listing a 0 if the tool did not annotate the protien as a CAZyme and a 1 if it did), BLAST_score_ratio (from creating the test sets), CAZyme_subject (the CAZyme if had the greatest sequence identitiy with in the genome when creating the test set).

**CAZy class annotation prediction:**. Evaluating per CAZy class indpendent of other CAZy classes, and the CAZy class multi-label classification.  
1. `class_predicted_classification_<date>.csv`, written in long form with the columns: Genomic_accession, Protein_accession, Prediction_tool, then one column per CAZy class (which are scored as 0 if not predicted to be in the class, and 1 if predicted to be in the CAZy class), Rand_index and Adjusted_rand_index
2. `class_classification_stats_per_test_set_<date>.csv`, written in long form with the columns: Genomic_accession, Prediction_tool, CAZy_class, Statistic_parameter, Statistic_value. Statistical parameters calcualted are specificity, sensitivity, precition, F1-score and accuracy.
3. `class_ground_truths_classifications_<date>.csv` containing the CAZy class annotations from CAZy for every protein across all input tests sets.
4. `class_stats_across_all_test_sets_<date>.csv` contains performance statistics calculated when pooling all test sets into a single test set. Contains the columns: Prediction_tool, CAZy_class, Specificity, Sensitivity, Precision, Fbeta_score, and accuracy.
5. `class_stats_per_test_set_<date>.csv` contains the calculated performance statistics when evaluating performance per test set. Written in long form, and contains the columns: Genomic_accession, Prediction_tool, CAZy_class, Statistic_parameter, Statistic_value.

**CAZy family annotation prediction:**. Evaluating per CAZy family indpendent of other CAZy families, and the CAZy family multi-label classification.  
1. `family_predicted_classification_<date>.csv`, written in long form with the columns: Genomic_accession, Protein_accession, Prediction_tool, then one column per CAZy family (which are scored as 0 if not predicted to be in the class, and 1 if predicted to be in the CAZy class), Rand_index and Adjusted_rand_index
2. `family_long_form_stats_df_<date>.csv`, written in long form with the columns: CAZy_family, Prediction_tool, Statistical_parameter, Statistic_value. The statistic patermeters included are:
  - Specificity
  - Sensitivity
  - Precision
  - Fbeta-score
  - Accuracy
  - Number of true negatives (TN) in the data set
  - Number of false negatives (FN) in the data set
  - Number of true positives (TP) in the data set
  - Number of false positivies (FP) in the data set
  - Sensitivity sample size (TP + FN)
4. `family_per_row_stats_df_<date>.csv`, with the columns: CAZy_family, Prediction_tool, Specificity, Sensitivity, Precision, Fbeta_score, Accuracy.
6. `family_ground_truth_classifications_<date>.csv` containing the CAZy family annotations from CAZy for every protein across all input tests sets.
7. `CAZy_fam_testset_freq_<date>.json` keyed by CAZy family and valued by interget listing the number of occurences of the respective CAZy family across all test sets.

By default `python3 calculate_stats_from_*.py` writes the output to the cwd. To specify a different output directory, add the `--output` flag to the command, followed by the path to the output directory. The output directory and its parent directories do not need to already exist, these can be created `pyrewton`.

#### Optional outputs

In addition to the default outputs, additional evaluations can be performed.

**Comparing performance per taxonomy group:**
The test sets may cover a range of taxonomy groups. `pyrewton` can be used to evaluate the perforamnce across all test sets as well as per user-defined taxonomy group.

Taxonomy groups can be defined using a `YAML` file, keyed by the name of the taxonomy group, and valued by a list of the genomic accessions of test sets. For example:
```yaml
"Bacteria":
    - GCA_000021645.1
    - GCA_000015865.1
    - GCA_000237085.1
    - GCA_016406125.1
    - GCA_013283915.1
"Eukaryote":
    - GCA_014784225.1
    - GCA_009017415.1
    - GCA_016861735.1
    - GCA_013426205.1
    - GCA_001592805.2
    - GCA_016767815.1
```
The taxonmy group file is provided to `pyrewton` using `--tax_groups` flag. For example:
```bash
python3 calculate_stats_from_*.py \
  predicted_cazymes \
  test_sets_dir \
  cazy.db \
  --tax_groups evaluation_tax_groups.yaml
```
The yaml file of CAZyme classifier combinations presented in Hobbs _et al_., 2021 is stored in `supplementary/test_set_tax_groups.yaml`.

Using the `--tax_groups` flag results in `pyrewton` producing additional output:

1. `binary_classification_tax_comparison_<date>.csv`, which is the same as the `binary_classification_evaluation_<date>.csv` dataframe with an additional 'Tax_group' column, listing the taxonomy group as defined in the `YAML` file for each genomic assembly.
2. `class_classification_tax_comparison_<date>.csv`, which is the same as the `class_predicted_classification_<date>.csv` dataframe with an additional 'Tax_group' column, listing the taxonomy group as defined in the `YAML` file for each genomic assembly.
3. `class_stats_all_test_sets_tax_comparison_<tax_group>_<time_stamp>.csv`, which is the same as the `class_stats_across_all_test_sets_<data>.csv` dataframe with an additional 'Tax_group' column and statistics calculated for only test sets belong to the given taxonomy group, as defined in the `YAML` file. One dataframe is generated for each taxonomy group in the `YAML` file.
4. `class_stats_per_test_set_tax_comparison_<tax_group>_<time_stamp>.csv`, which is the same as the `class_stats_per_test_set_<data>.csv` dataframe with an additional 'Tax_group' column and statistics calculated for only test sets belong to the given taxonomy group, as defined in the `YAML` file. One dataframe is generated for each taxonomy group in the `YAML` file.
5. `family_classification_tax_comparison_<date>.csv`, which is the same as the `family_predicted_classification_<date>.csv` dataframe but with an additional 'Tax_group' column, listing the taxonomy group as defined in the `YAML` file for each genomic assembly.
6. `family_per_row_stats_tax_comparison_<tax_group>_<date>.csv`, which is the same as the `family_per_row_stats_df_<date>.csv` dataframe but only contains the statistical parameters for test sets belonging to given taxonomy group, as defined in the `YAML` file. One dataframe is generated for each taxonomy group in the `YAML` file.
7. `family_long_form_stats_tax_comparison_<tax_group>_<date>.csv`, which is the same as the `family_long_form_stats_df_<date>.csv` dataframe but only contains the statistical parameters for test sets belonging to given taxonomy group, as defined in the `YAML` file. One dataframe is generated for each taxonomy group in the `YAML` file.

**Recombine classifiers:**
`pyrewton` can evaluate the performance of the concensus result of any three CAZyme classifiers evaluated by `pyrewton` (where the consensus result is defined as the result at least two of the three tools agree upon). For example, Hotpep in dbCAN, cam be substituted for eCAMI and/or CUPP, to determine if using the newer k-mer CAZyme classifiers improves the performance of dbCAN.

To define the three classifiers, use the ``--recombine`` flag and pass a text file with a unique combination of CAZy classifiers on each line. For example, the input file may contain:  
```bash
DIAMOND,HMMER,CUPP
DIAMOND,HMMER,eCAMI
```
to evaluate the impact of substituting Hotpep for newer k-mer-based CAZyme classifiers on dbCAN perforamnce.

Below is an example command for using `pyrewton` to recombine the CAZyme classifiers:
```bash
python3 calculate_stats_from_*.py \
  predicted_cazymes \
  test_sets_dir \
  cazy.db \
  --recombine cazyme_classifiers_combiations.txt
```

The text file of CAZyme classifier combinations presented in Hobbs _et al_., 2021 is stored in `supplementary/cazyme_classification_recombinations.txt`.

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
