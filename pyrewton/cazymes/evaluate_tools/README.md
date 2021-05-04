
# README: Independent Evaluation of CAZyme Classifiers

This module is for the independent evaluation of the CAZyme classifiers dbCAN, CUPP and eCAMI.

All Python and Rnotebooks for creating the test sets, invoking the classifiers, evaluating the performance and presenting the performance are included.

The Python scripts and Rnotebooks were used for creating the data for our poster presented at the [Microbiology Society Annual Conference 2021](https://doi.org/10.6084/m9.figshare.14370860.v3).

Read on for instructions on how to use the provided scripts to independetly evaluate the CAZyme classifiers.

If you have any difficulties please do raise an issue in this repository!

## Steps for evaluating CAZyme classifiers
1. Creating test sets
2. Invoking classifiers
3. Evaluating performance
4. Presenting the evaluation

## [1] Creating the test sets

In order to create the test sets [`cazy_webscraper`](https://doi.org/10.6084/m9.figshare.14370860.v3) was invoked to create a JSON file, keyed by GenBank protein accession and valued by a list of CAZy family annotations. This is necessary to determine which proteins are CAZymes and which proteins are non-CAZymes becuase each created test set includes an equal number of CAZymes and non-CAZymes.

The genomic assebmlies to be used to create the test sets must be stored in a YAML file. The YAML file must be keyed by NCBI taxonomy ID, and valued by a list of genomic assemblies. Each assembly must be represented as a list (in square brackets) of the genomic accession and the assemly name of a genomic assembly.

The YAML file `evaluation_genomic_assemblies_2021_03.yaml` was used in the evaluation for the Microbiology poster, and can be used as a template for future evaluations.

To create the test sets use the following command:
```
Python3 create_test_sets_from_dict.py <user_email> <path_to_yaml_file> <path_to_cazy_json_file> <path_to_output_dir>
```

The user email address is a requirement of NCBI Entrez, which is used for automating the retrieval of genomic assemblies.

The output directory does not already have to exist, as long as the parent directory exists the output directory will be created by `create_test_sets_from_dict.py`. If the output directory does exist you will need to use the flag `--force` to write to it. If you also do not want to delete the content already in the output directory (which is the default behaviour), use the `--nodelete` flag.

To change the number of CAZymes included in the test set use the `--sample_size` flag followed by the number of CAZymes to include. An equal number of non-CAZymes will be retrieved as well.

## [2] Invoking classifiers

To invoke the predictions tools use the following command:
```
python3 predict_cazymes.py <path_to_dir_containing_test_sets> --output <output_dir>
```

The same rules for the output directory used for `create_test_sets_from_dict.py` apply here.

This command need only be invoked once because it will invoke the prediction tools for each of the test sets included in the `<path_to_the_dir_containing_test_sets>`.

## [3] Evaluating performance

To perform the statistical evaluation of the prediction tools use the command:
```
python3 calculate_stats.py <path_to_output_dir_from_step_2> <path_to_cazy_json_file> --output <path_to_output_dir>
```

This script will calculate the statistical parameters of the performance, and generates a series of `CSV` files, which can be used in step 4 for visualise the data.

The Fbeta-score is calculated during the evaluation, be default beta is 1, but it can be changed using the command: `--beta <number>`.

For changing the parameters of the bootstrap resampling of the CAZyme/non-CAZyme predictions to observe typically expected ranges in performance of each prediction tool, the number of test sets included in the bootstrapping can be changed by using `--bs_sample_size` (default is 6). To change the number of rounds of bootstrap resampling per prediciton tool, per test set, can be changed by using the flag `--bs_resampling` (the default is 100). The default options were used for to create the data in out Microbiology Society conference poster.

The same rules for the output directory used for `create_test_sets_from_dict.py` apply here.

Many of the created `CSV` files stored the data in the long format, which is preferred format for parsing data within R, which is used in step 4. For more info on long and tidy data please refer to [Tidy Data](https://cran.r-project.org/web/packages/tidyr/vignettes/tidy-data.html).


### Dataframes created for the binary evaluation

'binary_classification_evaluation_{time_stamp}.csv'
Dataframe written in the long format, including the statistical evaluation of each of the prediction tools. It contains the columns: Statistic paramter, genomic assembly, prediction tool, statistic value

`bootstrap_accuracy_evaluation_{time_stamp}.csv'
Contains the results of bootstrap resampling the binary predictions, to evaluate the expected range in performance of each prediction tool.

`binary_classifications_{time_stamp}_{prediction.source}.csv`
One dataframe of this format is made for each test set, and includes the columns Protein accession, dbCAN, HMMER, Hotpep, CUPP, eCAMI and CAZy.  
Each row contains a unique protein. The columns named after a prediction tool include the CAZyme (1) / non-CAZyme (0) predicted classification, and the 'CAZy' column includes the CAZyme (1) and non-CAZyme (0) classification of the protein by CAZy.


### Dataframes created for the CAZy class and CAZy family multilabel classification

`{tool}_fam_prediction_df_{time_stamp}.csv`
These dataframes containg the CAZy family classifications for each protein.

`cazy_fam_fbeta_scores_{time_stamp}.csv` includes the columns: CAZy family, Prediction tool, Genomic accession and Fbeta score. The data is presented in the long (tidy) format.

`fam_stats_df_{time_stamp}.csv` contains the column: CAZy family, stat_paramter, prediction tool and stat_value, and is written in the long format.

`fam_stats_df_single_{time_stamp}.csv` contains the statistical evaluation for each cazy family on a single line per CAZy family and prediction tool. The columns are: CAZy family, prediction tool, specificity, recall, Fbeta-score, and accuracy.

`class_ground_truths_classifications_{time_stamp}.csv` contains CAZy's classification of the proteins into the CAZy classes.

`class_predictions_classifications_{time_stamp}.csv` contains the prediction tool CAZy class classifications for each protein.

`class_stats_across_all_test_sets_{time_stamp}.csv` contains the statistical parameter values for each CAZy class across all test sets.

`class_stats_per_test_set_{time_stamp}.csv` contains the statistical parameter values for each CAZy class, for each test set. This is better for seeing the variation of prediction tool performance.


## [4] Presenting the evaluation

To visualise the data of the statistical evaluation, we use R. We specifically use an R notebook which can be knitted into HTML and can be easily shared and viewed. In this directory we provide the Rnotebook used in our evaluation, and encourage others to use it as a template it for their own evaluations, and modify as they see fit. If you see areas of improvement for the R notebook please raise an issue at this repository!

The Rnotebook and the data used to create the figures for the poster presented at the Microbiology Society conference are stored in the directory `report`.