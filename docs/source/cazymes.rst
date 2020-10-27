
=========================
The cazymes module
=========================

The module ``cazymes`` contains submodules that automate the retrieval of protein sequences and associated annotations from UniProtKB,
invokes the CAZyme prediction tools dbCAN, CUPP and eCAMI to predict which proteins within a FASTA file are CAZymes, and evaluates
the performance of the CAZymes prediction tools.

| - :ref:`uniprot-def`
| - :ref:`genbank-mod-get-genbank-annotations-def`

.. _uniprot-def:

uniprot
----------------

The ``uniprot`` submodule automates the querying and parsing of data from [UniProtKB](https://www.uniprot.org/). The querying of UniProtKB is 
configuerable via a yaml file. For every protein retrieved from UniProtKB the protein sequence is written to a FASTA file (with a single FASTA file per query to UniProtKB), and the 
following data for each protein is written out to a summary dataframe (writing a dataframe per UniProt query):

* ‘NCBI Taxonomy ID’
* ‘Organism’
* ‘UniProtKB Entry ID’
* ‘UniProtKB Entry Name’
* ‘UniProtKB Protein Names’
* ‘EC number’
* ‘Length (Aa)’
* ‘Mass (Da)’
* ‘Domains’
* ‘Domain count’
* ‘UniProtKB Linked Protein Families’
* ‘Gene ontology IDs’
* ‘Gene ontology (molecular function)’
* ‘Gene ontology (biological process)’
* ‘Sequence’

**Configuring the querying of UniProt**

An example configuration yaml file is located within the directory pyrewton/cazymes/uniprot, called 'uniprot_cazy_ec_retrieval.yaml'.
Specific queries can be written under the 'queries' tag in a yaml file. These queries must be written in the [UniProt syntax](https://www.uniprot.org/help/text-search).
A list of UniProt query fields is available [here](https://www.uniprot.org/help/query-fields).

To restrict these queries to a specific set of candidate species add the NCBI taxonomy ID of the species under the 'tax_id' in the yaml file.
Specifically, this causes every query that is listed under the 'queries' tag to be performed for each species listed under the 'tax_id' tag,
saving time having to write out the same set of queries multiple times over with a different taxonomy ID. If no taxonomy ID is given under the 'tax_id' tag, 
and no taxonomy ID is included within the query under the 'queries' tag the search will no be restricted to a specific species.

Providing specific queries under the 'queries' tag of the configuration yaml file is also optional. If only the taxonomy ID of a species is given
and there are no queries under the 'queries' tag, all proteins for that given species will be retrieved from UniProtKB.

**If you are looking to retrieve all potential CAZymes from UniProtKB for a set of candidate species** the example yaml configuration file 'uniprot_cazy_ec_retrieval.yaml'
contains EC numbers retrieved from the [CAZy database](http://www.cazy.org/) and are strongly associated with CAZyme activity, and also contains a query
to retrieve all proteins that are catalogued in both UniProtKB and CAZy. Therefore, it is suggested this file is used as the configuration file by either changing the taxonomy
id to the your candidate species NCBI taxonomy IDs or removing all taxonomy IDs and retrieving all proteins which are linked to CAZy or annotated with one of the specified EC numbers
from UniProtKB.

Invoking the script from the command line
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

An example of basic operation is: 
``python3 get_uniprot_proteins <path_to_config_file.yaml>``


**Compulsary arguments**

``-u, --user`` - Although indicated as optional, Entrez requires an email address must be provided.
If not provided the programme will log this as an error and terminate.

**Optional arguments**

``-a, --fasta`` - Enable writing out FASTA files containing protein sequences

``-f, --force`` - Enable writting in specificed output directory if output directory already exists.

``-h, --help`` - Display help messages and exit, including listing arguments

``-l, --log`` - Specify name of log file (With extension). If only filename is given, log file will 
be written out to the current working directory, otherwise provide path including filename. If not 
option is given no log file will be written out, however, logs will still be printed to the terminal.

``-n, --nodelete`` - Enable not deleting files in existing output directory. If not enabled, output 
directory exists and writing in output directory is 'forced' then files in output directory will not 
be deleted, and new files will be written to the output directory.

``-o, --output`` - Specify output directory where all output dataframes and FASTA files are written to If not option is given output 
will be written to STDOUT.

``-v, --verbose`` - Enable verbose logging - changes logger level from WARNING to INFO.


.. _prediction-def:

prediction
----------------------

The ``prediction`` module is for the prediction of CAZymes and non-CAZymes from protein sequences within FASTA files. The module
is setup to take the output from from the ``uniprot`` and ``get_genbank_annotations`` submodules, although it will accept any FASTA file containing
protein sequences.

This module invokves three third-party CAZyme prediction tools that predict if a given protein sequence is a CAZyme or non-CAZyme, and then predicts 
the CAZy family of any predicted CAZymes. These tools include [dbCAN](https://github.com/linnabrown/run_dbcan), [CUPP](https://www.bioengineering.dtu.dk/english/researchny/research-sections/section-for-protein-chemistry-and-enzyme-technology/enzyme-technology/cupp), and [eCAMI](https://github.com/zhanglabNKU/eCAMI).

This module also evaluates the overall performance of each of the prediction tools, but also evaluates the performance of each prediction tool per input FASTA file. The latter allows the evaluation of each
prediction tool per candidate species if a each input FASTA file contains proteins from a single species/genomic assembly. In order to perform this evaluation a gold standard is required. CAZy is the gold standard
database for all identified and catalogued CAZymes, but does not offer a method of automated protein data retrieval, therefore, use this [CAZy webscrapper](https://github.com/HobnobMancer/cazy_webscraper) to automate
the retrieval of the necessary protein data from CAZy. Using this web scraper will ensure the data retrieved from CAZy is in the correct format for the statistical evaluation of the performnce of the dbCAN, CUPP and eCAMI.

**The output**

For every input FASTA file a respective output directory is created within the user specified output directory. In this way the output directory specified by the user becomes the parent directory of all output directories from the module.
Within each input FASTA files output directory there will be:

* The output directory from dbCAN (becuase dbCAN creates its own output directory)
* A fasta and log output files from CUPP
* An output .txt file from eCAMI
* A standardised output dataframe for each prediction tool (therefore, one output dataframe for dbCAN, DIAMOND, HMMER, Hotpep (which are invokved within dbCAN), CUPP and eCAMI)
* Results of the statistical evaluation of prediction tools performance
* Report presenting the evaluation of the prediction tools performance

The parent output directory (the output directory specified by the user) will also contain the results and report of the statistical evaluation of the prediction tools' overall performance across all input FASTA files passed to the module when the module was invoked.


Invoking the script from the command line
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This script is still under development.

.. note::
    get_genbank_annotations is still under development.
    Please see the `GitHub repository <https://github.com/HobnobMancer/pyrewton/tree/master>`_ for the latest developments.
