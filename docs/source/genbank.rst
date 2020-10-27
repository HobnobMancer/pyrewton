
=========================
The genbank module
=========================

The module ``genbank`` contains submodules that handle GenBank files. This includes the retrieval of GenBank files from
the NCBI Assembly database, and retrieval of protein annotations (writing out th protein sequences to FASTA files and protein data to a summary dataframe) from GenBank files.

| - :ref:`genbank-mod-get-ncbi-genomes-def`
| - :ref:`genbank-mod-get-genbank-annotations-def`

.. _genbank-mod-get-ncbi-genomes-def:

get_ncbi_genomes
----------------

``get_ncbi_genomes`` is a submodule of the ``genbank``, it takes a plain text file containing species scientific names
or NCBI taxonomy as input. The script finds the corresponding taxonomy ID or scientific name, as appropriate, as well
as retrieve all directly linked accession numbers. The submodule generates a dataframe containing 'Genus', 'Species',
'NCBI Taxonomy ID', and 'NCBI Accession Number'.

Invoking the script from the command line
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To invoke ``get_ncbi_genomes`` invoke the script ``get_ncbi_genomes.py``, for example, using the following command:   
``python3 get_ncbi_genomes <user_email> <-i path_to_input_.txt> <-o directory_to_store_assemblies> <-d species_dataframe_output_path``


**Compulsary arguments**

``-u, --user`` - Although indicated as optional, Entrez requires an email address must be provided.
If not provided the programme will log this as an error and terminate.

**Optional arguments**

``-d, --`dataframe`` - Specify output path for dataframe, which will be saved as a .csv file
(inclusion of file extensions is optional). If not provided dataframe will be written out to STDOUT.

``-f, --force`` - Enable writting in specificed output directory if output directory already exists.

``-g, --genbank`` - Enable or disable downloading of GenBank files.

``-h, --help`` - Display help messages and exit, including listing arguments

``-i, --input`` - Specify input filename (with extension) input file.If only the filename is supplied
``get_ncbi_genomes.py`` will only look in the curent working directory, otherwise provide path to
input file. If no option is given the default input is taken from STDIN.

``-l, --log`` - Specify name of log file (With extension). If only filename is given, log file will 
be written out to the current working directory, otherwise provide path including filename. If not 
option is given no log file will be written out, however, logs will still be printed to the terminal.

``-n, --nodelete`` - Enable not deleting files in existing output directory. If not enabled, output 
directory exists and writing in output directory is 'forced' then files in output directory will not 
be deleted, and new files will be written to the output directory.

``-o, --output`` - Specify filename (with extension) of output file. If not option is given output 
will be written to STDOUT.

``-r, --retries`` - Specifiy maximum number of retries before cancelling call to NCBI
if a network error is encountered. The default is a maximum of 10 retries. When
maximum is reached, a value of 'NA' is returned.

``-t, --timeout`` - Specify timeout limit of URL connection when downloading GenBank files.
Default is 10 seconds.

``-v, --verbose`` - Enable verbose logging - changes logger level from WARNING to INFO.

.. _genbank-mod-get-genbank-annotations-def:

get_genbank_annotations
----------------------

``get_genbank_annotations`` is a submodule of the ``genbank``, which reads NCBI GenBank files and
retrieves annotations of encoded Proteins. Specifically, ``get_genbank_annotations`` parses each genomic assembly 
contained within the input directory passed the module, during which it identifies a protein by its annotation as a 'CDS' feature.
The summary data (listed out below) for every protein is written out to a dataframe, with a unique protein in each row. The sequence of each
protein is also written out to a FASTA file. Proteins from the same genomic assembly are written out to the same FASTA file, therefore, a singel
FASTA file containing the sequences of all the proteins identified within the assembly is created per genomic assembly in the input directory.

Dataframe columns:
* 'Genus'
* 'Species'
* 'NCBI Taxonomy ID'
* 'NCBI Accession Number'
* 'NCBI Protein ID'
* 'Locus Tag'
* 'Gene Locus'
* 'NCBI Recorded Function'
* 'Protein Sequence'
* 'UniProtKB Entry ID'
* 'UniProt Entry Name'
* 'UniProtKB Protein Names'
* 'EC number'
* 'Length (Aa)'
* 'Mass (Da)'
* 'Domains'
* 'Domain count'
* 'UniProtKB Linked Protein Families'
* 'Gene ontology IDs'
* 'Gene ontology (molecular function)'
* 'Gene ontology (biological process)'

``get_genbank_annotations`` has been created as a subsequent or complementary module to ``get_ncbi_genomes`` and is designed to parse the outputs from ``get_ncbi_genomes``.

Invoking the script from the command line
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To invoke ``get_genbank_annotations`` invoke the script ``get_genbank_annotations.py``, for example, use the following command:   
``python3 get_genbank_annotations <path_to_input_df> <path_to_assemly_dir>``

**Compulsary arguments**

``--input_df`` - Path to the dataframe containing the accession numbers of the genomic assemblies. This dataframe is created by ``get_ncbi_genomes``.

``--genbank`` - Path to directory containing the genomic assemblies listed within the dataframe. Unless they have been moved, this is the output directory from ``get_ncbi_genomes``.

**Optional arguments**

``-h, --help`` - Display help messages and exit, including listing arguments

``-f, --force`` - Enable writting in specificed output directory if output directory already exists.

``-l, --log`` - Specify name of log file (With extension). If only filename is given, log file will 
be written out to the current working directory, otherwise provide path including filename. If not 
option is given no log file will be written out, however, logs will still be printed to the terminal.

``-n, --nodelete`` - Enable not deleting files in existing output directory. If not enabled, output 
directory exists and writing in output directory is 'forced' then files in output directory will not 
be deleted, and new files will be written to the output directory.

``-o, --output`` - Specify filename (with extension) of output file. If not option is given output 
will be written to STDOUT.

``-v, --verbose`` - Enable verbose logging - changes logger level from WARNING to INFO.

.. note::
    get_genbank_annotations is still under development.
    Please see the `GitHub repository <https://github.com/HobnobMancer/pyrewton/tree/master>`_ for the latest developments.
