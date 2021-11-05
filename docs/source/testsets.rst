
===============================================
Creating the test sets
===============================================
.. _test-set-label:

``pyrewton`` creates test sets from GenBank genomic assemblies. One test set is created per genomic assembly.

.. TIP::
    CAZy annotates protein sequences contained in the GenBank releases. Therefore, we recommend using 
    GenBank genomic assemblies instead of RefSeq assemblies. This is because it is possible some 
    protein sequences extracted from the genomic assemblies will have accession numbers not catalgoued 
    in CAZy and will be incorrectly identified as non-CAZymes by ``pyrewton``.

The Python scripts ``create_test_sets*.py`` are used to generate the test sets.

They are invoked using the same command structure:

.. code-block:: bash
    create_test_sets*.py \
    <user email> \
    <yaml file path> \
    <path to cazy json or db> \
    <path to output dir>

Two ``create_test_sets*.py`` scripts are provided.  
The script ``create_test_sets_from_db.py`` is used when retrieving CAZy family annotations from a 
local CAZyme database created using cazy_webscraper. The script ``create_test_sets_from_dict.py`` is 
used when retrieving CAZy family annotations from a JSON file keyed by GenBank accessions and 
valued by list of CAZy family annotations.

.. NOTE::
   ``cazy_webscraper`` is used to compile the local CAZyme database, and this must be compiled 
   before creating the test sets. Instructions for installing and using ``cazy_webscraper`` can 
   be found in:  
   - `GitPages <https://hobnobmancer.github.io/cazy_webscraper/>`_   
   - `pypi <https://pypi.org/project/cazy-webscraper/>`_   
   - `Read the Docts <https://cazy-webscraper.readthedocs.io/en/latest/?badge=latest>`_   
    in the GitHub repo, pypi and at Read the Docs.

.. TIP::
    For future evaluations we recommend using the local CAZyme database approach, because the 
    database includes logs of its creation, facilitating reproducibility of the evaluation.

Both ``create_test_sets*.py`` scripts take the same required and optional input arguments.  

**Required arguments:**
- User email address: Required by NCBI Entrez for downloading genomes to produce the test sets from
- YAML file path: Path to a YAMl file keyed by NCBI taxonomy IDs and valued by list of GenBank accessions, this defines the genomes to be downloaded and used to create the test sets
- Path to CAZy JSON or local CAZyme database db file
- Path to an output directory to write the output to - this output directory and its parents do not need to already exist

**Optional arguments:**
``--force``, ``-f`` - force overwriting existing output  
``--genomes``, ``-g`` - path to directory containing already downloaded genomes (saves redownloading the genomes if the ``pyrewton`` crashs when compiling the test sets)  
``--log``, ``-l`` - path to write out log file  
``--nodelete``, ``-n`` - do not delete existing output in existing output directory  
``--sample_size``, ``-s`` - number of CAZymes to be included in each test set. An equal number of non-CAZymes as CAZymes are automatically added to each test  
``--verbose``, ``-v`` - enable verbose logging

--------------------------
YAML of selected genomes
--------------------------

To tell ``pyrewton`` which genomes to download and compile test sets from, a YAML file is passed to 
``pyrewton``. This YAML file must have the structure of keyed by NCBI taxonomy IDs, and valued 
by list of GenBank accession data.

Specififically, each selected GenBank genomic assembly selected to be used to create a test set 
is represented as a list. The first element is the **genomic accession number**, the second element 
is the **assembly name**.

An example snippet of a YAML file is shown below:  
.. code-block:: yaml
    txid498019:
        - [GCA_003013715.2, ASM301371v2]
        - [GCA_008275145.1, ASM827514v1]
    txid573826:
        - [GCA_000026945.1, ASM2694v1]
    txid13502:
        - [GCA_011074865.2, ASM1107486v2]
    txid436017:
        - [GCA_000092065.1, ASM9206v1]

The YAML file used for the evaluation in Hobbs et al., 2021 can be found in the ``pyrewton`` 
`repo <https://github.com/HobnobMancer/pyrewton/tree/master/supplementary/mar_2021_eval>`_

-------------------------
Changing the sample size
-------------------------

The default test set size is 100 CAZymes and 100 non-CAZymes, which was used in Hobbs et al., 2021. 
To change the sample size, call the ``--sample_size`` flag and define the number of CAzymes to be 
included in the test set, an equal number of non-CAZymes wills be added to the test set. 
For example, ``--sample_size 200`` will produce test sets of 200 CAZymes and 200 non-CAZymes.

--------------------------
The output
--------------------------

In the output directory specified by the user, four directories and a plain text file are produced:

- ``alignment_scores``: contains ``.csv`` files of the BLAST all-versus-all scores of the selected CAZymes query against all non-CAZymes
- ``extracted_proteins_seqs``: FASTA files containing all extracted protein sequences from the downloaded genomic assemblies, one FASTA file is created per genome and contains all extracted protein sequences from that one genome
- ``genomes``: contains the downloaded genomic assemlbies in gbff.gz format
- ``test_sets``: contains the tests sets (FASTA) files for be used as input for each CAZyme classifier. One test set is created per genome.
- ``cazome_coverage_<time stamp>.txt``: Contains the following headers and data:
    - Genomic_accession: The accession of the genomic assembly
    - Total_proteins: Total number of proteins extracted from the genomic assembly
    - Total_CAZymes: Total number of CAZy annotated CAZymes extracted from the genomic assembly
    - Genome_CAZome_percentage: Percentage of the proteome contained within the CAZome
    - CAZome_coverage_percengate: Percentage of the identified CAZome included in the test set compiled from the genomic assembly
    - CAZyme_sample_size: Number of CAZymes included in the test set compiled from the genomic assembly
