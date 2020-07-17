
=============================
The annotations module
=============================

The module ``annotations`` contains submodules that handle cazyme annotations, including summarising the cazyme calls
distribution in a dataframe of cazymes and the retrieval of cazymes from UniProtKB for a given species.

| - :ref:`annotations-mod-get-uniprot-proteins`

.. _annotations-mod-get-uniprot-proteins:

get_uniprot_proteins
--------------------

``get_uniprot_proteins`` is a submodule of the ``annotations``, and retrieves all cazyme entries from UniProtKB
for each species passed to the submodule in a dataframe. The collected data is returned in a dataframe containing:

* 'Genus'
* 'Species'
* 'NCBI Taxonomy ID'
* 'UniProt entry ID'
* 'UniProt entry name'
* 'UniProt assigned protein names'
* 'EC Numbers'
* 'Length (Aa)'
* 'Mass (Da)'
* 'Domains'
* 'Domain count'
* 'UniProt linked protein families'
* 'GO IDs'
* 'GO molecular function'
* 'G0 biological process'
* 'Sequence'

Invoking the script from the command line
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Compulsary arguments**

``-d, --df_input`` - Path to input dataframe.

``-f, --force`` - Enable writting in specificed output directory if output directory already exists.

``-g, --genbank`` - Path to director containing GenBank files.

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
    get_uniprot_proteins is still under development.
    Please see the `GitHub repository <https://github.com/HobnobMancer/PhD_Project_Scripts/tree/master>`_ for the latest developments.
