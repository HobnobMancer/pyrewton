
================
Trouble Shooting
================

This section covers common errors expected to arise when invoking each module/submodule, and the probable causes.
If any other issues arise, please raise any issues with any of the programmers at the GitHub repository issues pages, by
following `the link <https://github.com/HobnobMancer/PhD_Project_Scripts/issues>`_.

| - :ref:`genbank-mod-get-ncbi-genomes`
| - :ref:`genbank-mod-get-cazyme-annotations`

genbank: get_nbci_genomes
=========================

.. _genbank-mod-get-ncbi-genomes:

This section deals with troubleshooting the genbank module's submodule get_ncbi_genomes,
for the retrieval of GenBank files from the NCBI Assembly database.

.. warning::
    The majority of issues will arise due to errors in the input file. Always ensure the input file does not contain any
    error or blank lines. Additionally, always ensure the correct path to the input file is provided.

**IOError**
This error will occur if there is a network issue when using Entrez to call to NCBI. The script will automatically
retry the call the set maximum number of times. If the maximum number of retries is met before connecting to NCBI without
encountering a network error, 'NA' is returned and stored in the dataframe.

**FileNotFoundError**
This error will occur is the incorrect path is provided as the input argument at the command line, or no input argument
is provided and STDIN contains no data. Ensure the path includes the file name, with extension. If this error occurs the
program will terminate.

**IndexError during scientific name retrieval**
This occurs when Entrez fails to retrieve a scientific name for the given taxonomy ID. This is potentially caused by a typo
in the taxonomy id provided in the input file. If this error occurs the string 'NA' will be returned.

**IndexError during taxonomy ID retrieval**
This occurs when Entrez fails to retrieve a taxonomy ID for the given scientific name. Returns 'NA'. This is potentially
caused by a typo in the species name in the input file, or a typo in a taxonomy ID 'NCBI:txid' prefix, causing the program
to misinterpret the ID as a species name and use it to try and retrieve a scientific name. If this error occurs the string
'NA' will be returned. If no taxonomy ID is available for the retrieval of accession numbers, the retrieval of accession
numbers is cancelled, and a value of 'NA' is returned.

**IndexError during assembly ID retrieval**
This occurs when Entrez fails to retrieve assembly IDs from NCBI. This may be because there are no directly linked assemblies
for the given taxonomy ID. Check the NCBI Taxonomy database to ensure there are 'directly' linked assemblies and not only
'subtree' assemblies. If this error occurs the program with exit the retrieve of the assembly IDs and not retrieve the
NCBI accession numbers, and return the string 'NA'. This allows for troubleshooting using on the specie(s) for which it is
required, to reduce demand on NCBI.

**RuntimeError during posting of assembly IDs to NCBI**
This error occurs when Entrez fails to post the retrieved assembly IDs, causing it to fail to retrieve the document summary.
This is potentially caused by incorrect formatting of the assembly IDs or the request is too large for Entrez/ NCBI.
If this is the case, repeat the procedure in batches. If this error occurs the program with exit the posting of the assembly
IDs and not retrieve the NCBI accession numbers, and return the string 'NA'. This allows for troubleshooting using on the
specie(s) for which it is required, to reduce demand on NCBI.

**RuntimeError or IndexError during retrieval of accession numbers**
This error occurs when Entrez fails to retrieve the document summary from NCBI. This is potentially caused by incorrect
formatting of the assembly IDs or the request is too large for Entrez/ NCBI. If this is the case, repeat the procedure
in batches. If this error occurs the program with exit the retrieval of the NCBI accession numbers, and return the
string 'NA'. This allows for troubleshooting using on the specie(s) for which it is required, to reduce demand on NCBI.

.. _genbank-mod-get-cazyme-annotations:

get_cazyme_annotations
======================

This section deals with troubleshooting the genbank module's submodule get_cazyme_annotations,
for retrieving cazyme annotations from GenBank files.
