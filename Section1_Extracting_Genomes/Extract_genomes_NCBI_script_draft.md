# Drafting out the script Extract_genomes_NCBI.py

Drafting out the script Extract_genomes_NCBI.py involves stating the overall aim of the script; how the script will be broken up into sections/overall flow/structure of the script; and planning out what functions (including what each function will aim to do).

??? indicates when there is an idea and not sure if it is a good idea
?_P_? indicates an idea but not sure it is in the correct place in the script

## Aims of Extract_genomes_NCBI.py

> LP: The first aim in this list should be to process a list of desired species/taxIDs
> LP: we should try to refer to the "numbers" as "IDs" or "accessions"

1. Extract taxonomic numbers from NCBI Taxonomy database
2. Use taxonomic numbers to pull down associated accession numbers from NCBI
3. Create a table of genus / species / taxonomic ID / accession number
4. Use the accession numbers to pull down genomic assemblies for each species
5. Write out the genomic assemblies with a standard file-naming system

### Section 1: Extract taxonomic numbers

> LP: in (2) I don't think you need to store the names - they're already stored in the file that's being read; also, if you use a Pandas dataframe you don't need to use a dictionary (the dataframe index plays the part of a dictionary key)
> LP: in (3) working with dataframes will suggest a natural way of iterating over the rows (one row per species/taxID)
> LP: in practice (5) will be done as part of (4): the link with species name is implied in the search, and finding the taxID is the same thing as "pulling it down"

1. Read the species list plain text file
2. Extract the species names and store in format that can be read/passed to Entrez
??? Pass species names to a dictionary so that the genus is the key and species is the value, or to a matrix, so that the names can be easily added to the output table in overall script aim 3.
3. Create a loop that will iterate over the names list (clearly going to need to also pass the names from the species list file to list then) and pass them to Entrez
4. Use Entrez search to find Taxonomic IDs
5. Pull down Taxonomic IDs and link it to the species name
    5a. This is starting to look like I will create a matrix, which I will slowly build up/add to

### Section 2: Pull down accession numbers

> LP: in (1) have a look at Pandas dataframes, I think they're what you think of as a "matrix"

1. Pass the taxnomic numbers to Entrez eLink
    1a. will need to iterate over a list of some sort so that each number/ID can be handled at a time
    1b. would a matrix faciltiate this? It's mutlple aligned arrays so it should do - look this up, especially the necessary syntax
2. Pull down the accession numbers found by Entrez eLink
    2a. store the accession numbers so that are associated with the taxonomic IDs
    2b. does sound like a matrix

### Section 3: Create table

If I build a matrix as I go along then this should be created already at this point and it will simply be the task at out putting the table in a human-readable format.

### Section 4: Use accession numbers to pull down assemblies

This will link with section 5 so that as soon as the assembly is pulled down it is written out to a file using a standard file format.

This will invovle using Entrez eFetch to pull down the NCBI Assembly database.

