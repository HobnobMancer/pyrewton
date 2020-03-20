# Drafting out the script Extract_genomes_NCBI.py

Drafting out the script Extract_genomes_NCBI.py involves stating the overall aim of the script; how the script will be broken up into sections/overall flow/structure of the script; and planning out what functions (including what each function will aim to do).

??? indicates when there is an idea and not sure if it is a good idea.
?_P_? indicates an idea but not sure it is in the correct place in the script.

## Aims of Extract_genomes_NCBI.py

1. Extract names of species from input file
2. Extract taxonomic numbers from NCBI Taxnomy database
3. Use taxonomic numbers to pull down associated accession numbers from NCBI
4. Create a table of genus / species / taxonomic ID / accession number
5. Use the accession numbers to pull down genomic assemblies for each species
6. Write out the genomic assemblies with a standard file-naming system

### Section 1: Extract names from Species_list.txt

**This will form the first function that will be implemented, called 'get_names_and_TaxIDs'.**

1. Use 'with open()' to open the plain text file containing the list of selected species: Species_list.txt
2. Species_list.txt will contain gensus-species names, comments and potentially taxanomic IDs as well; therefore, will need to distguish between species names and taxonomic IDs. Species names passed to one list and taxanomic names passed to another list.
3. Taxanomic IDs will start with 'NCBI:txid', therefore, create a for loop that will pass over the species_list and will check for the following conditions:
a. Is the first character '#' - comments will be distinguished by starting with '#' - if 'yes' do nothing, if 'no' move onto the next condition - might be best to use "if item[0] is not '#':" then nest the next two conditions within
b. Are the first nine characters 'NCBI:txid', if 'yes' pass item to taxonomic_id_list
c. If the first nine characters are not 'NCBI:txid' then pass the item to the genus_species_name_list

Example code structure to test first character to determine if taxonomic ID, species name of comment:
    for entry in input_list:
        if entry[0] is not '#':
            if entry[0:10] is 'NCBI:txid':
                taxonomic_id_list.append(entry)
            else:
                genus_species_name_list.append(entry)

**The lists genus_species_name_list and taxanomic_IDs_list will be defined as global lists**

### Section 2: Extract taxonomic numbers

1. Create a loop that will iterate over the genus_species_name_list and pass them to Entrez
2. Use Entrez search to find Taxonomic IDs
3. Pull down Taxonomic IDs and link it to the species name. This is starting to look like I will create a matrix, which I will slowly build up/add to

### Section 3: Pull down accession numbers

> LP: in (1) have a look at Pandas dataframes, I think they're what you think of as a "matrix"

1. Pass the taxnomic numbers to Entrez eLink
    1a. will need to iterate over a list of some sort so that each number/ID can be handled at a time
    1b. would a matrix faciltiate this? It's mutlple aligned arrays so it should do - look this up, especially the necessary syntax
2. Pull down the accession numbers found by Entrez eLink
    2a. store the accession numbers so that are associated with the taxonomic IDs
    2b. does sound like a matrix

### Section 4: Create table

1. Create empty dataframe
species_dataframe = pd.DataFrame({})
2. Create two empty lists for genus and species names
genus_list = []
species_list = []
3. Separate genus and species names:
for name in genus_species_name_list:
    genus_species = name.split(' ')
    genus_list.append(genus_species[0])
    species_list.append(genus_species[1])
4. Add columns with associated data to the dataframe
species_dataframe['Genus'] = genus_list
species_dataframe['Species'] = species_list
species_dataframe['Taxonomic ID'] = taxonomic_id_list
species_dataframe['Accession number'] = accession_number_list

And here lies the problem. This works if the input species list only contains genus-species names, however if the list includes taxonomic numbers they will be put into the table first, placing the taxonomic IDs out of place with their respective genus-species names. Also no genus-species name will be recorded for the taxonomic IDs taken from the input file.

Potentially, the taxonomic IDs taken from the inpu could be placed in one list, the taxonomic IDs pulled down placed in another and the taxonomic IDs taken from the input file then added to the list of taxonomic IDs that were pulled down. This would place the taxonomic IDs from the input file at the end of the list and thus not misalign the taxonomic IDs from the pull down from their associated genus-species names.

For the issue with the taxonomic IDs taken from the input file can be solved by either manual addition of the names afterwards, which mitigates the point of inputting taxonomic IDs to begin with. Alteranatively, state that only species names can be given and not taxonomic IDs, or use Entrez to pull down the associated species name and add this to the end of the genus_species_name_list. If using the method directly above this to add the input file taxonomic IDs to the end of the taxonomic IDs pulled down list then the associated species names of the taxonomic IDs from the input will be aligned in the dataframe.

Alternatively to the dataframe being created so late, a dataframe with the genus, species and taxonomic IDs could be created before pulling down the accession numbers. Then the taxonomic column can be passed to entrez, and then the list off accession numbers added to the dataframe as a new column.

### Section 5: Use accession numbers to pull down assemblies

This will link with section 5 so that as soon as the assembly is pulled down it is written out to a file using a standard file format.

This will invovle using Entrez eFetch to pull down the NCBI Assembly database.
