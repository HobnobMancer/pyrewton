# Drafting out the script `Extract_genomes_NCBI.py`

Drafting out the script Extract_genomes_NCBI.py involves stating the overall aim of the script; how the script will be broken up into sections/overall flow/structure of the script; and planning out what functions (including what each function will aim to do).

> LP: you can make your markers clearer by formatting them. I've done this below by adding backticks around them. Markdown interprets this as "code", and it's a neat way to indicate some kind of marker that should be understood as an intrusion into the flow of text.
>
> LP: it is best practice to use these backticks to indicate filenames, too. I've done this above and below, but you should handle the rest of the document.

`???` indicates when there is an idea and not sure if it is a good idea.
`?_P_?` indicates an idea but not sure it is in the correct place in the script.

## Aims of `Extract_genomes_NCBI.py`

1. Extract names of species from input file
2. Extract taxonomic numbers from NCBI Taxnomy database
3. Use taxonomic numbers to pull down associated accession numbers from NCBI
4. Create a table of genus / species / taxonomic ID / accession number
5. Use the accession numbers to pull down genomic assemblies for each species
6. Write out the genomic assemblies with a standard file-naming system

### Section 1: Extract names from Species_list.txt

**This will form the first function that will be implemented, called 'get_names_and_TaxIDs'.**

> LP: you can use the backticks to format code and filenames inline, within a sentence, too (as below)

1. Use `with open()` to open the plain text file containing the list of selected species: `Species_list.txt`
2. Species_list.txt will contain genus-species names, comments and potentially taxonomic IDs as well; therefore, will need to distinguish between species names and taxonomic IDs. Species names passed to one list and taxonomic names passed to another list.
3. Taxonomic IDs will start with 'NCBI:txid', therefore, create a for loop that will pass over the species_list and will check for the following conditions:
a. Is the first character '#' - comments will be distinguished by starting with '#' - if 'yes' do nothing, if 'no' move onto the next condition - might be best to use "if item[0] is not '#':" then nest the next two conditions within
b. Are the first nine characters 'NCBI:txid', if 'yes' pass item to taxonomic_id_list
c. If the first nine characters are not 'NCBI:txid' then pass the item to the genus_species_name_list

> LP: in the above, I think you won't need to maintain two lists, but you will need to process "<genus> <species>" differently to taxIDs.
>
> LP: what I mean by this is that, when you read a line, you will need to process *either* a "<genus> <species>" string, *or* a taxID. In either case, you will end up with an association between "<genus> <species>" and taxID. Both processing actions will return *similar* output for different input - they *converge*. Because they converge, you can treat the output of the functions as the same. That is, you can maintain a single list.
>
> LP: the way I think of this is as a single input, which is passed down one of two pathways (i.e. into one of two functions: `process_taxid()` or `process_name()`), to get the same output at the end: `(taxID, genus, species)`.

Example code structure to test first character to determine if taxonomic ID, species name of comment:

> LP: Markdown allows you to use "fenced blocks" for code. These are surrounded (start/end) by three backticks, as below. You can also name the language that's being written, so that it gets *syntax highlighting* - keywords are coloured differently from strings, etc. (to ease reading). VSCode is sensitive to this when it knows you're editing Markdown and will also syntax highlight live in the editor.

```python
for entry in input_list:
    if entry[0] is not '#':
        if entry[0:10] is 'NCBI:txid':
            taxonomic_id_list.append(entry)
        else:
            genus_species_name_list.append(entry)
```

> LP: The overall logic is sound. The think I'd change is where you are appending directly to two separate lists, instead call a function that returns `(taxID, genus, species)` and append that to a single common list. [NOTE: later you may want to change this to use *generators*, but we'll save that technical detail for now. Remind me later, though!]
>
> LP: The test for `NCBI:txid` will work, and it will be fast enough. I'd stil encourage you to use a *regular expression*, though. Using a regular expression will allow you to check the line for valid formatting, so that you can process something like `NCBI:txid469` or `NCBI:txid836475`, but discard bobbins like `NCBI:txidHaThisIsntARealTaxID`, without processing them. Regular expressions are incredibly useful, and crop up all over the place.

**The lists genus_species_name_list and taxanomic_IDs_list will be defined as global lists**

> LP: global variables are a "code smell" - they can indicate poorly-structured code. There are certainly times to use them, but this is not one of them ;)

### Section 2: Extract taxonomic numbers

1. Create a loop that will iterate over the genus_species_name_list and pass them to Entrez
2. Use Entrez search to find Taxonomic IDs
3. Pull down Taxonomic IDs and link it to the species name. This is starting to look like I will create a matrix, which I will slowly build up/add to

> LP: this could be broken out into a function (e.g. `process_name()` that takes a "<genus> <species>" name and returns `(taxid, genus, species)`)
>
> LP: having another function (e.g. `process_taxid()` that takes a taxID and returns the same data structure) would be useful.

### Section 3: Pull down accession numbers

> LP: in (1) have a look at Pandas dataframes, I think they're what you think of as a "matrix"

1. Pass the taxonomic numbers to Entrez eLink
    1a. will need to iterate over a list of some sort so that each number/ID can be handled at a time
    1b. would a matrix faciltiate this? It's mutlple aligned arrays so it should do - look this up, especially the necessary syntax
2. Pull down the accession numbers found by Entrez eLink
    2a. store the accession numbers so that are associated with the taxonomic IDs
    2b. does sound like a matrix

> LP: For (1) above, you can use `Entrez.efetch()` directly with the taxID as the database identifier to get genus/species names, as below. You could use this for the accession numbers, too.

```python
>>> with Entrez.efetch(db="Taxonomy", id="148305") as handle:
...     record = Entrez.read(handle)
... 
>>> record
[{'TaxId': '148305', 'ScientificName': 'Pyricularia grisea', 'OtherNames': {'GenbankSynonym': ['Magnaporthe grisea'], 'Misspelling': [], 'Name': [{'ClassCDE': 'authority', 'DispName': 'Magnaporthe grisea (T.T. Hebert) M.E. Barr, 1977'}, {'ClassCDE': 'authority', 'DispName': 'Pyricularia grisea Cooke ex Sacc., 1880'}, {'ClassCDE': 'type material', 'DispName': 'CBS 138707'}, {'ClassCDE': 'type material', 'DispName': 'CBS:138707'}, {'ClassCDE': 'misspelling', 'DispName': 'Magnaportha grisea'}], 'EquivalentName': [], 'Misnomer': [], 'GenbankAnamorph': [], 'Inpart': [], 'Synonym': [], 'CommonName': [], 'Teleomorph': [], 'Includes': [], 'Acronym': [], 'Anamorph': []}, 'ParentTaxId': '48558', 'Rank': 'species', 'Division': 'Plants and Fungi', 'GeneticCode': {'GCId': '1', 'GCName': 'Standard'}, 'MitoGeneticCode': {'MGCId': '4', 'MGCName': 'Mold Mitochondrial; Protozoan Mitochondrial; Coelenterate Mitochondrial; Mycoplasma; Spiroplasma'}, 'Lineage': 'cellular organisms; Eukaryota; Opisthokonta; Fungi; Dikarya; Ascomycota; saccharomyceta; Pezizomycotina; leotiomyceta; sordariomyceta; Sordariomycetes; Sordariomycetidae; Magnaporthales; Pyriculariaceae; Pyricularia', 'LineageEx': [{'TaxId': '131567', 'ScientificName': 'cellular organisms', 'Rank': 'no rank'}, {'TaxId': '2759', 'ScientificName': 'Eukaryota', 'Rank': 'superkingdom'}, {'TaxId': '33154', 'ScientificName': 'Opisthokonta', 'Rank': 'no rank'}, {'TaxId': '4751', 'ScientificName': 'Fungi', 'Rank': 'kingdom'}, {'TaxId': '451864', 'ScientificName': 'Dikarya', 'Rank': 'subkingdom'}, {'TaxId': '4890', 'ScientificName': 'Ascomycota', 'Rank': 'phylum'}, {'TaxId': '716545', 'ScientificName': 'saccharomyceta', 'Rank': 'no rank'}, {'TaxId': '147538', 'ScientificName': 'Pezizomycotina', 'Rank': 'subphylum'}, {'TaxId': '716546', 'ScientificName': 'leotiomyceta', 'Rank': 'no rank'}, {'TaxId': '715989', 'ScientificName': 'sordariomyceta', 'Rank': 'no rank'}, {'TaxId': '147550', 'ScientificName': 'Sordariomycetes', 'Rank': 'class'}, {'TaxId': '222544', 'ScientificName': 'Sordariomycetidae', 'Rank': 'subclass'}, {'TaxId': '639021', 'ScientificName': 'Magnaporthales', 'Rank': 'order'}, {'TaxId': '2528436', 'ScientificName': 'Pyriculariaceae', 'Rank': 'family'}, {'TaxId': '48558', 'ScientificName': 'Pyricularia', 'Rank': 'genus'}], 'CreateDate': '2001/01/09 13:56:00', 'UpdateDate': '2019/03/14 17:03:22', 'PubDate': '1999/03/08 10:57:00'}]
```

> LP: In (1a) above, you should be able to pass a list of search terms. See, e.g. [this part of the Biopython tutorial](https://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc115)
>
> LP: for (1b), think of the dataframe as storing your data; you could certainly generate a query list from a dataframe column

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

> LP: See above. Given a taxID you can still recover the genus and species names. You can have two functions that take different arguments, and return the same output `(taxid, genus, species)`. These would be a common input for generating the dataframe.

Potentially, the taxonomic IDs taken from the input could be placed in one list, the taxonomic IDs pulled down placed in another and the taxonomic IDs taken from the input file then added to the list of taxonomic IDs that were pulled down. This would place the taxonomic IDs from the input file at the end of the list and thus not misalign the taxonomic IDs from the pull down from their associated genus-species names.

For the issue with the taxonomic IDs taken from the input file can be solved by either manual addition of the names afterwards, which mitigates the point of inputting taxonomic IDs to begin with. Alteranatively, state that only species names can be given and not taxonomic IDs, or use Entrez to pull down the associated species name and add this to the end of the genus_species_name_list. If using the method directly above this to add the input file taxonomic IDs to the end of the taxonomic IDs pulled down list then the associated species names of the taxonomic IDs from the input will be aligned in the dataframe.

Alternatively to the dataframe being created so late, a dataframe with the genus, species and taxonomic IDs could be created before pulling down the accession numbers. Then the taxonomic column can be passed to entrez, and then the list off accession numbers added to the dataframe as a new column.

### Section 5: Use accession numbers to pull down assemblies

This will link with section 5 so that as soon as the assembly is pulled down it is written out to a file using a standard file format.

This will invovle using Entrez eFetch to pull down the NCBI Assembly database.
