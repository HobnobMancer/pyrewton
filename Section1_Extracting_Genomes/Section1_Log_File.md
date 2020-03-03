# Section 1: Extracting Genomes from NCBI
## Introduction to the notebook

### Overview of 'Extract Genomes'
Genomic sequences with and without annotations are stored in the NCBI database. The aim of this section of the project is to identify which fungi and omycete species will be most relevant to pull down for later examination to identify carbohydrate processing enzymes. Simultaneously, a Python script will be written to pull down the identified genomes from the NCBI database using Entrez from the BioPython module, in an automated manner. A test sample will be created to work on and test the script for pulling down genomes from NCBI.

During this time, time will be dedicated to examinging the different output files of the NCBI database and identifying the taxonomic numbers of the desired fungal and omycete species.

### Notes on pulling down genomes from NCBI
Some will have annotations some will not, those that do not will be annotated using computational methods in the next section of the project.

Each entity in the NCBI database could have two sequences: on a GeneBank format and the other a RefSeq. It is advisable to find take the former format

There was a planning session with Leighton on 28-02-2020, setting out a broad prelimnary overview of the computational section of the PhD project. From this an overflow flow diagram was created which will be updated with details and corrected path as time progresses.

Tasks:

- [X] 1) Create flow chart of overal computational plan and upload to github

- [ ] 2) Identify plant pathogenic and/or plant degrading fungal and omycete species
  - [ ] i) Look in the literature to identify suitable species
  - [ ] ii) Identify and store the taxanomic NCBI identifiers
  
- [ ] 3) Spend time looking around the NCBI database
  - [ ] i) What are the different formats of the files available?
  - [ ] ii) Which output file formats are going to be most relevant and useful to pull down?
  - [ ] iii) What are the differences between RefSeq and GeneBank assembles?
  
- [ ] 4) Write a short 10-15 line Python script for pulling down genomes from NCBI
  - [ ] i) Spend some time reading the documentation for Entrez.BioPython
  - [ ] ii) Spend some time practising using Entrez.BioPython
  - [ ] iii) Create a test dataset for using to test any code
  - [ ] iv) Write the script and upload to GitHub then create a pullrequest for Leighton to see and check
    
Task 1 completed on 28-02-2020, and uploaded to main area of the computational repository on GitHub

### Papers for identifying fungal and omycete speices most suitable for future analysis
https://www.hindawi.com/journals/ijg/2018/1974151/
https://www.sciencedirect.com/science/article/pii/S1749461318300289
https://www.ncbi.nlm.nih.gov/pubmed/25192611
https://www.ncbi.nlm.nih.gov/pubmed/30825514
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3852112/pdf/1471-2164-14-S5-S7.pdf -- similar idea to the computational part of the PhD project
https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-019-6358-x - Ascomycota focused
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4577117/ - Phytophthora focused

(Project was advertised as looking at fungal Rhynchosporium and Magnaporthegenera and oomycete  Phytophthora genus, get evidence to show that these should be included)

### Working list of fungi and omycete identified for pull down from NCBI
Trichoderma genus (38 assemblies on NCBI) - lit is riff of examples of their lignocelluluisc catalytic efficiency, limitation of being a well studied genus although studies seem to be limited to a few species (see Horta et al., 2018) - should I take all or select? Seeing as so well studied a few maybe take all? Seeing as a small collection are well studied then won't find many novel from them, hence take all as may find a few novel from all 38



## Working log

### 29-02-2020 Notes when looking at the NCBI fungal taxnomic database
Many of the fungal entries are not classified.
Papers from above indicate that root associated fungi typically possess fewew PCWDEs, therefore may not need to include them???

Trichoderma are ascomycota, and share the same species path as aspergillus (that are also commonly studied and used for their carbohydrate processing enzymes) to leotiomyceta -- should more focus be but onto one of the common taxonomy layers of Trichoderma and aspergillus? - Dikarya, Ascomycota, saccharomyceta, Pezizomycotina, Leotiomyceta
-- Penicillium is another set of well studied species for biomass degrdatation (A, P and T from Weinmann et al., 2013) and also from Dikarya, Ascomycota, saccharomyceta, Pezizomycotina, Leotiomyceta

rhynchosporium, madnaporthe are also in Dikarya, Ascomycota, saccharomyceta, Pezizomycotina, Leotiomyceta
Therefore there is definetly something about the : Dikarya, Ascomycota, saccharomyceta, Pezizomycotina, Leotiomyceta lineage: narrows it down to 1105 genomes / 2731 assemblies out of 2269 / 5721 assemblies

A lot of focus has been put on Ascomycota for industrially used fungi

There are 326 aspergillus assemblies - surely I can't take all of them. Maybe take some well characterised so as to included estabilished annotations to help with computationally identification of PCWDE and some of the most recently sequenced becuase these are less likely to be well annotated and thus more likely to find novel enzymes

Common taxanomy:
- [X] K
- [X] P
- [ ] C
- [ ] O
- [ ] F
- [ ] G
- [ ] S

Kingdom: Fungi
  Subkingdom: Dikarya
Phylum: Ascinycota
   No rank: saccharomyceta
   Subphylum: Pezizomycotina
   No rank: leoriomyceta

#### Table of taxanomy of fungal species of potential interest
|  No rank |   sordariomyceta  |   sordariomyceta  |   sordariomyceta  |       sordariomyceta      |                  |                  |
|:--------:|:-----------------:|:-----------------:|:-----------------:|:-------------------------:|:----------------:|:----------------:|
|   Class  |  Sordariomycetes  |  Sordariomycetes  |  Sordariomycetes  |       Leotiomycetes       |  Eurotiomycetes  |  Eurotiomycetes  |
| Subclass | Hypocreomycetidae | Sordariomycetidae | Sordariomycetidae |                           | Eurotiomycetidae | Eurotiomycetidae |
|   Order  |    Hypocreales    |    Sordariales    |   Magnaporthales  |         Helotiales        |    Eurotiales    |    Eurotiales    |
|  Family  |    Hypocreaceae   |    Sordariaceae   |  Magnaporthaceae  |                           |  Aspergillaceae  |  Aspergillaceae  |
|  No rank |                   |                   |                   | Helotiales incertae sedis |                  |                  |
|   Genus  |    Trichoderma    |     Neurospora    |    Magnaporthe    |       Rhynchosporium      |    Aspergillus   |    Penicillium   |
|  Species |                   |                   |                   |                           |                  |                  |

### 02-03-2020 Notes on going through indidivudal _Trichoderma spp._
There are some species that have barely been investigated, some have only been investigated in reference to human health but not for biomass degradation and others have been identified to degrade biomass during investigations to identify fungal species causing rot/degradation to a biomaterial/biomass of interest. These species I am interested in pulling down, especially becuase it appears that for these species the phenotype has been characterised by the specific enzymes involve in the degradation phenotype have not.

However, there is an issue here for some of the species of interst. Some of them do not have assemblies available on NCBI :grimacing: so if I want to go with them I need to see if Leighton has any connections that have non-publically available genomes assemlbies of these species. These species include: 
- _Trichoderma asperelloides_
- _T. pseduokoningii_ 
- _T. deliquescences_ 

Discoverd a new piece of software that can help wit graphically viewing genomic annotations: NCBI's Genome Workbench.

Another helpful website to reference is www.ncbi/nih.gov/assembly/model/

Found a couple of papers that might help explain how to select candidates for the PhD project: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1456863/pdf/tpc1801100.pdf, https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3630013/pdf/2042-5783-3-2.pdf

Need to spend more time looking at the different output file formats of NCBI and what each contains.

Going to spend the rest of the day looking at Entrez.BioPython and if I have time look at the other assemblies available for _Trichoderma spp._ for selection candidates for the NCBI assembly pulldown.

### 03-03-2020
##### Notes on ExtractGenomes_NCBI.py file
Added the getAccessionNumbers function to ExtractGenomes_NCBI.py file. This function extracts the accession numbers of the genomes to be pulled down from the NCBI Assembly database from a plain text file.

Added the storeNCBIdata function to ExtractGenomes_NCBI.py file. This function writes each pulled down Genbank file, pulled down from the NBCI Assembly database, to a new plain text file with the standard naming of 'data-of-pulldown_accession-Number_NCBI-database.txt'.

##### Notes on NCBI Assembly database file formats
Information taken from: https://www.ncbi.nlm.nih.gov/genome/doc/ftpfaq/#files

###### What is the file content within each specific assembly directory?
Assembly directories for all current assemblies, and for many previous assembly versions, include a core set of files and formats plus additional files relevant to the data content of the specific assembly. Directories for old assembly versions that predate the genomes FTP site reorganization contain only the assembly report, assembly stats assembly status files. All data files are named according to the pattern: [assembly accession.version]_[assembly name]_content.[format]

The entries below have the format: filename, download menu name in parentheses, description.

###### assembly_status.txt
A text file reporting the current status of this version of the assembly ("latest", "replaced", or "suppressed"). Any assembly anomalies are also reported.

###### *_assembly_report.txt (Assembly structure report)
Tab-delimited text file reporting the name, role and sequence accession.version for objects in the assembly. The file header contains meta-data for the assembly including: assembly name, assembly accession.version, scientific name of the organism and its taxonomy ID, assembly submitter, and sequence release date.

###### *_assembly_stats.txt (Assembly statistics report)
Tab-delimited text file reporting statistics for the assembly including: total length, ungapped length, contig scaffold counts, contig-N50, scaffold-L50, scaffold-N50, scaffold-N75 scaffold-N90.

###### *_assembly_regions.txt (Assembly regions report)
Provided for assemblies that include alternate or patch assembly units. Tab-delimited text file reporting the location of genomic regions and listing the alt/patch scaffolds placed within those regions.

###### *_assembly_structure directory
Contains AGP files that define how component sequences are organized into scaffolds and/or chromosomes. Other files define how scaffolds and chromosomes are organized into non-nuclear and other assembly-units, and how any alternate or patch scaffolds are placed relative to the chromosomes. Only present if the assembly has internal structure.

###### *_cds_from_genomic.fna.gz (CDS from genomic FASTA)
FASTA format of the nucleotide sequences corresponding to all CDS features annotated on the assembly, based on the genome sequence.

###### *_feature_count.txt.gz (Feature count)
Tab-delimited text file reporting counts of gene, RNA, CDS, and similar features, based on data reported in the *_feature_table.txt.gz file.

###### *_feature_table.txt.gz (Feature table)
Tab-delimited text file reporting locations and attributes for a subset of annotated features. Included feature types are: gene, CDS, RNA (all types), operon, C/V/N/S_region, and V/D/J_segment. Replaces the .ptt .rnt format files that were provided in the old genomes FTP directories.

###### *_genomic.fna.gz (Genomic FASTA)
FASTA format of the genomic sequence(s) in the assembly. Repetitive sequences in eukaryotes are masked to lower-case. The genomic.fna.gz file includes all top-level sequences in the assembly (chromosomes, plasmids, organelles, unlocalized scaffolds, unplaced scaffolds, and any alternate loci or patch scaffolds). Scaffolds that are part of the chromosomes are not included because they are redundant with the chromosome sequences; sequences for these placed scaffolds are provided under the assembly_structure directory.

###### *_genomic.gbff.gz (Genomic GenBank format)
GenBank flat file format of the genomic sequence(s) in the assembly. This file includes both the genomic sequence and the CONTIG description (for CON records), hence, it replaces both the .gbk .gbs format files that were provided in the old genomes FTP directories.

###### *_genomic.gff.gz (Genomic GFF)
Annotation of the genomic sequence(s) in Generic Feature Format Version 3 (GFF3). Sequence identifiers are provided as accession.version. Additional information about NCBI's GFF files is available at ftp://ftp.ncbi.nlm.nih.gov/genomes/README_GFF3.txt.

###### *_genomic.gtf.gz (Genomic GTF)
Annotation of the genomic sequence(s) in Gene Transfer Format Version 2.2 (GTF2.2). Sequence identifiers are provided as accession.version.

###### *_genomic_gaps.txt.gz (Genomic gaps)
Tab-delimited text file reporting the coordinates of all gaps in the top-level genomic sequences. The gaps reported include gaps specified in the AGP files, gaps annotated on the component sequences, and any other run of 10 or more Ns in the sequences.

###### *_protein.faa.gz (Protein FASTA)
FASTA format of the accessioned protein products annotated on the genome assembly.

###### *_protein.gpff.gz (Protein GenPept format)
GenPept format of the accessioned protein products annotated on the genome assembly.

###### *_rm.out.gz (RepeatMasker output)
RepeatMasker output; Provided for eukaryotes.

###### *_rm.run (RepeatMasker run info)
Documentation of the RepeatMasker version, parameters, and library (text format); Provided for eukaryotes.

###### *_rna.fna.gz (RNA FASTA)
FASTA format of accessioned RNA products annotated on the genome assembly; Provided for RefSeq assemblies as relevant (Note, RNA and mRNA products are not instantiated as a separate accessioned record in GenBank and are provided for some RefSeq genomes, most notably the eukaryotes.).

###### *_rna.gbff.gz (RNA GenBank format)
GenBank flat file format of RNA products annotated on the genome assembly; Provided for RefSeq assemblies as relevant.

###### *_rna_from_genomic.fna.gz (RNA from genomic FASTA)
FASTA format of the nucleotide sequences corresponding to all RNA features annotated on the assembly, based on the genome sequence.

###### *_translated_cds.faa.gz (Translated CDS)
FASTA sequences of individual CDS features annotated on the genomic records, conceptually translated into protein sequence. The sequence corresponds to the translation of the nucleotide sequence provided in the *_cds_from_genomic.fna.gz file.

###### *_wgsmaster.gbff.gz (WGS-master)
GenBank flat file format of the WGS master for the assembly (present only if a WGS master record exists for the sequences in the assembly).

###### annotation_hashes.txt
Tab-delimited text file reporting hash values for different aspects of the annotation data. The hashes are useful to monitor for when annotation has changed in a way that is significant for a particular use case and warrants downloading the updated records.

###### md5checksums.txt
File checksums are provided for all data files in the directory.
