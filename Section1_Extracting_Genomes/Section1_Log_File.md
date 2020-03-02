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
