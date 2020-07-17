# README for supplementary data

This directory contains supplementary data for the thesis accompanying the EASTBIO PhD Project.

## Contents

1. [CAZy EC numbers](#CAZy)
2. [Cazyme GO terms](#Cazyme)

## CAZy EC numbers

Contains EC numbers retrieved from the [CAZy database](http://www.cazy.org/#:~:text=The%20CAZy%20database%20describes%20the%20families%20of%20structurally-related,enzymes%20that%20degrade%2C%20modify%2C%20or%20create%20glycosidic%20bonds.)

Incomplete EC numbers retrieved from CAZy are listeted at the top of each column. The rows underneath contain the identifer (fourth) digit of the EC number when retrieved from CAZy as a complete EC number.

EC numbers were retrieved on 2020-07-09, the most recent CAZy update had been on the 2020-06-26.

## Cazyme GO terms

Constitutes too files. The file ending in 'all' contains all GO function terms retrieved from [QuickGO](https://www.ebi.ac.uk/QuickGO/), with duplicates excluded, therefore, each row contains a unique GO term. The file ending 'used_in_uniprotkb' contains all the GO function terms retrieved from QuickGO that are used at least once in UniProtKB to annotate a gene product, i.e. the file contains all entries from 'potential_cazyme_inferring_GO_terms_all.csv' that have more than 0 annotations count.

GO terms were retrieved on 2020-07-12 and 2020-07-13, the most recent [GO (Gene Ontology)](http://geneontology.org/#:~:text=Gene%20Ontology%20Resource%20%7C%20The%20Gene%20Ontology%20%28GO%29,at%20the%20molecular%2C%20cellular%20and%20tissue%20system%20levels.) update had been on the 2020-06-01.
