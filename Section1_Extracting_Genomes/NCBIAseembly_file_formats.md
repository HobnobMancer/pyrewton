## Notes on NCBI Assembly database file formats
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
