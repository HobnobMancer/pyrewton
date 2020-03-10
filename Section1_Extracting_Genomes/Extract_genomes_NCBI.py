#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Script pulls down assembles from NCBI database

from Bio import Entrez
import datetime

# Fill out pull down form:
Entrez.email = 'eemh1@st-andrews.ac.uk'  # enter email address
date = datetime.datetime.now()
date_of_pulldown = date.strftime("%Y-%m-%d")
NCBI_database = 'NCBIAssembly'  # name the database as a string


# Define function to extract accession numbers from plain text file
def get_accession_numbers(accession_numbers_file):
    """Extract accession numbers from a plain text file.

    This function extract the accession numbers of
    genomic sequences stored in the NCBI assembly database.
    The input file must be a plain text file in
    the same folder as this script.
    Additionally, the input file must not contain
    any blank/empty lines and only contain the accession numbers.
    The accession numbers will be returned as a
    list with each number being assigned an individual element
    """

    # inputFILE = open(r" - path for file containing accession numbers - ")
    accession_numbers = []
    with open("input_file.txt", "r") as file_handle:
        accession_numbers = file_handle.readlines
    return(accession_numbers)


# Create function download as save assembly data from NCBI
def extract_assembly(accession_number):
    """Pull down genomic sequences from the NCBI assembly database.

    This function pulls down the data from the specified NCBI
    database in the specified format.
    For information about Entrez.BioPython module's
    keywords for rettype and retmode (which define
    the output format of the pulled down data, see:
    https://biopython.org/DIST/docs/api/Bio.Entrez-module.html
    https://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc111
    >> Chapter 9
    This function uses the Entrez eFetch handel which retrieves
    data from NCBI, for more information see:
    https://www.ncbi.nlm.nih.gov/books/NBK25499/#chapter4.EFetch
    
    
    The downloaded assembly is written to an output file,
    with the standard file name of 'date - accession number
    - NCBI database'.
    """

    # The arguments rettype="gb" and retmode="text" pull down records in the GenBank format
    rettype_argmnt = "gb"
    retmode_argmnt = "text"

    handle = Entrez.efetch(db=NCBI_database, id=accession_number, rettype=rettype_argmnt, retmode=retmode_argmnt)
    assembly = handle.read()

    # Write data to ouput file
    file_format = '.txt' # define the file extention for the output file
    file_name = date_of_pulldown + '_' + accession_number + '_' + NCBI_database + file_format

    with open(file_name, "w") as output_handle:
        output_handle.write(assembly)
    
    handle.close()



# Main script body

# Retrieve accession numbers from input pile
accession_numbers = get_accession_numbers('NCBI_accession_numbers.txt')

# Pull down genomic sequences from NCBI and store in individual files
for record in accession_numbers:
    extract_assembly(record)
