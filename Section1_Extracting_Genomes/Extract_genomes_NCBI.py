#!/usr/bin/env python
# -*- coding: utf-8 -*-

# This script is a general script for pulling down genomic assembles from the NCBI Assembly database,
# identified by their taxanomic asscession numbers provided in a separate plain text file.

from Bio import Entrez

# Fill out pull down form:
Entrez.email = eemh1@st-andrews.ac.uk # enter email of person performing the pulldown
date_of_pulldown = '03-03-2020' # enter the start date of executing this script, as a string
NCBI_database = 'NCBIAssembly' # name the database from which data will be pulled from, as a string




# Create function to continually access desired NCBI database, and pull down data in desired output format.
def extract_assembly(accession_number):
    """This function pulls down the data from the specified NCBI database in the specified format.
    For information about Entrez.BioPython module's keywords for rettype and retmode (which define
    the output format of the pulled down data, see:
    https://biopython.org/DIST/docs/api/Bio.Entrez-module.html
    https://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc111 >> Chapter 9
    This function uses the Entrez eFetch handel which retrieves data from NCBI, for more information see:
    https://www.ncbi.nlm.nih.gov/books/NBK25499/#chapter4.EFetch"""

    # The arguments rettype="gb" and retmode="text" will pull down the record in the GenBank format.
    rettype_argmnt = "gb"
    retmode_argmnt = "text"

    handel = Entrez.efetch(db=NCBI_database, id=accession_number, rettype=rettype_argmnt, retmode=retmode_argmnt)
    print(handel.read())
    


# Define function to extract accession numbers of datasets of interest in NCBI database.
def get_accession_numbers(accession_numbers_file):
    """Extracts the accession numbers from the 'inputFILE' containing the accession numbers.
    The input file must be a plain text file in the same folder as this script.
    Additionally, the input file must not contain any blank/empty lines and only contain the accession numbers.
    The accession numbers will be returned as a list with each number being assigned an individual element"""

    # inputFILE = open(r" - path for file containing accession numbers - ")
    input_file = open(r"/home/em/PhD_Project_Scripts/Section_1_Extracting_Genomes/NCBI_AccessionNumbers_inputFILE.txt", "rt")
    accesssion_numbers = input_file.readlines()
    input_file.close()
    return(accesssion_numbers)


# Define function to store pulled down data from NCBI in a individual data file for that associated accession number
def store_ncbi_data(working_accession_number, NCBI_data):
    """This function stores the data pulled down from the NCBI database in a new file, named in the predefined manner.
    File naming system: 'data-of-pulldown_Accession#_NCBI-database-name.
    Example file name: 02-03-2020_GCA00014548_NCBIAssembly
    The file extension will also be defined using this function.
    The created output file will be stored in the same directory as this scritp"""
    # date-of-pulldown  and NCBIdatabase (the database from which data will be pulled) are defined at the top the script
    data_accession_number = NCBI_data # the accession number of the NCBI data of interest is passed to the function

    file_format = '.txt' # define the file extention for the output file
    file_name = date_of_pulldown + '_' + working_accession_number + '_' + NCBI_database + file_format

    output_file = open(file_name, "w")    # create file
    output_file.write(NCBI_data)          # write out pulled down data
    output_file.close()                  # close file



# add code to print out each pulled down genome to a new text file, named after it's ID with a standard name
# formal as well

# Extract
