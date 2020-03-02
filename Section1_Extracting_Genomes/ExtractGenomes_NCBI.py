# This script is a general script for pulling down genomic assembles from the NCBI Assembly database,
# identified by their taxanomic asscession numbers provided in a separate plain text file.

from Bio import Entrez
Entrez.email = eemh1@st-andrews.ac.uk

def extractAssembly(database, ID, rettype_argmnt, retmode_argmnt)
handel = Entrez.efetch(db = database, id = ID, rettype = rettype_argmnt, retmode = retmode_argmnt)
print(handle.read())
# add code to print out each pulled down genome to a new text file, named after it's ID with a standard name formal as well


# The arguments rettype="gb" and retmode="text" let us download this record in the GenBank format.
