# This script is a general script for pulling down genomic assembles from the NCBI Assembly database,
# identified by their taxanomic asscession numbers provided in a separate plain text file.

from Bio import Entrez
Entrez.email = eemh1@st-andrews.ac.uk


# Create function to continually access desired NCBI database, and pull down data in desired output format.
def extractAssembly(database, ID, rettype_argmnt, retmode_argmnt):
    """This function pulls down the data from the specied NCBI database in the specified format.
    For information about Entrez.BioPython module's keywords for rettype and retmode (which define
    the output format of the pulled down data, see:
    https://biopython.org/DIST/docs/api/Bio.Entrez-module.html
    https://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc111 >> Chapter 9
    This function uses the Entrez eFetch handel which retrieves data from NCBI, for more information see:
    https://www.ncbi.nlm.nih.gov/books/NBK25499/#chapter4.EFetch
    """
    
    handel = Entrez.efetch(db=database, id=ID, rettype=rettype_argmnt, retmode=retmode_argmnt)
    print(handel.read())
    
# Define function to extract accession numbers of datasets of interest in NCBI database.
# The accession numbers must be stored in a plain text file in the same folder as this script.
# Additionally, the plain text file must contain no blank/empty lines and only the accession numbers.
def getAccessionNumbers(accessionNumbersFile):
    
    
# add code to print out each pulled down genome to a new text file, named after it's ID with a standard name
# formal as well

# Extract 





# The arguments rettype="gb" and retmode="text" let us download this record in the GenBank format.
