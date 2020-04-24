# Predicated errors and bugs in the Python Extract_genomes_NCBI.py file

## Predicated bugs/errors

_Solution to all_ `IOerror: network error` _is to retry._
_Do this automatically by writing in the script to retry or_\
_make the user manually restart the script? The later seems_\
_like a better idea. Entrez already does a try of three,_
_does this count for network errors?_
_Could a retry system in the script be written so the retries_\
_or spaced out by time, with default values in the script,_\
_but the user can change the max number of tries and time_\
_between the tries if they want?_

1. Calling parse_input_file()
Expected error type: ?
Potential errors more likely to occur as a result of one of the functions\
being called within parse_line_file() encountering an area. Therefore,\
the error type will be very broad.

2. Calling collate_accession_numbers()
Expected error type: ?
Potential errors more likely to occur as a result of one of the functions\
being called within collate_accession_numbers() encountering an area.\
Therefore, the error type will be very broad.

3. Errors in parse_input_file()
3a. Calling the input file.
Expected error: file not found
Error message recieved: `FileNotFoundError: [Errno 2] No such file or directory: '<filename>'`
Potential cause of error: file does not exist in cwd, or not in directory given by input\
or a typo in the input of filename and/or path.

3b. Calling get_genus_species_name()

i. 
Expected error: index error
Error message recieved: `IndexError: list index out of range`
Program terminates when trying to retrieve name from retrieved entrez record
Potential cause of error: typo in species name leading to no entry from entrez containing\
the information of a species being pulled down
In error log, include which line was being processed to make it easier for use to double\
check the probable cause.

ii.
Expected error: IOError
Cause of error: Network error when using Entrez to call to NCBI

3c. Calling get_tax_id()

i.
Expected error: index error
Error message recieved: `IndexError: list index out of range`
Program terminates when trying to retrieve tax id from retrieved entrez record
Potential cause of error: typo in species name or taxonomy id (so that taxonomy id it fails\
taxonomy id check and script processes it as a species name), leading to no entry from entrez\
containing the information of a species being pulled down
In error log, include which line was being processed to make it easier for use to double\
check the probable cause.

ii.
Expected error: `IOError`
Cause of error: Network error when using Entrez to call to NCBI

3d. Creating dataframe
Potential error of empty entries if nothing was retrieved from the entrez call.
_Look up if pandas can be used to check for this_

4. Errors in collate_accession_numbers()
4a. Calling get_accession_numbers()
These areas will mainly arise from the Entrez calls

4ai. Using Entrez to retrieve assembly ids
i.
Error message recieved: `IndexError: list index out of range`
Cause of error: If taxonomy ID is not recognised, terminates programme
Edit so prints out taxonomy ID with the issue

ii.
Expected error: `IOError`
Cause of error: Network error

4aii. Using Entrez to post assembly ids
i
Expected error: `RuntimeError(value)`
Error message: `RuntimeError: cannot get document summary`
Cause of error: incorretly formated assembly ids into query, or too large post\
try reducing total number of species
_Should I add an option for reduced speed, add it as an args question,_\
_so that if yes, if statement checks it here and if yes adds a 1 sec_\
_sleep timer, to prevent over demand of network that causes you to be_\
_kicked off._

ii
Expected error: `IOError`
Cause of error: Network error

4aiii. Using Entrez to retrieve accession numbers per assembly id
i.
Expected error: `RuntimeError(value)`
Error message: `RuntimeError: cannot get document summary`
Cause of error: incorretly formated assembly ids into query, or too large post\
try reducing total number of species, or too mant requests per second
_Should I add an option for reduced speed, add it as an args question,_\
_so that if yes, if statement checks it here and if yes adds a 1 sec_\
_sleep timer, to prevent over demand of network that causes you to be_\
_kicked off._

ii.
Expected error: `IOError`
Cause of error: Network error when using Entrez to call to NCBI
