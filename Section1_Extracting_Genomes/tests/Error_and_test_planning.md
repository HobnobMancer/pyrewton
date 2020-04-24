# Predicated errors and bugs in the Python Extract_genomes_NCBI.py file

## Predicated bugs/errors

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
Expected error: index error
Error message recieved: `IndexError: list index out of range`
Program terminates when trying to retrieve name from retrieved entrez record
Potential cause of error: typo in species name leading to no entry from entrez containing\
the information of a species being pulled down
In error log, include which line was being processed to make it easier for use to double\
check the probable cause.

3c. Calling get_tax_id()
Expected error: index error
Error message recieved: `IndexError: list index out of range`
Program terminates when trying to retrieve tax id from retrieved entrez record
Potential cause of error: typo in species name or taxonomy id (so that taxonomy id it fails\
taxonomy id check and script processes it as a species name), leading to no entry from entrez\
containing the information of a species being pulled down
In error log, include which line was being processed to make it easier for use to double\
check the probable cause.

3d. Creating dataframe
Potential error of empty entries if nothing was retrieved from the entrez call.
_Look up if pandas can be used to check for this_

4. Errors in collate_accession_numbers()
