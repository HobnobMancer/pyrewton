def parse_dbcan_output(overview_path, hotpep_file, fasta_file):
    """Parse the data in the 'overview.txt' files from dbCAN.
    
    The 'overview.txt' file contains the CAZyme predicion data from
    HMMER, Hotpep and DIAMOND. A dbCAN consensus result is retrieved
    by finding predictions that at least to of the prediction tools
    agree upon.
    
    Each protein from the input FASTA file containing all query protein
    sequences parsed by dbCAN is represented by a CazymeProteinprediction
    instance. Each predicted CAZyme domains in each predicted CAZyme is
    represented by a CazymeDomain instance, which are sotred within their
    parent CAZyme's CazymeProteinPrediction instance. A CazymeProteinPrediction
    instance is created for each protein in the inpu FASTA file parsed
    by dbCAN and for each CAZyme prediction tool: dbCAN, HMMER, Hotpep and
    DIAMOND.
    
    :param overview_path: Path, path to dbCAN output overview.txt file
    :param hotpep_file: Path, path to Hotpep.out file
    :param fasta_file: Path, path to fasta file containing query protein sequences
          parsed by dbCAN
    
    Return 4 lists, each containing the CazymeProteinPrediction instances for
    HMMER, Hotpep, DIAMOND and dbCAN.
    """
    # open the overview.txt file
    with open(overview_path, "r") as ofh:
        overview_file = ofh.read().splitlines()
    
    # create an empty dictionary to store the CazymeProteinPrediction instances
    # {protein_accesion: CazymeProteinPrediction_instance}
    overview_dict = {}

    # skip the the first line in the file, becuase it contains the 'column' headings
    for line in tqdm(overview_file[1:], desc="Parsing overview.txt"):
        # separate the line content
        line = line.split("\t")
        
        protein_accession = line[0]
        
        # create a CazymeProteinPrediction instance each for HMMER, Hotpep and DIAMOND
        hmmer_prediction = get_hmmer_prediction(line[1], protein_accession)
        hotpep_prediction = get_hotpep_prediction(line[2], protein_accession)
        diamond_prediction = get_diamond_prediction(line[3], protein_accession)
        
        # get the dbCAN consensus result
        dbcan_prediction = get_dbcan_consensus(hmmer_prediction, hotpep_prediction, diamond_prediction, protein_accession, line[-1])
        
        # add the CazymeProteinPrediction instances to the overview_dict
        try:
            overview_dict[protein_accession]
            prediction_outputs = [["HMMER", hmmer_prediction], ["Hotpep", hotpep_prediction], ["DIAMOND", diamond_prediction], ["dbCAN", dbcan_prediction]]
            
            for prediction_tool in prediction_outputs:
                try:
                    overview_dict[protein_accession][prediction_tool[0]]
                    print(
                        f"WARNING -- Protein {protein_accession} already catalogued with {prediction_tool[0]}"
                    )
                except KeyError:  # raised if CazymeProteinPrediciton instance not stored for prediction_tool
                    overview_dict[protein_accession][prediction_tool[0]] = prediction_tool[1]
        
        except KeyError:  # raised if protein has not been added yet
            overview_dict[protein_accession] = {
                "HMMER": hmmer_prediction,
                "Hotpep": hotpep_prediction,
                "DIAMOND": diamond_prediction,
                "dbCAN": dbcan_prediction,
            }
    
    # get EC numbers predicted by Hotpep
    print("proteins in overview_dict pre EC=", len(list(overview_dict.values())))
    overview_dict = get_hotpep_ec_numbers(overview_dict, hotpep_file)
    print("proteins in overview_dictafter EC, pre non=", len(list(overview_dict.values())))
    
    # add proteins predicted to be non-CAZyme to the overview_dict
    overview_dict = get_non_cazymes(fasta_file, overview_dict)
    print("proteins in overview_dict after non=", len(list(overview_dict.values())))

    # now to retrieve a list each of CazymeProteinPredictions for HMMER, Hotpep, DIAMOND and dbCAN
    print("proteins in overview_dict=", len(list(overview_dict.values())))
    print("accessions=\n", list(overview_dict.keys()))
    hmmer_predictions, hotpep_predictions, diamond_predictions, dbcan_predictions = separate_predictions(overview_dict)
    
    return hmmer_predictions, hotpep_predictions, diamond_predictions, dbcan_predictions
        

def get_hmmer_prediction(hmmer_data, protein_accession):
    """Retrieve HMMER prediciton data from the dbCAN overview.txt file.
    
    :param hmmer_data: str, output data from HMMER written in the overview.txt file.
    :param protein_accession: str
    
    Return a CazymeProteinPrediction instance.
    """
    # check if the protein was classified as a non-CAZyme
    if hmmer_data == "-":
        cazyme_classification = 0  # non-CAZyme
        protein = CazymeProteinPrediction("HMMER", protein_accession, cazyme_classification)
        return protein
    
    # create a CazymeProteinPrediction instance to represent the CAZyme
    cazyme_classification = 1
    parent_cazyme = CazymeProteinPrediction("HMMER", protein_accession, cazyme_classification)

    # HMMER can predicted multiple CAZyme domaisn for a CAZyme
    # HMMER separates each of the predicted CAZyme domains by additiona symbols
    parent_cazyme.cazyme_domains = get_hmmer_cazyme_domains(hmmer_data.split("+"), protein_accession)
    
    return parent_cazyme

    
def get_hmmer_cazyme_domains(predicted_domains, protein_accession):
    """Retrieve the CAZyme domains predicted by HMMER.
    
    :param predicted_domains: list of strings, each string containing a CAZyme domain prediction
    :param protein_accession: str
    
    Return list of CazymeDomain instances. Each instance represents a unique CAZyme domain.
    """
    cazyme_domains = {}  # {"cazyfam-cazysubfam": CazymeDomain_instance}

    for domain in predicted_domains:
        # separate the predicted CAZy family from the domain range
        domain = domain.split("(")
        
        # get the predicted CAZy family(subfamily)
        domain_name = domain[0]  # CAZy (sub)family
        cazy_family, cazy_subfamily = get_hmmer_cazy_family(domain_name)
        
        if cazy_family is None:  # unexpected data format found by get_hmmer_cazy_family()
            continue
        
        # get the predicted amino acid domain range
        domain_range = domain[1][:-1]  # drop the terminal ")"
        # check it matches the expected format and thus is the predicted domain range
        try:
            re.match(r"\d+?-\d+", domain_range).group()
        except AttributeError:  # raised if not correct format
            print(
                f"WARNING -- UNKNOWN DATA TYPE OF {domain_range} FOR {protein_accession}\n"
                f"Not adding domain range for the domain {cazy_family}-{cazy_subfamily}"
            )
            domain_range = np.nan

        # create identifier for the CAZyme domain to check if it has already been catalogued
        domain_identifier = f"{cazy_family}-{cazy_subfamily}"
    
        # check if the CAZyme domain has already been cataglouged
        try:
            existing_domain = cazyme_domains[domain_identifier]

            # CAZyme domain already catalogued, check if the retrieved domain range is the same
            if domain_range not in existing_domain.domain_range:
                # append the new domain_range
                existing_domain.domain_range.append(domain_range)
                print("multiple domain ranges=", existing_domain, existing_domain.domain_range, protein_accession)

        except KeyError:  # raised if the CAZyme domain has not yet been catalogued
            cazyme_domains[domain_identifier] = CazymeDomain(
                prediction_tool="HMMER",
                protein_accession=protein_accession,
                cazy_family=cazy_family,
                cazy_subfamily=cazy_subfamily,
                domain_range=[domain_range],
            )

    return list(cazyme_domains.values())


def get_hmmer_cazy_family(domain_name):
    """Retrieve the name of the CAZyme domain, this is the predicted CAZy family and subfamily.
    
    If a subfamily is not predicted, the CAZy subfamily is set as a null value.
    If the data formate does not match the expected format, the CAZy family will
    be set by None, which will be used to not add the current working domain to
    the parent CazymeProteinPrediction instance.
    
    :param domain_name: str, contains the CAZy (sub)family
    
    Return two strings, the CAZy family and the CAZy subfamily.
    """
    # check if unusal CAZy family format or subfamily was given
    if domain_name.find("_") != -1:
        try:
            re.match(r"\D{2,3}\d+?_\D", domain_name).group()  # check unusal CAZy family formating
            cazy_family = domain_name[:domain_name.find("_")]
            cazy_subfamily = np.nan

        except AttributeError:  # raised if not an usual CAZy family format
            try:
                # check if its a CAZy subfamily
                re.match(r"\D{2,3}\d+?_\d+", domain_name).group()
                cazy_family = domain_name[:domain_name.find("_")]
                cazy_subfamily = np.nan

            except AttributeError:  # raised if it doesn't match the expected CAZy subfamily format
                print(
                    f"WARNING -- UNKNOWN DATA TYPE OF {domain_name} FOR {protein_accession}\n"
                    "Not adding as a domain."
                )
                return None, None

    else:
        # potentially have a CAZy family
        try:
            re.match(r"\D{2,3}\d+", domain_name).group()
            cazy_family = domain_name
            cazy_subfamily = np.nan

        except AttributeError:  # raised if doesn't match expected CAZy family
            print(
                f"WARNING -- UNKNOWN DATA TYPE OF {domain_name} FOR {protein_accession}\n"
                "Not adding as a domain."
            )
            return None, None
    
    return cazy_family, cazy_subfamily


def get_hotpep_prediction(hotpep_data, protein_accession):
    """Retrieve Hotpep prediciton data from the dbCAN overview.txt file.
    
    :param hmmer_data: str, output data from Hotpep written in the overview.txt file.
    :param protein_accession: str
    
    Return a CazymeProteinPrediction instance.
    """
    # check if the protein was classified as a non-CAZyme
    if hotpep_data == "-":
        cazyme_classification = 0  # non-CAZyme
        protein = CazymeProteinPrediction("Hotpep", protein_accession, cazyme_classification)
        return protein
    
    # create a CazymeProteinPrediction instance to represent the CAZyme
    cazyme_classification = 1
    parent_cazyme = CazymeProteinPrediction("Hotpep", protein_accession, cazyme_classification)

    # Hotpep can predicted multiple CAZyme domaisn for a CAZyme
    # Hotpep separates each of the predicted CAZyme domains by additiona symbols
    parent_cazyme.cazyme_domains = get_hotpep_cazyme_domains(hotpep_data.split("+"), protein_accession)
    
    return parent_cazyme


def get_hotpep_cazyme_domains(predicted_domains, protein_accession):
    """Retrieve the CAZyme domains predicted by Hotpep.
    
    :param predicted_domains: list of strings, each string containing a CAZyme domain prediction.
    
    Return list of CazymeDomain instances. Each instance represents a unique CAZyme domain.
    """
    cazyme_domains = {}  # {"cazyfam-cazysubfam": CazymeDomain}
    
    for domain in predicted_domains:
        # remove k-mer group number
        domain = domain.split("(")[0]
        
        # check if a CAZy subfamily was predicted
        if domain.find("_") != -1:  # found a CAZy subfamily
            try:
                # check if its a CAZy subfamily
                re.match(r"\D{2,3}\d+?_\d+", domain).group()
                cazy_family = domain[:domain.find("_")]
                cazy_subfamily = np.nan

            except AttributeError:  # raised if it doesn't match the expected CAZy subfamily format
                print(
                    f"WARNING -- UNKNOWN DATA TYPE OF {domain} FOR {protein_accession}\n"
                    "Not adding as a domain."
                )
                continue
        
        else:  # found a CAZy family
            try:  # check it is a CAZy family
                re.match(r"\D{2,3}\d+", domain).group()
                cazy_family = domain
                cazy_subfamily = np.nan

            except AttributeError:  # raised if doesn't match expected CAZy family
                print(
                    f"WARNING -- UNKNOWN DATA TYPE OF {domain} FOR {protein_accession}\n"
                    "Not adding as a domain."
                )
                continue
        
        new_domain = CazymeDomain("Hotpep", protein_accession, cazy_family, cazy_subfamily)
        
        try:
            cazyme_domains[f"{cazy_family}-{cazy_subfamily}"]
            print(
                f"WARNING -- DUPLICATE CAZYME DOMAIN PREDICTIONS OF CAZy family = {cazy_family},\n"
                f" CAZy subfamily = {cazy_subfamily}, in protein {protein_accession}.\n"
                "Only one CazymeDomain instance being added to the respective CazymeProteinPrediction instance."
            )
        except KeyError:
            cazyme_domains[f"{cazy_family}-{cazy_subfamily}"] = new_domain
    
    return list(cazyme_domains.values())
        

def get_hotpep_ec_numbers(overview_dict, hotpep_file):
    """Retrieve EC numbers from the Hotpep.out file, and add them to the CAZyme domains.
    
    :param overview_dict: dict, keyed by protein accession, valued by dictionary of 
                            predicion_tool:CazymeDomain
    :param hotpep_file: Path, path to Hotpep.out file
    
    Return overview_dict
    """
    with open (hotpep_file, "r") as fh:
        hotpep_file_data = fh.read().splitlines()
    
    for line in tqdm(hotpep_file_data[1:], desc="Parsing Hotpep.out"):  # skip the first line containing headings
        # separate out the data
        line = line.split("\t")
        
        protein_accession = line[2].strip()
        
        # retrieve the EC numbers
        ec_numbers = line[-1]
        ec_numbers = ec_numbers.split(",")

        # strip any remaining white space becuase the separating comma
        # is sometimes followed by a space, and remove the kmer cluster number
        index = 0
        for index in range(len(ec_numbers)):
            ec_numbers[index] = ec_numbers[index].strip()
            # remove k-mer group number (EC:kmer_group#)
            ec_numbers[index] = ec_numbers[index].split(":")[0]
        
            if ec_numbers[index] == "NA":
                ec_numbers[index] = np.nan
                continue

            try:
                re.match(r"\d+?\.(\d+?|-)\.(\d+?|-)\.(\d+?|-)", ec_numbers[index])  # missing digits listed as '-'

            except AttributeError:  # raised if not in the expected EC number format
                print(
                    f"WARNING -- unexpected data format of '{ec_numbers[index]}' for protein "
                    f"{protein_accession}.\n"
                    "Not adding this data to the CAZyme's CazymeProteinPrediction instance."
                )
                ec_numbers.remove(ec_numbers[index])
        
        # get CAZy family and subfamily of the current working domain
        cazy_family = line[0]
        if cazy_family.find("_") != -1:
            cazy_subfamily = cazy_family
            cazy_family = cazy_family[:cazy_family.find("_")]
        else:
            cazy_subfamily = np.nan

        domain_identifer = [cazy_family, cazy_subfamily]
        
        parent_cazyme = overview_dict[protein_accession]["Hotpep"]
        parent_domains = parent_cazyme.cazyme_domains
        
        for domain in parent_domains:
            if domain_identifer == [domain.cazy_family, domain.cazy_subfamily]:
                domain.ec_numbers = ec_numbers
        
        parent_cazyme.cazyme_domains = parent_domains
        overview_dict[protein_accession]["Hotpep"] = parent_cazyme
        
    return overview_dict


def get_diamond_prediction(diamond_data, protein_accession):
    """Retrieve DIAMOND prediction data from dbCAN overview.txt file.
    
    :param diamond_data, str, output data from DIAMOND written in the overview.txt
    :param protein_accession, str
    
    Return a CazymeProteinPrediction instance.
    """
    # check if the protein was classified as a non-CAZyme
    if diamond_data == "-":
        cazyme_classification = 0  # non-CAZyme
        protein = CazymeProteinPrediction("DIAMOND", protein_accession, cazyme_classification)
        return protein
    
    # create a CazymeProteinPrediction instance to represent the CAZyme
    cazyme_classification = 1
    parent_cazyme = CazymeProteinPrediction("DIAMOND", protein_accession, cazyme_classification)
    
    # DIAMOND can predict multiple CAZyme domains per CAZyme, each domain is separated by '+'
    parent_cazyme.cazyme_domains = get_diamond_predicted_domains(diamond_data.split("+"), protein_accession)
 
    return parent_cazyme


def get_diamond_predicted_domains(predicted_domains, protein_accession):
    """Retrieve the CAZyme domains predicted by DIAMOND, from the overview.txt file.
    
    :param predicted_domains: list of strings, each string contains a unique CAZyme domain prediction
    :param protein_accession: str
    
    Return list of CazymeDomain instances, one instance per unique CAZyme domain.
    """
    cazyme_domains = {}  # store CazymeDomain instances, keyed by CAZyme domain "fam-subfam"
    
    for predicted_domain in predicted_domains:
        # check if a subfamily was predicted
        try:
            re.match(r"\D{2,3}\d+?_\d+", predicted_domain).group()
            cazy_family = predicted_domain[:predicted_domain.find("_")]
            cazy_subfamily = predicted_domain
        
        except AttributeError:  # raised if predicted_domain is not a CAZy subfamily
            try:
                # check if a CAZy family was predicted
                re.match(r"\D{2,3}\d+", predicted_domain).group()
                cazy_family = predicted_domain
                cazy_subfamily = np.nan
            
            except AttributeError:  # raised if not a CAZy family
                print(
                    f"WARNING -- Unexpected data formate for '{predicted_domain}' for\n"
                    "DIAMOND in the overview.txt file. Not adding domain to the CAZyme."
                )
            continue
        
        new_domain = CazymeDomain("DIAMOND", protein_accession, cazy_family, cazy_subfamily)
        
        # check if the domain has been parsed before
        try:
            cazyme_domains[f"{cazy_family}-{cazy_subfamily}"]
            print(
                f"WARNING -- duplicate CAZyme domains predicted for protein {protein_accession} "
                " by DIAMOND.\n"
                f"CAZyme domain (CAZy family = {cazy_family}, CAZy subfamily = {cazy_subfamily}),\n"
                "added to the protein once as a single CazymeDomain instance"
            )
            continue
        
        except KeyError:
            cazyme_domains[f"{cazy_family}-{cazy_subfamily}"] = new_domain

    
    return list(cazyme_domains.values())


def get_dbcan_consensus(hmmer_prediction, hotpep_prediction, diamond_prediction, protein_prediction, no_of_tools):
    """Retrieve the consensus CAZyme and CAZyme domain prediction for dbCAN.
    
    Consensus is defined as a result that at least two of the tools (HMMER, Hotpep and DIAMOND)
    agree upon.
    
    :param hmmer_prediction: CazymeProteinPrediction instance containing the HMMER prediction
    :param hotpep_prediction: CazymeProteinPrediction instance containing the Hotpep prediction
    :param diamond_prediction: CazymeProteinPrediction instance containing the DIAMOND prediction
    :param protein_accession: str
    :param no_of_tools: str, data from the '#ofTools' column from the overview.txt file, it states
        the number of tools that list the protein as a CAZyme
    
    Return CazymeProteinPrediction instance to represent the dbCAN consensus prediction
    """
    # first check if the consensus was that the protein was a non-CAZyme
    if no_of_tools == "1":
        # non-CAZyme
        cazyme_classification = 0  # non-CAZyme
        protein = CazymeProteinPrediction("dbCAN", protein_accession, cazyme_classification)
        return protein
    
    # get the consensus CAZyme domain predicitons
    
    # first get each of the predicted CAZyme domains for each tool
    hmmer_domains = get_dbcan_cazyme_domains(hmmer_prediction)
    hotpep_domains = get_dbcan_cazyme_domains(hotpep_prediction)
    diamond_domains = get_dbcan_cazyme_domains(diamond_prediction)

    # get the predicted CAZyme domains that are predicted by at two tools
    hmmer_hotpep_consensus = list(set(hmmer_domains) & set(hotpep_domains))
    hotpep_diamond_consensus = list(set(hotpep_domains) & set(diamond_domains))
    hmmer_diamond_consensus = list(set(hmmer_domains) & set(diamond_domains))

    # get the CAZyme domains predicted by all prediction tools
    all_consensus = list(set(hmmer_domains) & set(hotpep_domains) & set(diamond_domains))

    # combine all consensus results
    combined_consensus = (
        hmmer_hotpep_consensus + hotpep_diamond_consensus + hmmer_diamond_consensus + all_consensus
    )
    # remove duplicate entries
    combined_consensus = list(set(combined_consensus))
    
    # create a CazymeProteinPrediction instance to represent the CAZyme
    cazyme_classification = 1  # CAZyme
    parent_cazyme = CazymeProteinPrediction("dbCAN", protein_accession, cazyme_classification)
    
    # create a CazymeDomain instance for each domain
    cazyme_domains = []  # store CazymeDomain instances
    for domain in combined_consensus:
        domain = domain.split("**")
        cazy_family = domain[0]

        if domain[1] == "nan":
            cazy_subfamily = np.nan
        else:
            cazy_subfamily = domain[1]

        new_domain = CazymeDomain("dbCAN", protein_accession, cazy_family, cazy_subfamily)

        cazyme_domains.append(new_domain)
    
    # add predicted CAZyme domains to the instance representing the CAZyme
    parent_cazyme.cazyme_domains = cazyme_domains
    
    return parent_cazyme
    
    
def get_dbcan_cazyme_domains(cazyme_prediction):
    """Retrieve the predicted CAZyme domains, representing each domains as a string 'fam_subfam'.
    
    :param cazyme_prediction: CazymeProteinPrediction instance.
    
    Return a list of strings with each string representing a single predicted CAZyme domain.
    """
    cazyme_domain_instances = cazyme_prediction.cazyme_domains
    
    cazyme_domains = []  # store list representing CAZyme domains [fam(str), subfam(str)]
    
    # if a prediction tool predicts a protein is a non_CAZyme, non_cazymes == np.nan (a float)
    if type(cazyme_domain_instances) is not float:
        for domain in cazyme_domain_instances:
            cazy_family = domain.cazy_family
            cazy_subfamily = domain.cazy_subfamily
            # separate the family and subfamily '**' becuase these do not appear in the either
            cazyme_domains.append(f"{cazy_family}**{cazy_subfamily}")
    
    return cazyme_domains


def get_non_cazymes(fasta_path, overview_dict):
    """Add proteins that dbCAN identified as non-CAZymes to the collection of CazymeProteinPrediction instances.
    
    :param fasta_path: Path, FASTA file used as input for CUPP.
    :param overview_dict: dict, key=protein_accession, value=dict valued by prediction tool,
          valued by CazymeProteinPrediction intance
    
    Return a dictionary keyed by protein accessions and valued by dictionary of tool:CazymeProteinPrediction.
    """
    # open the FASTA path
    with open(fasta_path) as fh:
        fasta = fh.read().splitlines()

    for line in tqdm(fasta, desc="Adding non-CAZymes"):
        if line.startswith(">"):
            protein_accession = line[1:].strip()
            
            # check if the protein has been listed as a CAZyme by CUPP
            try:
                overview_dict[protein_accession]
            except KeyError:
                # raised if protein not in cupp_predictions, inferring proten was not labelled as CAZyme
                cazyme_classification = 0
                overview_dict[protein_accession] = {}
                for prediction_tool in ["HMMER", "Hotpep", "DIAMOND", "dbCAN"]:
                    overview_dict[protein_accession][prediction_tool] = CazymeProteinPrediction(
                        prediction_tool,
                        protein_accession,
                        cazyme_classification,
                    )

    return overview_dict


def separate_predictions(overview_dict):
    """Create separate lists of CazymeProteinPrediction intances for each prediction tool.
    
    :param overview_dict: dict containing all CazymeProteinPrediciton instances
        {protein_accession: {tool: CazymeProteinPrediction instance}}
    
    Return 4 lists of CazymeProteinPrediction instances    
    """
    hmmer_predictions = []
    hotpep_predictions = []
    diamond_predictions = []
    dbcan_predictions = []
    for protein_accession in list(overview_dict.keys()):
        try:
            hmmer_predictions.append(overview_dict[protein_accession]["HMMER"])
        except KeyError:
            print(
                f"Could not find HMMER prediction for {protein_accession}"
            )

        try:
            hotpep_predictions.append(overview_dict[protein_accession]["Hotpep"])
        except KeyError:
            print(
                f"Could not find HMMER prediction for {protein_accession}"
            )

        try:
            diamond_predictions.append(overview_dict[protein_accession]["DIAMOND"])
        except KeyError:
            print(
                f"Could not find HMMER prediction for {protein_accession}"
            )

        try:
            dbcan_predictions.append(overview_dict[protein_accession]["dbCAN"])
        except KeyError:
            print(
                f"Could not find HMMER prediction for {protein_accession}"
            )
    
    return hmmer_predictions, hotpep_predictions, diamond_predictions, dbcan_predictions
