# Parsing Output from CAZyme Prediction Tools

Each of the CAZyme prediction tools, dbCAN, CUPP and eCAMI, store their output in different output file types and store the data in different formats. The submodule `pyrewton.cazymes.prediction.parse` contains Python scripts that parse and standardise the output from each of the CAZyme prediction tools.

The retrieved prediction data is stored in two data models:

**Proteins** are represented by the `CazymeProteinPrediction` class, which contains the `cazyme_classification` binary attribute, which represents if the protein was predicted to be a CAZyme (and awarded a score of 1) or a non-CAZyme (awarded a score of 0).

Each `CazymeProteinPrediction` possesses a list of `cazyme_domains`, with length 0 for non-CAZymes and length greater than or equal to 1 for CAZymes. Each element in `cazyme_domains` represents a single, unique predicted CAZyme domain, which is represented by a `CazymeDomain` instance. Each predicted CAZyme domain is identified by a unique CAZy family - CAZy subfamily combination.

CAZy family and CAZy subfamily predicted annotations are retrieved per predicted CAZyme domain from all CAZyme prediction tools.

In addition to the CAZy family and subfamily annotations, the following data is retrieved from each CAZyme prediction tool:

**dbCAN**  
_Nothing extra_  

**HMMER**  
- Domain range (in amino acids)  

**Hotpep**  
- EC numbers  

**DIAMOND**  
_Nothing extra_  

**CUPP**
- Domain range (in amino acids)
- EC numbers  

**eCAMI**
- EC numbers  
