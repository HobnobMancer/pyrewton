.. pyrewton documentation master file, created by
   sphinx-quickstart on Fri Jul 17 10:28:58 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

=====================================
Welcome to pyrewton's documentation!
=====================================

| 
| Version v0.1.3 2021/11/05
| DOI: 10.5281/zendo.3876218 
| `GitHub repository <https://github.com/HobnobMancer/PhD_Project_Scripts/tree/master>`_
| 

|DOI| |Funding| |PhD licence| |CircleCI| |codecov| |Documentation
Status| |Python| |Research|

.. |DOI| image:: https://zenodo.org/badge/DOI/10.5281/zenodo.3876218.svg
   :target: https://doi.org/10.5281/zenodo.3876218
.. |Funding| image:: https://img.shields.io/badge/Funding-EASTBio-blue
   :target: http://www.eastscotbiodtp.ac.uk/
.. |PhD licence| image:: https://img.shields.io/badge/Licence-MIT-green
   :target: https://github.com/HobnobMancer/PhD_Project_Scripts/blob/master/LICENSE
.. |CircleCI| image:: https://circleci.com/gh/HobnobMancer/pyrewton.svg?style=shield
   :target: https://circleci.com/gh/HobnobMancer/pyrewton
.. |codecov| image:: https://codecov.io/gh/HobnobMancer/pyrewton/branch/master/graph/badge.svg
   :target: https://codecov.io/gh/HobnobMancer/pyrewton
.. |Documentation Status| image:: https://readthedocs.org/projects/pyrewton/badge/?version=latest
   :target: https://pyrewton.readthedocs.io/en/latest/?badge=latest
.. |Python| image:: https://img.shields.io/badge/Python-v3.7.---orange
   :target: https://www.python.org/about/
.. |Research| image:: https://img.shields.io/badge/Bioinformatics-Protein%20Engineering-ff69b4
   :target: http://www.eastscotbiodtp.ac.uk/eastbio-student-cohort-2019

-----------
Overview
-----------

``pyrewton`` is a Python3 package for the identification of Carbohydrate Active enZymes (CAZymes) 
from candidate species, providing the user a complete CAZyome (all CAZymes encoded within a genome) 
for each candidate species, and the independent and comprehensive evaluation of the CAZyme prediction tools /
classifier [dbCAN](https://github.com/linnabrown/run_dbcan), 
[CUPP](https://www.bioengineering.dtu.dk/english/researchny/research-sections/section-for-protein-chemistry-and-enzyme-technology/enzyme-technology/cupp), 
and [eCAMI](https://github.com/zhanglabNKU/eCAMI).

``pyrewton`` provides a reproducible method for independently and comprehensively evaluating CAZyme classifiers. 
Specifically, evaluating the ability of CAZyme classifiers to:
- Differentiate between CAZymes and non-CAZymes
- Predict CAZy class annotations
- Predict multi-label CAZy class annotations
- Predict CAZy family annotations
- Predict multi-label CAZy family annotations

``pyrewton`` can also be used to compile a comprehensize CAZome database. ``pyrewton`` downloads all 
GenBank genomic assemblies for a list of candidate species provided as NCBI taxonomy IDs or scientific names. 
Protein sequences are then extracted from the assemblies and queried against a local CAZyme database (which **must** 
be created using [``cazy_webscraper``](https://hobnobmancer.github.io/cazy_webscraper/)). CAZy family annotations 
are retrieved from the local CAZyme database and stored in a local CAZome SQLite3 database, build by ``pyrewton``. 
Proteins extracted from the genomic assemblies and not found in the local CAZyme database are parsed by dbCAN to 
faciltiate annotating the comprehensive CAZome. ``pyrewton`` can add additional protein information (PDB accession, 
EC numbers and protein names) from [UniProt](https://www.uniprot.org/) to the local CAZome database, 
to create a comprehensive CAZome 
database.

Supplemenatry for the published independent evaluation of dbCAN, CUPP and eCAMI presented in Hobbs _et al._, 2021 is 
provided via `GitPages <https://hobnobmancer.github.io/pyrewton/>`_

------------
Contents
-------------

.. toctree::
   :maxdepth: 3

   evaluating
   test-sets
   trouble-shooting
   notebooks
   license

--------------
Citation
--------------

If you use ``pyrewton`` for indpendent benchmarking CAZyme classifers and/or compiling a local CAZome 
SQLite3 database, please cite our work:

   Hobbs, E. E. M., Gloster, T. M., Chapman, S., Pritchard, L. (2021): Microbiology Society Annual Conference 2021. figshare. Poster. https://doi.org/10.6084/m9.figshare.14370836.v

If using ``pyrewton`` for indpendent benchmarking CAZyme classifers do not forget to cite the tools you 
are evaluating.

---------------
Requirements
---------------

Python version 3.7+
Miniconda3 managed microenvironment, incorporated code checkers are included in list form in 'requirements.txt'.
Miniconda3 environment file is also available in the GitHub repository: 'environment.yml'.

--------------
Installation
--------------

1. Create a new virtual environment. (To install Conda please see the Conda `documentation <https://docs.conda.io/projects/conda/en/latest/user-guide/install/>`_)
.. code-block:: bash
   conda create -n pyrewton

2. Clone the ``pyrewton`` repository
.. code-block:: bash
   git clone https://github.com/HobnobMancer/pyrewton.git
   cd pyrewton  # navigate into the repo root

3. Install ``pyrewton`` 
.. code-block:: bash
   pip3 install -e .

4. Install the CAZyme classifiers into the correct directories
.. code-block:: bash
   python3 . cpt -p .

.. TIP::
   Following this method ensures all requirments are installed, as well as installing the CAZyme 
   classifers into the correct directories. This is essential becuase the CAZyme classifiers use 
   hard coded paths to access their respective datasets. ``pyrewton`` is required to navigate to the 
   correct directory to invoke these classifiers, so that the classifers can access their respective 
   datasets.

---------------
Notebooks
---------------

Jupyter notebook environments were crated, documenting how pyrewton was used during the EastBIO 2019-2023 PhD Project,
the GitHub pages for which are `available here <https://hobnobmancer.github.io/pyrewton/>`_. These can be used as examples for how to use pyrewton in research.

------------------------------
Help, Contribute and Support
------------------------------

Many of the common errors expected to arise during the operation of the scripts provided in this repository are
covered in this documentation, including the probable causes of these issues.

Please raise any issues with any of the programmers at the GitHub repository issues pages, by
following the `link <https://github.com/HobnobMancer/PhD_Project_Scripts/issues>`_.

.. note::
   pyrewton is still in development, and further functionalities are being added.
   Please see the `GitHub repository <https://github.com/HobnobMancer/pyrewton/tree/master>`_ for the latest developments.
