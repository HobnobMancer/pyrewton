.. pyrewton documentation master file, created by
   sphinx-quickstart on Fri Jun 19 16:19:25 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

===================================
Welcome to pyrewton's documentation
===================================

.. image:: docu_cover.jpg
   :width: 60%
   :align: left
   :alt: Logos of EastBIO, BBSRC, University of St Andrews and James Huttton Institute

| 
| Version v0.1.1 2020/06/04
| DOI: 10.5281/zendo.3876218 
| `GitHub repository <https://github.com/HobnobMancer/PhD_Project_Scripts/tree/master>`_
| 



Overview
========

pyrewton is a Python3 programme, run at the command line and free to use under the MIT License.
The programme can be used to identifying carbohydrate processing enzyme engineering candidates
for advanced biocatalysis. Specifically, pyrewton supports:

Downloading of all genomic assemblies (as GenBank files .gbff) from the
`NCBI Assembly database <https://www.ncbi.nlm.nih.gov/assembly>`_ associated with each species passed to the programme

Summarising the cazyme annotations within GenBank files

Retrieving all cazymes from `UniProtKB <https://www.uniprot.org/>`_ for each species
passed to the script

*More detailed documentation for each module is linked to in the contents table below, including links to documentation to help with trouble shooting.*

Contents
--------

.. toctree::
   :maxdepth: 2

   genbank
   annotations
   trouble-shooting
   notebooks

Installation
============

| The easiest way to install pyrewton is to use ``pip``:
| ``pip3 install -e <path to pyrewton setup.py file>``.
| This will install all required Python packages and dependencies.

Requirements
------------

Python version 3.7+
Miniconda3 managed microenvironment, incorporated code checkers are included in list form in 'requirements.txt'.
Miniconda3 environment file is also available in the GitHub repository: 'environment.yml'.

Notebooks
=========

Jupyter notebook environments were crated, documenting how pyrewton was used during the EastBIO 2019-2023 PhD Project,
the GitHub pages for which are `available here <https://hobnobmancer.github.io/PhD_Project_Scripts/>`_. These can be used as examples for how to use pyrewton in research.

Contribute and Support
======================

Many of the common errors expected to arise during the operation of the scripts provided in this repository are
covered in this documentation, including the probable causes of these issues.

Please raise any issues with any of the programmers at the GitHub repository issues pages, by
following the `link <https://github.com/HobnobMancer/PhD_Project_Scripts/issues>`_.

.. note::
   pyrewton is still in development, and further functionalities are being added.
   Please see the `GitHub repository <https://github.com/HobnobMancer/PhD_Project_Scripts/tree/master>`_ for the latest developments.
