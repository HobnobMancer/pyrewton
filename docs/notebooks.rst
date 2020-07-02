
=========
Notebooks
=========

Jupyter notebook environments have been used to facilitate explanation of the function and operation of pyrewton. These notebooks explicitly describe how the scripts were invoked to during the EastBIO PhD project (‘Identificating Engineering Candidate for Advanced Biocatalysis in Biofuel Production’), and were written to supplement the thesis. However, these notebooks can also be used as an example as to invoking ``pyrewton`` to identify cazyme canidated for engineering.
The notebooks all containing the final output generated when invoked during the EastBIO PhD project.

.. note::
    These notebooks contain exerts of code to facilitate the explanation of the code code architecture and function, as well as details provided on how the code was implemented during the project. Therefore, the code in many of the code cells is not runnable and replication or exploration of the data should be performed using the original scripts provided within the repository.

Accessing the notebooks
-----------------------

The notebooks can be accessed via the GitHub pages created for the host `repository <https://hobnobmancer.github.io/PhD_Project_Scripts/>`_, or directly via the list below.

Accessing the notebooks in browser
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

`01_downloading_genbank_files <https://hobnobmancer.github.io/PhD_Project_Scripts/Notebooks/01_Downloading_GenBank_Files.html>`_
This refers to invoking get_ncbi_genomes to download all directly linked GenBank files for each species passed to the program.

`02_retrieving_genbank_annotations <https://hobnobmancer.github.io/PhD_Project_Scripts/Notebooks/02_CAZy_GenBank_Summary.html>`_
This refers to invoking get_cazyme_annotations to retrieve all cazyme annotations from GenBank files.


Accessing the notebooks via the terminal
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

From the top level of the repository directory on your system, start the Jupyter notebook server by issuing the command: ``Jupyter notebook``. This will open a new browser window or tab containing the Jupyter homepage, with a listing of all files and directories as laid out in the repository.
Navigate to the ‘notebooks/’ directory, containing all Jupyter notebooks for the repository. To open a notebook simple click on the title of the notebook.
For more information please see the Jupyter notebook `Quick-start guide <https://jupyter-notebook-beginner-guide.readthedocs.io/en/latest/?fbclid=IwAR1yIwkYCDjcw5FJZ7CfKES3l72HubqGYGcFrVrUKwWZoYh4NHy3VVu0AgQ>`_ and this `tutorial <https://www.tutorialspoint.com/jupyter/jupyter_quick_guide.htm>`_.
