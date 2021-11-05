
============================================
Independent Evaluation of CAZyme Classifiers
============================================

At the present ``pyrewton`` supports the independent and comprehensive evaluation of the CAZyme 
classifiers `dbCAN <https://academic.oup.com/nar/article/46/W1/W95/4996582>`_, 
`CUPP <https://biotechnologyforbiofuels.biomedcentral.com/articles/10.1186/s13068-019-1436-5>`_, and 
`eCAMI <https://pubmed.ncbi.nlm.nih.gov/31794006/>`_

Evaluation (or benchmarking) CAZyme classifiers is handled by the ``cazymes/evaluate_tools`` submodule.

The evaluation process is split up into 4 stage:
1. Creating the test sets :ref:`test-set-label`
2. Invoking the CAZyme classifiers :ref:`test-set-label`
3. Calculating statistical parameters :ref:`test-set-label`
4. Presenting the results :ref:`test-set-label`

Supplementary material for `Hobbs et al., 2021 <https://figshare.com/articles/poster/Microbiology_Society_Annual_Conference_2021/14370836>`_ 
is hosted in the GitHub repository to faciltiate reproduction of the presented evaluation, as well 
as be used as template for future evaluations. The supplementary material can be found 
 `here <https://github.com/HobnobMancer/pyrewton/tree/master/supplementary>`_ .
