#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author:
# Emma E. M. Hobbs
#
# Contact
# eemh1@st-andrews.ac.uk
#
# Emma E. M. Hobbs,
# Biomolecular Sciences Building,
# University of St Andrews,
# North Haugh Campus,
# St Andrews,
# KY16 9ST
# Scotland,
# UK
#
# The MIT License
"""Evaluate accuracy of binary classification of CAZyme and non-CAZymes
by dbCAN, CUPP and eCAMI.

:cmd_args:

:func __:
"""
from sklearn import metrics

labels_true = [0, 0, 0, 1, 1]
labels_pred = [0, 0, 1, 0, 1]

result = metrics.adjusted_rand_score(labels_true, labels_pred)
print(result)

result_1 = metrics.f1_score(labels_true, labels_pred, average="micro")
print(result_1)
