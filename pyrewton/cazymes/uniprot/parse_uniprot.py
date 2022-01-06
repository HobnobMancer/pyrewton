#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) University of St Andrews 2020-2021
# (c) University of Strathclyde 2020-2021
# (c) James Hutton Institute 2020-2021
#
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
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
"""Retrieve data from UniProt for CAZymes in a CAZome database"""


import re

import numpy as np


def get_substrate_binding_site_data(data, data_index):
    position = data[data_index].strip()
    position = int(position.replace("BINDING ", ""))
    
    substrate_binding_site = {
        'position': position,
        'note': None,
        'evidence': None,
    }
    
    # check if any more data relating to this metal binding site
    data_index += 1
    
    index = data_index
    limit_index = data_index
    
    for limit_index in range(len(data[data_index:])):
        new_data = data[index].strip()
        
        if new_data.startswith("BINDING"):
            # new binding site
            break
        
        elif new_data.startswith("/note"): 
            note = new_data.replace('/note="', '')
            note = note.replace('"', '')
            
            if note != "Substrate":
                if substrate_binding_site['note'] is None:
                    substrate_binding_site['note'] = note
                else:
                    new_note = f"{substrate_binding_site['note']} {new_data}"
                    substrate_binding_site['note'] = new_note
        
        elif new_data.startswith("/evidence"):
            evidence = new_data.replace('/evidence="', '')
            evidence = evidence.replace('"', '')
            
            if substrate_binding_site['evidence'] is None:
                substrate_binding_site['evidence'] = evidence
            else:
                new_evidence = f"{substrate_binding_site['evidence']} {evidence}"
                substrate_binding_site['evidence'] = new_evidence
        
        else:  # part of note that covered multiple lines
            new_data = new_data.replace('"', '')
            
            if len(new_data) != 0:
                if substrate_binding_site['note'] is None:
                    substrate_binding_site['note'] = new_data
                else:
                    new_note = f"{substrate_binding_site['note']} {new_data}"
                    substrate_binding_site['note'] = new_note
                
        index += 1
        if index == len(data):
            break  # came to end of list
    
    substrate_binding_tuple = (
        substrate_binding_site['position'],
        substrate_binding_site['note'],
        substrate_binding_site['evidence'],
    )
    
    return substrate_binding_tuple, index


def get_glycosylation_data(data, data_index):
    position = data[data_index].strip()
    position = int(position.replace("CARBOHYD ", ""))
    
    glycosylation_site = {
        'position': position,
        'note': None,
        'evidence': None,
    }
    
    # check if any more data relating to this metal binding site
    data_index += 1
    
    index = data_index
    limit_index = data_index
    
    for limit_index in range(len(data[data_index:])):
        new_data = data[index].strip()
        
        if new_data.startswith("CARBOHYD"):
            # new binding site
            break
        
        elif new_data.startswith("/note"): 
            note = new_data.replace('/note="', '')
            note = note.replace('"', '')
            
            if glycosylation_site['note'] is None:
                glycosylation_site['note'] = note
            else:
                new_note = f"{glycosylation_site['note']} {new_data}"
                glycosylation_site['note'] = new_note
        
        elif new_data.startswith("/evidence"):
            evidence = new_data.replace('/evidence="', '')
            evidence = evidence.replace('"', '')
            
            if glycosylation_site['evidence'] is None:
                glycosylation_site['evidence'] = evidence
            else:
                new_evidence = f"{glycosylation_site['evidence']} {evidence}"
                glycosylation_site['evidence'] = new_evidence
        
        else:  # part of note that covered multiple lines
            new_data = new_data.replace('"', '')
            if len(new_data) != 0:
                if glycosylation_site['note'] is None:
                    glycosylation_site['note'] = new_data
                else:
                    new_note = f"{glycosylation_site['note']} {new_data}"
                    glycosylation_site['note'] = new_note
                
        index += 1
        if index == len(data):
            break  # came to end of list
    
    glycosylation_tuple = (
        glycosylation_site['position'],
        glycosylation_site['note'],
        glycosylation_site['evidence'],
    )
    
    return glycosylation_tuple, index


def get_temp_range(value, term):
    """Return the temp range listed for the type of term, identified by the 'term'"""
    lower_temp, upper_temp = None, None
    # 
    try:
        temp = re.search(rf'{term} (above|over) (\d+?\.\d+?) degrees Celsius', value).group()
        lower_temp = temp
        upper_temp = None

    except AttributeError:
        try:
            temp = re.search(rf'{term} (above|over) (\d+?) degrees Celsius', value).group()
            lower_temp = temp
            upper_temp = None
        
        except AttributeError:
            try:
                temp = re.search(rf'{term} below (\d+?\.\d+?) degrees Celsius', value).group()
                lower_temp = None
                upper_temp = temp
        
            except AttributeError:
                try:
                    temp = re.search(rf'{term} below (\d+?) degrees Celsius', value).group()
                    lower_temp = None
                    upper_temp = temp
                
                except AttributeError:
                    try:
                        temp = re.search(rf'{term} is (((\d+?|\d+?\.\d+?))|(\d+?|\d+?\.\d+?)(-| to )(\d+?|\d+?\.\d+?)) degrees Celsius', value).group()
                        if temp.find("to") != -1:
                            lower_temp = temp.split(" to ")[0]
                            upper_temp = temp.split(" to ")[1]
                        
                        else:
                            lower_temp = temp.split("-")[0]
                            upper_temp = temp.split("-")[1]
                            
                    except AttributeError:
                        try:
                            temp = re.search(rf'{term} is (\d+?\.\d+?) degrees Celsius', value).group()
                            lower_temp = temp
                            upper_temp = temp
                            
                        except AttributeError:
                            try:
                                temp = re.search(rf'{term} is (\d+?) degrees Celsius', value).group()
                                lower_temp = temp
                                upper_temp = temp
                            
                            except AttributeError:
                                lower_temp, upper_temp = None, None

    return lower_temp, upper_temp


def get_temp_data(results_table):
    temp_dependence_data = set()  # tuples (lower stable, upper stable, lower loss of activity, upper loss of activity, note, evidence)

    for value in results_table['Temperature dependence']:
        try:
            if np.isnan(value):
                continue
        except TypeError:
            pass

        if value.startswith('BIOPHYSICOCHEMICAL PROPERTIES') is False:
            continue
        
        # attempt to retrieve optimal temp range
        lower_optimum, upper_optimum = get_temp_range(value, "Optimum temperature is ")

        # attempt to retrieve thermostable range
        lower_stable, upper_stable = get_temp_range(value, "Thermostable to ")
        
        # attempt to retrieve loss of activity range
        lower_loss, upper_loss = get_temp_range(value, "Complete loss of activity ")
        
        if lower_optimum is None and \
        upper_optimum is None and \
        lower_stable is None and \
        upper_stable is None and \
        lower_loss is None and \
        upper_loss is None:
            continue
        
        # retrieve the evidence
        evidence = re.search(r"\{\.+?\})", value).group()
        if evidence == "":
            evidence = None
        else:
            evidence = evidence.replace("{","")
            evidence = evidence.replace("}","").strip()
        
        note = value[:value.find("{")]

        temp_dependence_data.add(
            (lower_optimum, upper_optimum, lower_stable, upper_stable, lower_loss, upper_loss, note, evidence)
        )
    
    return temp_dependence_data


def get_optimum_ph(results_table):
    """Retrieve the optimum pH for the protein's activity."""
    optimum_ph_data = set()  # set of tuples (pH, note, evidence)

    for value in results_table['pH dependence']:
        try:
            if np.isnan(value):
                continue
        except TypeError:
            pass

        if value.startswith('BIOPHYSICOCHEMICAL PROPERTIES'):
            optimum_phs = re.findall(r'Optimum pH is .*? {.*?}', value)

            for optimum in optimum_phs:
                ph = None
                
                try:
                    ph = re.search(r"(((\d?\.\d?)|\d?)-((\d?\.\d?)|\d?))", optimum).group()
                    lower_ph = ph.split("-")[0]
                    upper_ph = ph.split("-")[1]
                    
                except AttributeError:  # raised if not a range, but a single temp is given
                    ph = re.search(r"(\d?\.\d?)", optimum).group()
                    if ph == '':  # if a int not a float is provided, an empty str is returned
                        ph = re.search(r"\d?", optimum_phs_8[0]).group()
                    lower_ph = ph
                    upper_ph = ph
                
                if ph is not None:
                    if optimum.find("with") != -1:
                        note = re.search(r"with .*? {", optimum).group()
                        note = note.replace("{", "").strip()
                        if note.endswith("."):
                            note = note[:-2]
                    else:
                        note = None
                        
                    evidence = optimum[optimum.find("{")+1:optimum.find("}")]
                    
                else:
                    continue
                
                # not all pH values are strored as floats in UniProt
                # standardised the datatype for the local db
                optimum_ph_data.add( (float(lower_ph), float(upper_ph), note, evidence) )

    return optimum_ph_data
