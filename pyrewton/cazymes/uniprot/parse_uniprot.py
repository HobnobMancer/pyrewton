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
