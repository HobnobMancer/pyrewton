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


import logging
import re

import numpy as np


def get_sites_data(results_row, column_name, line_starter, get_data_func, protein_id):
    """Generic func for retrieving data for a specific item from a UniProt results df.
    
    :param results_row: pandas df, df of UniProt query results
    :param column_name: str, name of the column to retrieve data from
    :param line_starter: str, substring to identify new data in a cell
    :param get_data_func: func, func for retrieving and parsing data
    :param protein_id: int, local CAZome db ID for the source protein
    
    Return set(), one item per row to be inserted
    """
    sites = set()

    value = results_row[column_name]

    try:
        if np.isnan(value):
            return sites  # no data to retrieve
    except TypeError:
        pass

    data = value.split(";")

    index = 0
    data_index = 0
    for index in range(len(data)):

        if data[data_index].strip().startswith(line_starter):  # new row to be inserted in the db
            site, data_index = get_data_func(data, data_index, protein_id)
            sites.add(site)

        if data_index == len(data):
            break

    return sites


def get_substrate_binding_site_data(data, data_index, protein_id):
    """"""
    position = data[data_index].strip()
    position = int(position.replace("BINDING ", ""))
    
    substrate_binding_site = {
        'position': position,
        'substrate': None,
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
                if substrate_binding_site['substrate'] is None:
                    substrate_binding_site['substrate'] = note
                else:
                    new_note = f"{substrate_binding_site['substrate']} {new_data}"
                    substrate_binding_site['substrate'] = new_note
            
            else:
                substrate_binding_site['substrate'] = note
        
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
                if substrate_binding_site['substrate'] is None:
                    substrate_binding_site['substrate'] = new_data
                else:
                    new_note = f"{substrate_binding_site['substrate']} {new_data}"
                    substrate_binding_site['substrate'] = new_note
                
        index += 1
        if index == len(data):
            break  # came to end of list
    
    substrate_binding_tuple = (
        protein_id,
        substrate_binding_site['position'],
        substrate_binding_site['substrate'],
        substrate_binding_site['evidence'],
    )
    
    return substrate_binding_tuple, index


def get_glycosylation_data(data, data_index, protein_id):
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
        protein_id,
        glycosylation_site['position'],
        glycosylation_site['note'],
        glycosylation_site['evidence'],
    )
    
    return glycosylation_tuple, index


def get_temp_data(results_row, protein_id):
    temp_dependence_data = set()  # tuples (lower stable, upper stable, lower loss of activity, upper loss of activity, note, evidence)

    value = results_row['Temperature dependence']
    try:
        if np.isnan(value):
            return temp_dependence_data  # no data to retrieve
    except TypeError:
        pass

    if value.startswith('BIOPHYSICOCHEMICAL PROPERTIES') is False:
        return temp_dependence_data  # no data to retrieve
    
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
        return temp_dependence_data  # no data to retrieve
    
    # retrieve the evidence
    evidence = re.search(r"\{\.+?\})", value).group()
    if evidence == "":
        evidence = None
    else:
        evidence = evidence.replace("{","")
        evidence = evidence.replace("}","").strip()
    
    note = value[:value.find("{")]

    temp_dependence_data.add(
        (
            protein_id,
            lower_optimum,
            upper_optimum,
            lower_stable,
            upper_stable,
            lower_loss,
            upper_loss,
            note,
            evidence,
        )
    )

    return temp_dependence_data


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


def get_optimum_ph(results_row, protein_id):
    """Retrieve the optimum pH for the protein's activity."""
    optimum_ph_data = set()  # set of tuples (pH, note, evidence)

    value = results_row['pH dependence']

    try:
        if np.isnan(value):
            return optimum_ph_data
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
                    ph = re.search(r"\d?", optimum).group()
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
            optimum_ph_data.add( (protein_id, float(lower_ph), float(upper_ph), note, evidence) )

    return optimum_ph_data


def get_citations(results_row, protein_id):
    """Retrieve publication citations for the protein."""
    citations = set()  # tuples, one tuple per PubMed ID

    for value in results_row['PubMed ID']:

        try:
            if np.isnan(value):
                continue
        except TypeError:
            pass

        try:
            data = value.split(";")
        except AttributeError:
            citations.add((value,))
            continue

        index = 0

        for pubmed_id in data:
            try:
                pubmed_id = int(pubmed_id)
                citations.add( (protein_id, pubmed_id,) )
            except ValueError:
                continue
            

    return citations


def get_active_site_data(results_row, protein_id):
    active_sites = set()  # store tuples (aa position, evidence)
    site_types = set()
    activities = set()
    
    value = results_row['Active site']

    try:
        if np.isnan(value):
            return active_sites, site_types, activities
    except TypeError:
        pass

    data = value.split(";")

    # items group in pairs, the first is the Aa site position, the second the evidence
    index = 0
    data_index = 0

    while index < len(data):
        if data[data_index].strip().startswith("ACT_SITE"):  # new metal binding site
            active_site, new_site_types, new_activities, data_index = extract_active_site_data(
                data,
                data_index,
                protein_id,
            )
            active_sites.add(active_site)
            site_types = site_types.union(new_site_types)
            activities = activities.union(new_activities)

        if data_index == len(data):
            break
            
    return active_sites, site_types, activities


def extract_active_site_data(data, data_index, protein_id):
    position = data[data_index].strip()
    position = int(position.replace("ACT_SITE ", "").strip())
    
    active_site = {
        'position': position,
        'type': None,  # e.g. Nucleotphile or Proton donor
        'activity': None,  # e.g. for esterase activity
        'note': None,  # misc
        'evidence': None,
    }
    
    site_types = set()
    activities = set()
    
    # check if any more data relating to this metal binding site
    data_index += 1
    
    index = data_index
    limit_index = data_index
    
    for limit_index in range(len(data[data_index:])):
        new_data = data[index].strip()
        
        if new_data.startswith("ACT_SITE"):
            # new active site
            break
        
        elif new_data.startswith("/note"): 
            note = new_data.replace('/note="', '')
            site_type = note.replace('"', '')
            
            active_site['type'] = site_type
            site_types.add( (site_type,) )
            
        elif new_data.startswith("for") or new_data.endswith("activity"):
            activity = new_data.replace('"', '')
            activity = activity.replace("for", "")
            
            active_site['activity'] = activity.strip()
            activities.add( (activity,) )
        
        elif new_data.startswith("/evidence"):
            evidence = new_data.replace('/evidence="', '')
            evidence = evidence.replace('"', '')
            
            if active_site['evidence'] is None:
                active_site['evidence'] = evidence
            else:
                new_evidence = f"{active_site['evidence']} {evidence}"
                active_site['evidence'] = new_evidence
        
        else:  # part of note that covered multiple lines
            new_data = new_data.replace('"', '')
            
            if len(new_data) != 0:
                if active_site['note'] is None:
                    active_site['note'] = new_data
                else:
                    new_note = f"{active_site['note']} {new_data}"
                    active_site['note'] = new_note
                
        index += 1
        if index == len(data):
            break  # came to end of list
    
    active_site_tuple = (
        protein_id,
        active_site['position'],
        active_site['type'],
        active_site['activity'],
        active_site['note'],
        active_site['evidence'],
    )
    
    return active_site_tuple, site_types, activities, index


def get_metal_binding_sites(results_row, protein_id):
    metal_binding_sites = set()  # store tuples (position, ion, ion_number, note, evidence)
    metals = set()  # store tuples (ion,)

    value = results_row['Metal binding']

    try:
        if np.isnan(value):
            return metal_binding_sites, metals
    except TypeError:
        pass

    data = value.split(";")

    index = 0
    data_index = 0
    for index in range(len(data)):

        if data[data_index].strip().startswith("METAL"):  # new metal binding site
            new_site, new_metals, data_index = get_metal_binding_data(data, data_index, protein_id)
            metal_binding_sites.add(new_site)
            metals = metals.union(new_metals)

        if data_index == len(data):
            break

    return metal_binding_sites, metals


def get_metal_binding_data(data, data_index, protein_id):
    logger = logging.getLogger(__name__)

    position = data[data_index].strip()
    position = int(position.replace("METAL ", ""))
    
    metal_binding_site = {
        'position': position,
        'ion': None,  # name of metal ion, e.g. Calcium
        'ion_number': None,  # if multiple ions of the same name can bind, this distinguishes between them
        'note': None,
        'evidence': None,
    }
    metals = set()
    
    # check if any more data relating to this metal binding site
    data_index += 1
    
    index = data_index
    limit_index = data_index
    
    for limit_index in range(len(data[data_index:])):
        new_data = data[index].strip()
        
        if new_data.startswith("METAL"):
            # new binding site
            break
        
        elif new_data.startswith("/note"): 
            note = new_data.replace('/note="', '')
            
            ion = note.split(" ")[0]
            ion = ion.replace('"', '')
            try:
                ion_number = note.split(" ")[1]
                ion_number = int(ion_number.replace('"', ''))
            except IndexError:
                ion_number = None
            except ValueError:
                if ion_number.startswith("/note"):
                    ion = ion_number.replace("/note", "").strip()
                else:
                    if metal_binding_site['note'] is None:
                        metal_binding_site['note'] = ion_number
                    else:
                        metal_binding_site['note'] += f" {ion_number}"
            
            metal_binding_site['ion'] = ion
            metals.add( (ion,) )
        
        elif new_data.startswith("/evidence"):
            evidence = new_data.replace('/evidence="', '')
            evidence = evidence.replace('"', '')
            
            if metal_binding_site['evidence'] is None:
                metal_binding_site['evidence'] = evidence
            else:
                new_evidence = f"{metal_binding_site['evidence']} {evidence}"
                metal_binding_site['evidence'] = new_evidence
        
        else:  # part of note that covered multiple lines
            new_data = new_data.replace('"', '')
            
            if len(new_data) != 0:
                if metal_binding_site['note'] is None:
                    metal_binding_site['note'] = new_data
                else:
                    new_note = f"{metal_binding_site['note']} {new_data}"
                    metal_binding_site['note'] = new_note

        index += 1
        if index == len(data):
            break  # came to end of list
    
    metal_binding_tuple = (
        protein_id,
        metal_binding_site['position'],
        metal_binding_site['ion'],
        metal_binding_site['ion_number'],
        metal_binding_site['note'],
        metal_binding_site['evidence'],
    )
    
    return metal_binding_tuple, metals, index


def get_cofactor_data(results_row, protein_id):
    cofactor_data = get_sites_data(
        results_row,
        'Cofactor',
        'COFACTOR',
        get_cofactors,
        protein_id,
    )

    new_molecules = set()

    for data_tuple in cofactor_data:
        new_molecules.add( (data_tuple[1],) )

    return cofactor_data, new_molecules


def get_cofactors(data, data_index, protein_id):
    cofactor_molecule = data[data_index].strip()
    cofactor_molecule = cofactor_molecule.replace("COFACTOR: Name=", "")
    
    cofactor = {
        'cofactor': cofactor_molecule,
        'note': None,
        'evidence': None,
    }
    
    # check if any more data relating to this metal binding site
    data_index += 1
    
    index = data_index
    limit_index = data_index
    
    for limit_index in range(len(data[data_index:])):
        new_data = data[index].strip()
        
        if new_data.startswith("COFACTOR"):
            # new binding site
            break
        
        elif new_data.startswith("/note"): 
            note = new_data.replace('/note="', '')
            note = note.replace('"', '')
            
            if cofactor['note'] is None:
                cofactor['note'] = note
            else:
                new_note = f"{cofactor['note']} {new_data}"
                cofactor['note'] = new_note
        
        elif new_data.startswith("/evidence") or new_data.startswith('Xref'):
            evidence = new_data.replace('/evidence="', '')
            evidence = evidence.replace('Xref=', '')
            evidence = evidence.replace('"', '')
            
            if cofactor['evidence'] is None:
                cofactor['evidence'] = evidence
            else:
                new_evidence = f"{cofactor['evidence']} {evidence}"
                cofactor['evidence'] = new_evidence
        
        else:  # part of note that covered multiple lines
            new_data = new_data.replace('"', '')
            if len(new_data) != 0:
                if cofactor['note'] is None:
                    cofactor['note'] = new_data
                else:
                    new_note = f"{cofactor['note']} {new_data}"
                    cofactor['note'] = new_note
                
        index += 1
        if index == len(data):
            break  # came to end of list
    
    cofactor_tuple = (
        protein_id,
        cofactor['cofactor'],
        cofactor['note'],
        cofactor['evidence'],
    )
    
    return cofactor_tuple, index


def get_pdb_ecs(results_row, protein_id, column_name):
    value = results_row[column_name]

    accessions = set()
    protein_relationships = set()

    try:
        if np.isnan(value):
            return protein_relationships, accessions
    except TypeError:
        pass

    data = value.split(";")

    for item in data:
        if item != "":
            accessions.add( (item.strip(),) )
            protein_relationships.add( (protein_id, item.strip()) )

    return protein_relationships, accessions
