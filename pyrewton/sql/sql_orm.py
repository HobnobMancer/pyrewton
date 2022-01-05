#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) University of St Andrews 2020-2021
# (c) University of Strathclyde 2020-2021
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
#
"""Define the local CAZome db ORM"""


import logging
import re

import sqlite3

from sqlalchemy import (
    Boolean,
    Column,
    Float,
    ForeignKey,
    Index,
    Integer,
    PrimaryKeyConstraint,
    String,
    Table,
    UniqueConstraint,
    MetaData,
    create_engine,
    event,
    exc,
)
from sqlalchemy.engine import Engine
from sqlalchemy.ext.compiler import compiles
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relationship, sessionmaker
from sqlalchemy.sql.expression import BinaryExpression, func, literal
from sqlalchemy.sql.operators import custom_op


metadata_obj = MetaData()
Base = declarative_base()
Session = sessionmaker()


# Enable regular expression searching of the database
class ReString(String):
    """Enchanced version of standard SQLAlchemy's :class:`String`.

    Supports additional operators that can be used while constructing filter expressions.
    """
    class comparator_factory(String.comparator_factory):
        """Contains implementation of :class:`String` operators related to regular expressions."""
        def regexp(self, other):
            return RegexMatchExpression(self.expr, literal(other), custom_op('~'))

        def iregexp(self, other):
            return RegexMatchExpression(self.expr, literal(other), custom_op('~*'))

        def not_regexp(self, other):
            return RegexMatchExpression(self.expr, literal(other), custom_op('!~'))

        def not_iregexp(self, other):
            return RegexMatchExpression(self.expr, literal(other), custom_op('!~*'))


class RegexMatchExpression(BinaryExpression):
    """Represents matching of a column againsts a regular expression."""


@compiles(RegexMatchExpression, 'sqlite')
def sqlite_regex_match(element, compiler, **kw):
    """Compile the SQL expression representing a regular expression match for the SQLite engine."""
    # determine the name of a custom SQLite function to use for the operator
    operator = element.operator.opstring
    try:
        func_name, _ = SQLITE_REGEX_FUNCTIONS[operator]
    except (KeyError, ValueError) as e:
        would_be_sql_string = ' '.join((compiler.process(element.left),
                                        operator,
                                        compiler.process(element.right)))
        raise exc.StatementError(
            f"unknown regular expression match operator: {operator} {would_be_sql_string} {e}"
        )

    # compile the expression as an invocation of the custom function
    regex_func = getattr(func, func_name)
    regex_func_call = regex_func(element.left, element.right)
    return compiler.process(regex_func_call)


@event.listens_for(Engine, 'connect')
def sqlite_engine_connect(dbapi_connection, connection_record):
    """Listener for the event of establishing connection to a SQLite database.

    Creates the functions handling regular expression operators
    within SQLite engine, pointing them to their Python implementations above.
    """
    if not isinstance(dbapi_connection, sqlite3.Connection):
        return

    for name, function in SQLITE_REGEX_FUNCTIONS.values():
        dbapi_connection.create_function(name, 2, function)


# Mapping from the regular expression matching operators
# to named Python functions that implement them for SQLite.
SQLITE_REGEX_FUNCTIONS = {
    '~': ('REGEXP',
          lambda value, regex: bool(re.match(regex, value))),
    '~*': ('IREGEXP',
           lambda value, regex: bool(re.match(regex, value, re.IGNORECASE))),
    '!~': ('NOT_REGEXP',
           lambda value, regex: not re.match(regex, value)),
    '!~*': ('NOT_IREGEXP',
            lambda value, regex: not re.match(regex, value, re.IGNORECASE)),
}


#
# Relationship/Linker tables
#


# linker table between Proteins and PDB structures
proteins_pdbs = Table(
    "Proteins_Pdbs",
    Base.metadata,
    Column("protein_id", Integer, ForeignKey("Proteins.protein_id")),
    Column("pdb_id", Integer, ForeignKey("Pdbs.pdb_id")),
    PrimaryKeyConstraint("protein_id", "pdb_id"),
)


# linker table between Proteins and PDB structures
proteins_ecs = Table(
    "Proteins_Ecs",
    Base.metadata,
    Column("protein_id", Integer, ForeignKey("Proteins.protein_id")),
    Column("ec_id", Integer, ForeignKey("Ec_numbers.ec_id")),
    PrimaryKeyConstraint("protein_id", "ec_id"),
)


#
# Add tax tables
#


class Taxonomy(Base):
    """Describes the taxonomic data for a source organism. Data retrieved from NCBI"""
    __tablename__ = "Taxonomies"
    
    __table_args__ = (
        UniqueConstraint("genus", "species", "ncbi_tax_id"),
        Index("tax_index", "taxonomy_id", "ncbi_tax_id")
    )
    
    taxonomy_id = Column(Integer, primary_key=True)  # autoincremental
    ncbi_tax_id = Column(Integer)  # excludes the NCBI:txid prefix
    genus = Column(String)
    species = Column(String)

    assembly = relationship("Assembly", back_populates="tax_assembly")

    def __str__(self):
        return (
            f"-Source organism, genus={self.genus}, species={self.species}, "
            f"ncbi={self.ncbi_tax_id}, id={self.taxonomy_id}-"
        )

    def __repr__(self):
        return (
            f"<Class Taxonomy, genus={self.genus}, species={self.species}, "
            f"ncbi={self.ncbi_tax_id}, id={self.taxonomy_id}>"
        )


class Assembly(Base):
    """Represent a genomic assembly from which CAZymes (protein sequences) were retrieved."""
    __tablename__ = "Assemblies"
    
    __table_args__ = (
        UniqueConstraint("assembly_accession"),
        Index("assembly_index", "assembly_accession")
    )
    
    assembly_id = Column(Integer, primary_key=True)
    taxonomy_id = Column(Integer, ForeignKey("Taxonomies.taxonomy_id"))
    assembly_accession = Column(String)
    
    tax_assembly = relationship("Taxonomy", back_populates="assembly")
    assem_proteins = relationship("Protein", back_populates="assemblies")
    
    def __str__(self):
        return (
            f"-Genomic assembly, assembly_accession={self.assembly_accession}, id={self.assembly_id}-"
        )

    def __repr__(self):
        return (
            f"<Class Assembly: assembly_accession={self.assembly_accession}, "
            f"id={self.assembly_id}, tax={self.taxonomy_id}>"
        )


#
# add domain tabs
#


class Domain(Base):
    """Represent a domain/module within a protein."""
    __tablename__ = "Domains"
    
    __table_args__ = (
        Index("domain_index", "domain_id", "protein_id", "classifier_id", "family_id"),
    )
    
    domain_id = Column(Integer, primary_key=True)
    protein_id = Column(Integer, ForeignKey("Proteins.protein_id"))
    classifier_id = Column(Integer, ForeignKey("Classifiers.classifier_id"))
    family_id = Column(Integer, ForeignKey("CazyFamilies.family_id"))
    domain_range = Column(String)
    
    dom_proteins = relationship("Protein", back_populates="domains")
    classifier = relationship("Classifier", back_populates="parent_domain")
    family = relationship("CazyFamily", back_populates="parent_domain")
    
    def __str__(self):
        return f"-Domain, id={self.domain_id}, protein={self.protein_id}-"

    def __repr__(self):
        return f"<Domain, id={self.domain_id}, protein={self.protein_id}>"

        
class Classifier(Base):
    """Represent a CAZyme classifier.
    
    Classifiers: CAZy, dbCAN (consensus prediction), HMMER,
    Hotpep and DIAMOND
    """
    __tablename__ = "Classifiers"
    
    __table_args__ = (
        Index("classifier_index", "classifier_id", "classifier"),
    )
    
    classifier_id = Column(Integer, primary_key=True)
    classifier = Column(String)
    version = Column(String)
    cazy_training_set_date = Column(String)  # date of CAZy release
    
    parent_domain = relationship("Domain", back_populates="classifier")

    def __str__(self):
        return f"-Classifier, id={self.classifier_id}, name={self.classifier}-"

    def __repr__(self):
        return f"<Classifier, id={self.classifier_id}, name={self.classifier}>"

    
class CazyFamily(Base):
    """Represent the source genomic assembly of a protein."""
    __tablename__ = "CazyFamilies"
    
    __table_args__ = (
        UniqueConstraint("family"),
        Index("fam_index", "family_id", "family"),
    )
    
    family_id = Column(Integer, primary_key=True)
    family = Column(String)

    parent_domain = relationship("Domain", back_populates="family")

    def __str__(self):
        return f"-CAZy Fam, id={self.family_id}, fam={self.family}-"

    def __repr__(self):
        return f"<CAZy Fam, id={self.family_id}, fam={self.family}>"


#
# Protein table
#


class Protein(Base):
    """Represent a single protein."""
    __tablename__ = "Proteins"
    
    __table_args__ = (
        UniqueConstraint("protein_id"),
        Index("protein_index", "protein_id", "genbank_accession"),
    )
    
    protein_id = Column(Integer, primary_key=True)
    assembly_id = Column(Integer, ForeignKey("Assemblies.assembly_id"))
    genbank_accession = Column(String)
    protein_name = Column(String)  # retrieved from UniProt
    uniprot_id = Column(String)  # retrieved from UniProt
    mass = Column(Integer)  # retrieved from UniProt
    length = Column(Integer)  # length of Aa seq
    sequence = Column(String)
    
    # data from CAZy, dbCAN and genomes
    domains = relationship("Domain", back_populates="dom_proteins")
    assemblies = relationship("Assembly", back_populates="assem_proteins")
    # info about the possibility of transmembrane regions
    transmembranes = relationship("Transmembrane", back_populates="trans_proteins")
    # data retrieved from UniProt
    active_sites = relationship("ActiveSite", back_populates="act_proteins")
    substrate_sites = relationship("SubstrateBindingSite", back_populates="sub_proteins")
    metal_sites = relationship("MetalBindingSite", back_populates="met_proteins")
    cofactors = relationship("Cofactor", back_populates="cof_proteins")
    glycosylations = relationship("Glycosylation", back_populates="gly_proteins")
    temperatures = relationship("Temperature", back_populates="temp_proteins")
    optimal_phs = relationship("OptimalPH", back_populates="ph_proteins")
    data_citations = relationship("Citation", back_populates="cite_proteins")
    # many-to-many relationships
    pdbs = relationship(
        "Pdb",
        secondary=proteins_pdbs,
        back_populates="pdb_proteins",
        lazy="dynamic",
    )
    ec_numbers = relationship(
        "Ec_number",
        secondary=proteins_ecs,
        back_populates="ec_proteins",
        lazy="dynamic",
    )

    def __str__(self):
        return f"-Protein, id={self.protein_id}, name={self.protein_name}-"

    def __repr__(self):
        return f"<Protein, id={self.protein_id}, name={self.protein_name}>"


#
# Table for transmembrane data
#


class Transmembrane(Base):
    """Describe information about the possible presence of a transmembrane region in a protein"""
    __tablename__ = "Transmembranes"
    __table_args__ = (
        UniqueConstraint("transmembrane_id", "protein_id"),
        Index("transmembrane_id", "protein_id")
    )

    transmembrane_id = Column(Integer, primary_key=True)
    protein_id = Column(Integer, ForeignKey("Proteins.protein_id"))
    uniprot_transmembrane = Column(Boolean)
    tmhmm_transmembrane = Column(Integer)
    signal_p_transmembrane = Column(Boolean)

    trans_proteins = relationship("Protein", back_populates="transmembranes")

    def __str__(self):
        return (
            f"-Transmembrane region uniprot={self.uniprot_transmembrane},"
            f"tmhmm={self.tmhmm_transmembrane},signal_p={self.signal_p_transmembrane}-"
        )

    def __repr__(self):
        return (
            f"<Transmembrane region uniprot={self.uniprot_transmembrane},"
            f"tmhmm={self.tmhmm_transmembrane},signal_p={self.signal_p_transmembrane}>"
        )


#
# Add tabs for storing data from UniProt
#


class Pdb(Base):
    """Describe a PDB accession number of protein structure."""
    __tablename__ = "Pdbs"
    __table_args__ = (
        UniqueConstraint("pdb_accession"),
        Index('pdb_id', "pdb_accession")
    )

    pdb_id = Column(Integer, primary_key=True)
    pdb_accession = Column(String)

    pdb_proteins = relationship(
        "Protein",
        secondary=proteins_pdbs,
        back_populates="pdbs",
        lazy="dynamic",
    )

    def __str__(self):
        return f"-PDB accession={self.pdb_accession}, id={self.pdb_id}-"

    def __repr__(self):
        return f"<Class Pdb accession={self.pdb_accession}, id={self.pdb_id}>"


class Ec_number(Base):
    """Describe an EC number."""
    __tablename__ = "Ec_numbers"
    __table_args__ = (
        UniqueConstraint("ec_number"),
        Index('ec_id', "ec_number")
    )

    ec_id = Column(Integer, primary_key=True)
    ec_number = Column(String)

    ec_proteins = relationship(
        "Protein",
        secondary=proteins_ecs,
        back_populates="ec_numbers",
        lazy="dynamic",
    )

    def __str__(self):
        return f"-EC, id={self.ec_id}, number={self.ec_number}-"

    def __repr__(self):
        return f"<EC, id={self.ec_id}, number={self.ec_number}>"


class ActiveSite(Base):
    """Represent the active sites (Aa) in the enzymes."""
    __tablename__ = "ActiveSites"
    
    __table_args__ = (
        UniqueConstraint("protein_id", "activesite_id"),
        Index("active_site_index", "protein_id", "activesite_id", "position"),
    )
    
    protein_id = Column(Integer, ForeignKey("Proteins.protein_id"))
    activesite_id = Column(Integer, primary_key=True)
    position = Column(Integer)
    type_id = Column(Integer, ForeignKey("SiteTypes.type_id"))
    activity_id = Column(Integer, ForeignKey("AssociatedActivities.activity_id"))
    note = Column(ReString)
    evidence = Column(ReString)
    
    act_proteins = relationship("Protein", back_populates="active_sites")
    active_site_types = relationship("SiteType", back_populates="parent_active_sites")
    associated_activities = relationship("AssociatedActivity", back_populates="source_active_sites")
    
    def __str__(self):
        return f"-ActiveSite, protein={self.protein_id}, position={self.position}, id={self.activesite_id}-"

    def __repr__(self):
        return (
            f"<Class ActiveSite: protein={self.protein_id}, position={self.position},\n"
            f"id={self.activesite_id}, type_id={self.type_id}, activity_id={self.activity_id},\n"
            f"note={self.note}, evidence={self.evidence}>"
        )


class ActiveSiteType(Base):
    """Represent the type of active sites, e.g. nucleophile or base donor."""
    __tablename__ = "SiteTypes"
    
    __table_args__ = (
        UniqueConstraint("site_type"),
        Index("active_index", "site_type", "type_id")
    )
    
    type_id = Column(Integer, primary_key=True)
    site_type = Column(ReString)
    
    parent_active_sites = relationship("ActiveSite", back_populates="active_site_types")

    def __str__(self):
        return f"-ActiveSiteType, type_id={self.type_id}, site_type={self.site_type}-"

    def __repr__(self):
        return f"<ActiveSiteType, type_id={self.type_id}, site_type={self.site_type}>"


class AssociatedActivity(Base):
    """Represent the type of active sites, e.g. nucleophile or base donor."""
    __tablename__ = "AssociatedActivities"
    
    __table_args__ = (
        UniqueConstraint("associated_activity"),
        Index("associated_index", "associated_activity", "activity_id")
    )
    
    activity_id = Column(Integer, primary_key=True)
    associated_activity = Column(ReString)
    
    source_active_sites = relationship("ActiveSite", back_populates="associated_activities")

    def __str__(self):
        return (
            f"-AssociatedActivity, activity_id={self.activity_id}, associated_activity={self.associated_activity}-"
        )

    def __repr__(self):
        return (
            f"<AssociatedActivity, activity_id={self.activity_id}, associated_activity={self.associated_activity}>"
        )


class SubstrateBindingSite(Base):
    """Represent the substrate binding sites (Aa) in the enzymes."""
    __tablename__ = "SubstrateBindingSites"
    
    __table_args__ = (
        UniqueConstraint("protein_id", "substratesite_id"),
        Index("sub_site_index", "protein_id", "substratesite_id", "position")
    )
    
    protein_id = Column(Integer, ForeignKey("Proteins.protein_id"))
    substratesite_id = Column(Integer, primary_key=True)
    position = Column(Integer)
    note = Column(ReString)  # Typically the substrate, although may be burried within text
    evidence = Column(ReString)
    
    sub_proteins = relationship("Protein", back_populates="substrate_sites")
    
    def __str__(self):
        return f"-SubstrateBindingSite, protein={self.protein_id}, position={self.position}, id={self.substratesite_id}-"

    def __repr__(self):
        return (
            f"<SubstrateBindingSite, protein={self.protein_id}, position={self.position}, id={self.substratesite_id}>"
        )


class MetalBindingSite(Base):
    """Represent the metal binding sites (Aa) in the enzymes."""
    __tablename__ = "MetalBindingSites"
    
    __table_args__ = (
        UniqueConstraint("position", "ion", "ion_number"),
        Index("metal_site_index", "protein_id", "metalsite_id", "position")
    )
    
    protein_id = Column(Integer, ForeignKey("Proteins.protein_id"))
    metalsite_id = Column(Integer, primary_key=True)
    position = Column(Integer)
    ion = Column(Integer, ForeignKey("Metals.ion_id"))
    ion_number = Column(Integer)
    note = Column(ReString)  # Typically the substrate, although may be burried within text
    evidence = Column(ReString)
    
    met_proteins = relationship("Protein", back_populates="metal_sites")
    metals = relationship("Metal", back_populates="metal_sites")
    
    def __str__(self):
        return f"-MetalBindingSite, protein={self.protein_id}, position={self.position}, id={self.metalsite_id}-"

    def __repr__(self):
        return (
            f"<MetalBindingSite, protein={self.protein_id}, position={self.position}, id={self.metalsite_id}>"
        )

class Metal(Base):
    """Represent the metal binding sites (Aa) in the enzymes."""
    __tablename__ = "Metals"
    
    __table_args__ = (
        UniqueConstraint("ion"),
        Index("metal_index", "ion_id", "ion")
    )
    
    ion_id = Column(Integer, primary_key=True)
    ion = Column(String)
    
    metals = relationship("Metal", back_populates="metal_sites")
    
    def __str__(self):
        return f"-Metal, protein={self.ion_id}, ion={self.ion}-"

    def __repr__(self):
        return f"<Metal, protein={self.ion_id}, ion={self.ion}>"


class Cofactor(Base):
    """Represent the binding sites (Aa) of enzyme cofactors."""
    __tablename__ = "Cofactors"
    
    __table_args__ = (
        Index("cofactor_index", "protein_id", "cofactor_id", "molecule_id"),
    )
    
    protein_id = Column(Integer, ForeignKey("Proteins.protein_id"))
    cofactor_id = Column(Integer, primary_key=True)
    molecule_id = Column(Integer, ForeignKey("CofactorMolecules.molecule_id"))
    note = Column(ReString)  # Typically the substrate, although may be burried within text
    evidence = Column(ReString)
    
    cof_proteins = relationship("Protein", back_populates="cofactors")
    molecules = relationship("CofactorMolecule", back_populates="cofactor_site")
    
    def __str__(self):
        return f"-Cofactor, protein={self.protein_id}, cofactor_id={self.cofactor_id}-"

    def __repr__(self):
        return (
            f"<Cofactor, protein={self.protein_id}, cofactor_id={self.cofactor_id}, id={self.molecule_id}>"
        )

class CofactorMolecule(Base):
    """Represent enzyme, Cofactor molecules"""
    __tablename__ = "CofactorMolecules"
    
    __table_args__ = (
        UniqueConstraint("molecule"),
        Index("molecule_index", "molecule_id", "molecule"),
    )
    
    molecule_id = Column(Integer, primary_key=True)
    molecule = Column(String)

    cofactor_site = relationship("Cofactor", back_populates="molecules")
    
    def __str__(self):
        return f"-CofactorMolecule, id={self.molecule_id}, molecule={self.molecule}-"

    def __repr__(self):
        return f"<CofactorMolecule, cofactor_id={self.molecule_id}, cofactor={self.molecule}>"


class Glycosylation(Base):
    """Represent the glycosylation of an enzyme, by representing the sugar bound to an enzyme."""
    __tablename__ = "Glycosylations"
    
    __table_args__ = (
        UniqueConstraint("protein_id", "glycosylation_id"),
        Index("glycos_index", "protein_id", "glycosylation_id")
    )
    
    protein_id = Column(Integer, ForeignKey("Proteins.protein_id"))
    glycosylation_id = Column(Integer, primary_key=True)
    note = Column(ReString)  # Typically a description of the sugar architecture
    evidence = Column(ReString)
    
    gly_proteins = relationship("Protein", back_populates="glycosylations")
    
    def __str__(self):
        return f"-Glycosylation, protein={self.protein_id}-"

    def __repr__(self):
        return f"<Glycosylation, protein={self.protein_id}, note={self.note}>"


class Temperature(Base):
    """Represent the optimal, non-optimal and thermostable temperature ranges for an enzymes."""
    __tablename__ = "Temperatures"
    
    __table_args__ = (
        UniqueConstraint("protein_id", "temperature_id"),
        Index("temp_index", "protein_id", "temperature_id")
    )
    
    protein_id = Column(Integer, ForeignKey("Proteins.protein_id"))
    temperature_id = Column(Integer, primary_key=True)
    lower_optimal = Column(Float)
    upper_optimal = Column(Float)
    lower_thermostable = Column(Float)
    upper_thermostable = Column(Float)
    lower_lose_activity = Column(Float)  # Activity lose may start at the lower temp
    upper_lose_activity = Column(Float)  # complete or greater lose is achieved at the upper temp
    note = Column(ReString)  # A copy of the full text retrieved from UniProt, in case parsing fails
    evidence = Column(ReString)
    
    temp_proteins = relationship("Protein", back_populates="temperatures")
    
    def __str__(self):
        str_repr = f"-Temperature, protein={self.protein_id},"
        
        if self.lower_optimal is not None:
            str_repr += f" lower_optimal={self.lower_optimal}"
        if self.upper_optimal is not None:
            str_repr += f" upper_optimal={self.upper_optimal}"
        if self.lower_thermostable is not None:
            str_repr += f" lower_thermostable={self.lower_thermostable}"
        if self.upper_thermostable is not None:
            str_repr += f" upper_thermostable={self.upper_thermostable}"
        if self.lower_lose_activity is not None:
            str_repr += f" lower_lose_activity={self.lower_lose_activity}"
        if self.upper_lose_activity is not None:
            str_repr += f" upper_lose_activity={self.upper_lose_activity}"
        str_repr += "-"
        
        return str_repr

    def __repr__(self):

        str_repr = f"<Temperature, protein={self.protein_id} id={self.temperature_id}\n"
        
        if self.lower_optimal is not None:
            str_repr += f" lower_optimal={self.lower_optimal}"
        if self.upper_optimal is not None:
            str_repr += f" upper_optimal={self.upper_optimal}"
        if self.lower_thermostable is not None:
            str_repr += f" lower_thermostable={self.lower_thermostable}"
        if self.upper_thermostable is not None:
            str_repr += f" upper_thermostable={self.upper_thermostable}"
        if self.lower_lose_activity is not None:
            str_repr += f" lower_lose_activity={self.lower_lose_activity}"
        if self.upper_lose_activity is not None:
            str_repr += f" upper_lose_activity={self.upper_lose_activity}"
        str_repr += ">"
        
        return str_repr


class OptimalPH(Base):
    """Represent the optiml pH range of an enzyme."""
    __tablename__ = "OptimalPHs"
    
    __table_args__ = (
        UniqueConstraint("protein_id", "ph_id"),
        Index("ph_index", "protein_id", "ph_id")
    )
    
    protein_id = Column(Integer, ForeignKey("Proteins.protein_id"))
    ph_id = Column(Integer, primary_key=True)
    lower_pH = Column(Float)
    upper_pH = Column(Float)
    note = Column(ReString)  # entire string from UniProt in case parsing fails
    evidence = Column(ReString)
    
    ph_proteins = relationship("Protein", back_populates="optimal_phs")
    
    def __str__(self):
        return f"-Optimal pH range, protein={self.protein_id}-"

    def __repr__(self):
        return f"<OptimalPH, protein={self.protein_id}, note={self.note}>"


class Citation(Base):
    """Represent the literature sources of protein data."""
    __tablename__ = "Citations"
    
    __table_args__ = (
        UniqueConstraint("protein_id", "citation_id"),
        Index("citation_index", "protein_id", "citation_id")
    )
    
    protein_id = Column(Integer, ForeignKey("Proteins.protein_id"))
    citation_id = Column(Integer, primary_key=True)
    citation = Column(Integer)  # PubMed Id
    
    cite_proteins = relationship("Protein", back_populates="data_citations")
    
    def __str__(self):
        return f"-Citation, protein={self.protein_id}, id={self.citation_id}-"

    def __repr__(self):
        return (
            f"<Citation, protein={self.protein_id}, id={self.citation}, citation={self.citation}>"
        )


def get_cazome_db_connection(db_path, args):
    """Create open session to local CAZy SQL database.
    
    :param db_path: Path, path to local CAZome database
    
    Return an open database session.
    """
    logger = logging.getLogger(__name__)
    logger.info("Opening session to an existing local database")

    engine = create_engine(f"sqlite+pysqlite:///{db_path}", echo=args.sql_echo, future=True)
    Base.metadata.create_all(engine)
    Session.configure(bind=engine)  # allows for calls to session later on when required
    connection = engine.connect()
    
    return connection
