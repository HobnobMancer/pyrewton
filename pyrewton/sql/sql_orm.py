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
    Column,
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


class Ncbi_Tax_Id(Base):
    """Describes an NCBI taxonomy ID."""
    __tablename__ = "Ncbi_Tax_Ids"
    __table_args__ = (
        UniqueConstraint("ncbi_tax_id"),
        Index("ncbi_id", "ncbi_tax_id")
    )

    ncbi_id = Column(Integer, primary_key=True)  # autoincremental
    ncbi_tax_id = Column(String)

    species = relationship("Taxonomy", back_populates="ncbi_ids")

    def __str__(self):
        return f"-NCBI, NCBI:txid={self.ncbi_tax_id}, id={self.ncbi_id}-"

    def __repr__(self):
        return f"<NCBI, NCBI:txid={self.ncbi_tax_id}, id={self.ncbi_id}>"


class Taxonomy(Base):
    """Describes the taxonomic data for a source organism. Data retrieved from NCBI"""
    __tablename__ = "Taxs"
    
    __table_args__ = (
        UniqueConstraint("genus", "species", "ncbi_tax_id"),
        Index("tax_index", "taxonomy_id", "ncbi_tax_id")
    )
    
    taxonomy_id = Column(Integer, primary_key=True)  # autoincremental
    ncbi_id = Column(Integer, ForeignKey("Ncbi_Tax_Ids.ncbi_id"))
    genus = Column(String)
    species = Column(String)
    ncbi_tax_id = Column(Integer)  # excludes the NCBI:txid prefix

    ncbi_ids = relationship("Ncbi_Tax_id", back_populates="species")
    assembly = relationship("Assembly", back_populates="tax_assembly")

    def __str__(self):
        return f"-Source organism, genus={self.genus}, species={self.species}, ncbi={self.ncbi_tax_id}, id={self.taxonomy_id}-"

    def __repr__(self):
        return f"<Class Taxonomy, genus={self.genus}, species={self.species}, ncbi={self.ncbi_tax_id}, id={self.taxonomy_id}>"


class Assembly(Base):
    """Represent a genomic assembly from which CAZymes (protein sequences) were retrieved."""
    __tablename__ = "Asssemblies"
    
    __table_args__ = (
        UniqueConstraint("assembly_accession"),
        Index("assembly_index", "assembly_accession")
    )
    
    assembly_id = Column(Integer, primary_key=True)
    assembly_accession = Column(String)
    taxonomy_id = Column(Integer, ForeignKey("Taxs.taxonomy_id"))
    
    tax_assembly = relationship("Taxonomy", back_populates="assembly")
    
    def __str__(self):
        return (
            f"-Genomic assembly, assembly_accession={self.assembly_accession}, id={self.assembly_id}-"
        )

    def __repr__(self):
        return (
            f"<Class Assembly: assembly_accession={self.assembly_accession}, id={self.assembly_id}, tax={self.taxonomy_id}>"
        )


def get_db_connection():
    return
