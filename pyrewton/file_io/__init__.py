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
"""Manage input and output directory handling."""

import shutil


def make_output_directory(output, logger, force, nodelete):
    """Create output directory for genomic files.

    :param args: Namespace with cmd-line arguments
    :param logger: logger object

    args should contain attributes: output (path to output directory), nodelete (boolean
    to indicate whether an existing output directory should be removed, if present), and
    force (boolean to indicate whether an existing directory should be overwritten)

    Raises FileExistsError if an attempt is made to create a directory that already
    exists, without a suitable args.force/args.nodelete combination
    """
    if force is True:
        logger.warning(
            "Output directory %s exists, nodelete is %s", output, nodelete,
        )
        if nodelete and output.exists():
            logger.warning("Not deleting directory %s", output)
        elif output.exists():
            logger.warning("Deleting directory %s", output)
            shutil.rmtree(output)

    logger.info("Creating directory %s", output)
    try:
        output.mkdir(exist_ok=force)
    except FileExistsError:
        logger.warning(
            "Out directory already exists. New directory not made, writing to existing directory"
        )


def write_out_dataframe(dataframe, logger, outdir, force):
    """Write out dataframe to output directory.

    :param species_table: pandas dataframe
    :param logger: logger object
    :param outdir: cmd-args, Path, output directory
    :param force: booleon, cmd-line argument to enable/disable over writing of existing files
    """
    # Check if overwrite of existing directory will occur
    logger.info("Checking if output directory for dataframe already exists")
    if outdir.exists():
        if force is False:
            logger.warning(
                "Specified directory for dataframe already exists.\nExiting writing out dataframe."
            )
            return ()
        else:
            logger.warning(
                "Specified directory for dataframe already exists.\nForced overwritting enabled."
            )
    logger.info("Writing out species dataframe to directory")

    dataframe.to_csv(outdir)


def write_out_pre_named_dataframe(dataframe, df_name, logger, outdir, force):
    """Write out dataframe to output directory.

    :param dataframe: pandas dataframe
    :param df_name: str, name of dataframe
    :param logger: logger object
    :param outdir: cmd-args, Path, output directory
    :param force: booleon, cmd-line argument to enable/disable over writing of existing files
    """
    # Check if overwrite of existing directory will occur
    logger.info("Checking if output directory for dataframe already exists")
    output_path = outdir / f"{df_name}.csv"
    if output_path.exists():
        if force is False:
            logger.warning(
                "Specified directory for dataframe already exists.\nExiting writing out dataframe."
            )
            return ()
        else:
            logger.warning(
                "Specified directory for dataframe already exists.\nForced overwritting enabled."
            )
    logger.info("Writing out species dataframe to directory")

    dataframe.to_csv(output_path)
