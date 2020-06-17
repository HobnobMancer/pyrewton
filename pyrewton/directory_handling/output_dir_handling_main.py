#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author:
# Emma E. M. Hobbs

# Contact
# eemh1@st-andrews.ac.uk

# Emma E. M. Hobbs,
# Biomolecular Sciences Building,
# University of St Andrews,
# North Haugh Campus,
# St Andrews,
# KY16 9ST
# Scotland,
# UK

# The MIT License
"""Contains generic functions for handling output directories.

:func make_output_directory: create directory for files to be written to
:func write_out_dataframe: write out data frame to csv file
"""

import shutil
import sys


def make_output_directory(output, logger, force, nodelete):
    """Create output directory for genomic files.

    Check if directory indicated for output existed already.
    If so check if force overwrite enabled. If not terminate programme.
    If so, check if deletion of exiting files was enabled.
    If so, exiting files in output directory are deleted.
    Create output directory, expecting error if already exists.

    :param output: Path, path to output directory
    :param logger: logger object
    :param force: bool, cmd-line args to enable/disable over writing existing directory
    :param nodelete: bool, cmd-line args to enable/disable deleting of existing files

    Return nothing.
    """
    logger.info("Checking if specified output directory exists.")
    # If output directory specificed at cmd-line, check output directory does not already exist
    if output.exists():
        if force is False:
            logger.critical(
                (
                    "Output directory already exists and forced overwrite not enabled.\n"
                    "Terminating program."
                )
            )
            sys.exit(1)
        else:
            if nodelete is False:
                logger.warning(
                    (
                        "Output directory already exists and forced complete overwrite enabled.\n"
                        "Deleting existing content in outdir."
                    )
                )
                # delete existing content in outdir
                shutil.rmtree(output)
            else:
                logger.warning(
                    (
                        "Output directory already exists and forced addition of files "
                        "to outdir enables."
                    )
                )
    # Recursively make output directory (directory may exist if
    # --force==True and --nodelete==True, so exist_ok is required)
    output.mkdir(exist_ok=True)
    return


def write_out_dataframe(species_table, logger, outdir, force, nodelete):
    """Write out dataframe to output directory.

    :param species_table: pandas dataframe
    :param logger: logger object
    :param outdir: cmd-args, Path, output directory
    :param force: booleon, cmd-line argument to enable/disable over writing of existing files
    :param nodelete: boolean, cmd-line args to enable/disable deleting of existing files in outdir

    return Nothing.
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

    species_table.to_csv(outdir)

    return ()
