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
"""Download genomic assemblies from NCBI"""


import logging
import re
import time

from pathlib import Path
from socket import timeout
from tqdm import tqdm
from urllib.error import HTTPError, URLError
from urllib.request import urlopen

from Bio import Entrez


def get_genomic_assembly(assembly_accession, suffix="genomic.gbff.gz"):
    """Coordinate downloading Genomic assemmbly from the NCBI Assembly database.

    :param assembly_accession: str, accession of the Genomic assembly to be downloaded
    :param suffix: str, suffix of file

    Return path to downloaded genomic assembly.
    """
    # compile url for download
    genbank_url, filestem = compile_url(assembly_accession, suffix)

    # create path to write downloaded Genomic assembly to
    out_file_path = "_".join([filestem.replace(".", "_"), suffix])
    out_file_path = Path(out_file_path)

    # download GenBank file
    download_file(genbank_url, out_file_path, assembly_accession, "GenBank file")

    return out_file_path


def compile_url(accession_number, suffix, ftpstem="ftp://ftp.ncbi.nlm.nih.gov/genomes/all"):
    """Retrieve URL for downloading the assembly from NCBI, and create filestem of output file path

    :param accession_number: str, asseccion number of genomic assembly

    Return str, url required for download and filestem for output file path for the downloaded assembly.
    """
    # search for the ID of the record
    with entrez_retry(
        10, Entrez.esearch, db="Assembly", term="GCA_000021645.1[Assembly Accession]", rettype='uilist',
    ) as handle:
        search_record = Entrez.read(handle)

    # retrieve record for genomic assembly
    with entrez_retry(
        10,
        Entrez.esummary,
        db="assembly",
        id=search_record['IdList'][0],
        report="full",
    ) as handle:
        record = Entrez.read(handle)

    assembly_name = record["DocumentSummarySet"]["DocumentSummary"][0]["AssemblyName"]

    escape_characters = re.compile(r"[\s/,#\(\)]")
    escape_name = re.sub(escape_characters, "_", assembly_name)

    # compile filstem
    filestem = "_".join([accession_number, escape_name])

    # separate out filesteam into GCstem, accession number intergers and discarded
    url_parts = tuple(filestem.split("_", 2))

    # separate identifying numbers from version number
    sub_directories = "/".join(
        [url_parts[1][i : i + 3] for i in range(0, len(url_parts[1].split(".")[0]), 3)]
    )

    # return url for downloading file
    return (
        "{0}/{1}/{2}/{3}/{3}_{4}".format(
            ftpstem, url_parts[0], sub_directories, filestem, suffix
        ),
        filestem,
    )


def download_file(genbank_url, out_file_path, accession_number, file_type):
    """Download file.

    :param genbank_url: str, url of file to be downloaded
    :param out_file_path: path, output directory for file to be written to
    :param accession_number: str, accession number of genome
    :param file_type: str, denotes in logger file type downloaded

    Return nothing.
    """
    logger = logging.getLogger(__name__)
    # Try URL connection
    try:
        response = urlopen(genbank_url, timeout=45)
    except (HTTPError, URLError, timeout) as e:
        logger.error(
            f"Failed to download {file_type} for {accession_number}", exc_info=1,
        )
        return

    if out_file_path.exists():
        logger.warning(f"Output file {out_file_path} exists, not downloading")
        return

    # Download file
    file_size = int(response.info().get("Content-length"))
    bsize = 1_048_576
    try:
        with open(out_file_path, "wb") as out_handle:
            # Using leave=False as this will be an internally-nested progress bar
            with tqdm(
                total=file_size,
                leave=False,
                desc=f"Downloading {accession_number} {file_type}",
            ) as pbar:
                while True:
                    buffer = response.read(bsize)
                    if not buffer:
                        break
                    pbar.update(len(buffer))
                    out_handle.write(buffer)
    except IOError:
        logger.error(f"Download failed for {accession_number}", exc_info=1)
        return

    return


def entrez_retry(retries, entrez_func, *func_args, **func_kwargs):
    """Call to NCBI using Entrez.

    :param retries: int, maximum number of retries excepted if network error encountered
    :param entrez_func: function, call method to NCBI
    :param *func_args: tuple, arguments passed to Entrez function
    :param ** func_kwargs: dictionary, keyword arguments passed to Entrez function

    Returns record.
    """
    logger = logging.getLogger(__name__)
    record, retries, tries = None, retries, 0

    while record is None and tries < retries:
        try:
            record = entrez_func(*func_args, **func_kwargs)

        except IOError:
            # log retry attempt
            if tries < retries:
                logger.warning(
                    f"Network error encountered during try no.{tries}.\nRetrying in 10s",
                    exc_info=1,
                )
                time.sleep(10)
            tries += 1

    if record is None:
        logger.error(
            "Network error encountered too many times. Exiting attempt to call to NCBI"
        )
        return

    return record
