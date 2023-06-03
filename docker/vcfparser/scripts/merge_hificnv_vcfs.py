#!/usr/bin/env python3
"""
This script aggregates multiple HiFiCNV output VCFs into a single 
gzipped TSV and optionally anonymizes the data by:
(1) randomizing the order of sample-specific data on a per-variant basis and
(2) changing the sample names to an arbitrary counter starting with the provided prefix.

It also adds a metadata line to the VCF file to reflect use of this script.

VCFs should be sorted and optionally gzipped.
"""

__version__ = "0.1.0"

from loguru import logger
import argparse
from codecs import open
from parser import VCFParser
import sys
import os
import random


def parse_args(args):
    """Parse command line parameters"""
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        "vcfs",
        nargs="*",
        help="list of optionally gzipped VCF files to be aggregated and anonymized",
    )
    parser.add_argument(
        "--anonymize_prefix",
        dest="anonymize_prefix",
        help="sample name prefix for header, triggers anonymization",
    )
    parser.add_argument("--outfile", dest="outfile", help="output file name")
    parser.add_argument(
        "--loglevel",
        dest="loglevel",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        default="INFO",
        help="set the logging level",
    )
    parser.add_argument(
        "--version",
        action="version",
        version="%(prog)s (version {version})".format(version=__version__),
    )

    return parser.parse_args()


def merge_hificnv_vcfs(
    vcf_files: list,
    anonymize_prefix: str = None,
    outfile: str = None,
    meta_string: str = "merge_hificnv_vcfs.py",
):
    """Combine multiple VCFs into single VCF and optionally anonymize samples"""

    # check if outfile exists
    if outfile:
        if not outfile.endswith(".vcf"):
            raise ValueError("Output file should end in .vcf")
        if os.path.isfile(outfile):
            raise IOError(f"Output file {outfile} already exists.")
        else:
            file = open(outfile, "w")
    else:
        logger.info("No output file provided, printing to stdout")
        file = sys.stdout

    if len(vcf_files) < 2:
        raise ValueError("Please provide at least two VCF files to aggregate")

    logger.success("Successfully terminated")


if __name__ == "__main__":
    """This is executed when run from the command line."""

    args = parse_args(sys.argv[1:])

    logger.configure(handlers=[{"sink": sys.stderr, "level": args.loglevel}])

    merge_hificnv_vcfs(
        args.vcfs, args.anonymize_prefix, args.outfile, " ".join(sys.argv)
    )
