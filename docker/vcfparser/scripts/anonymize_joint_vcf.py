#!/usr/bin/env python3
"""
This script anonymizes multi-sample VCF files by:
(1) randomizing the order of sample-specific data on a per-variant basis and
(2) changing the sample names to an arbitrary counter starting with the provided prefix.

It also adds a metadata line to the VCF file to reflect use of this script.
"""

__version__ = "0.1.0"

from loguru import logger
import sys
import argparse
import os
from codecs import open
from parser import VCFParser


def parse_args(args):
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument("vcf", help="optionally gzipped VCF file to be randomized")
    parser.add_argument("cohort_prefix", help="sample name prefix for header line")
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


def anonymize_join_vcf(
    vcf: str,
    cohort_prefix: str,
    outfile: str = None,
    meta_string: str = "anonymize_joint_vcf.py",
):
    """Parse VCF, change sample names, update metadata, and print shuffled VCF"""

    if outfile:
        if os.path.isfile(outfile):
            raise IOError(f"Output file {outfile} already exists.")

    vcf = VCFParser(vcf)

    # change sample names
    if (n_samples := len(vcf.samples)) < 2:
        raise IOError(
            "This VCF file only contains one sample, so anonymization is not useful."
        )
    else:
        vcf.change_sample_names(
            [f"{cohort_prefix}_{i}" for i in range(1, n_samples + 1)]
        )

    # add metadata line
    vcf.add_meta("commandline", meta_string)

    if outfile:
        logger.info(f"Writing anonymized VCF to {outfile}")
        with open(outfile, mode="w", encoding="utf-8", errors="strict") as out:
            # print metadata and header line
            out.write(vcf.print_header())

            # print each record with samples shuffled
            for variant in vcf(shuffle_samples=True):
                out.write(str(variant) + "\n")
    else:
        # print header
        print(vcf.print_header())

        # print each record with samples shuffled
        for variant in vcf(shuffle_samples=True):
            print(variant)

    logger.success("Successfully terminated")


if __name__ == "__main__":
    """This is executed when run from the command line."""

    args = parse_args(sys.argv[1:])

    logger.configure(handlers=[{"sink": sys.stderr, "level": args.loglevel}])

    anonymize_join_vcf(args.vcf, args.cohort_prefix, args.outfile, " ".join(sys.argv))
