#!/usr/bin/env python3
"""
This script aggregates multiple TRGT output VCFs into a single 
gzipped TSV and optionally anonymizes the data by:
(1) randomizing the order of sample-specific data on a per-variant basis and
(2) changing the sample names to an arbitrary counter starting with the provided prefix.

It also adds a metadata line to the VCF file to reflect use of this script.

TRGT VCFs should be sorted and optionally gzipped.
"""

__version__ = "0.1.0"

from loguru import logger
import argparse
from codecs import open
from parser import VCFParser
import sys
import os
import random


def catch(func, *args, exception, handle, message, **kwargs):
    """Execute a function and raise an exception if it fails"""
    try:
        return func(*args, **kwargs)
    except exception:
        raise handle(message) from None


def release(func, *args, exception, handle, message, **kwargs):
    """Execute a function and raise an exception if it succeeds"""
    try:
        func(*args, **kwargs)
        raise handle(message)
    except exception:
        pass


def flatten(lst):
    """Flatten a list of lists"""
    return [item for sublist in lst for item in sublist]


def parse_args(args):
    """Parse command line parameters"""
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        "vcfs",
        nargs="*",
        help="list of optionally gzipped TRGT VCF files to be aggregated and anonymized",
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


def merge_trgt_vcfs(
    vcf_files: list,
    anonymize_prefix: str = None,
    outfile: str = None,
    meta_string: str = "merge_trgt_vcfs.py",
):
    """Combine multiple TRGT VCFs into single VCF and optionally anonymize samples"""

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

    # parse all VCFs, create iterators for all vcfs except the first
    vcfs = []
    vcf_iterators = []
    for vcf_file in vcf_files:
        vcfs.append(VCFParser(vcf_file))
        if len(vcfs) > 1:
            vcf_iterators.append(iter(vcfs[-1]))

    # make sure metadata are the same
    for vcf in vcfs[1:]:
        if vcf.metadata != vcfs[0].metadata:
            raise ValueError("VCFs have different metadata")

    # print metadata and of first vcf
    vcfs[0].add_meta("commandline", meta_string)
    print(vcfs[0].print_metadata(), file=file)

    # print headerline
    samples = vcfs[0].samples + flatten([vcf.samples for vcf in vcfs[1:]])
    if anonymize_prefix:
        samples = [f"{anonymize_prefix}_{i}" for i in range(1, len(samples) + 1)]
    headerline = "\t".join(vcfs[0].header_columns + samples)
    print(headerline, file=file)

    # iterate through first vcf, get corresponding records from other vcfs
    for sample1_record in vcfs[0]:
        # go to next record for each vcf iterator, if any are empty, raise IndexError
        other_sample_records = [
            catch(
                next,
                vcf_iterator,
                exception=StopIteration,
                handle=IndexError,
                message="VCFs have different numbers of variants",
            )
            for vcf_iterator in vcf_iterators
        ]

        # make sure info and format fields are the same for each variant
        for other_sample_record in other_sample_records:
            if sample1_record.info != other_sample_record.info:
                raise ValueError(
                    "VCFs have different INFO fields for the same locus. Were the VCFs sorted?"
                )
            if sample1_record.format != other_sample_record.format:
                raise ValueError(
                    "VCFs have different INFO fields for the same locus. Were the VCFs sorted?"
                )

        alts = sample1_record.alts + flatten(
            [other_sample_record.alts for other_sample_record in other_sample_records]
        )
        alts = sorted(list(set(alts)))

        # if not all ref calls, update genotypes to reflect full list of alt alleles
        if "." in alts and len(alts) > 1:
            alts.remove(".")
            sample1_record.update_genotypes(ref=sample1_record.ref, alts=alts)
            [
                other_sample_record.update_genotypes(ref=sample1_record.ref, alts=alts)
                for other_sample_record in other_sample_records
            ]
            sample_data = sample1_record.sample_data + [
                other_sample_record.sample_data[0]
                for other_sample_record in other_sample_records
            ]
            if anonymize_prefix:
                sample_data = random.sample(sample_data, len(sample_data))

            # ['CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT',samples]
            aggregated_record = [
                sample1_record.chrom,
                sample1_record.pos,
                sample1_record.id,
                sample1_record.ref,
                ",".join(alts),
                ".",
                ".",
                sample1_record.info,
                sample1_record.format,
            ] + sample_data
            print("\t".join(aggregated_record), file=file)

    # make sure the first vcf doesn't have fewer variants than the others
    for vcf_iterator in vcf_iterators:
        release(
            next,
            vcf_iterator,
            exception=StopIteration,
            handle=IndexError,
            message="VCFs have different numbers of variants",
        )

    logger.success("Successfully terminated")


if __name__ == "__main__":
    """This is executed when run from the command line."""

    args = parse_args(sys.argv[1:])

    logger.configure(handlers=[{"sink": sys.stderr, "level": args.loglevel}])

    merge_trgt_vcfs(
        args.vcfs, args.anonymize_prefix, args.outfile, " ".join(sys.argv)
    )
