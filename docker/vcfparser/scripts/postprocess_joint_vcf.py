#!/usr/bin/env python3
"""
This script fixes VCF ploidy and optionally anonymizes multi-sample VCF files by:
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
    parser.add_argument(
        "--anonymize_prefix",
        dest="anonymize_prefix",
        help="sample name prefix for header, triggers anonymization",
    )
    parser.add_argument(
        "--hap_bed",
        dest="hap_bed",
        type=str,
        help="BED file with haploid regions, required for fixing ploidy",
    )
    parser.add_argument(
        "--sample_sexes",
        dest="sexes",
        type=str,
        nargs="+",
        help="List of sample_name+sex (e.g. sample1+M sample2+female), \
            required for fixing ploidy",
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


def postprocess_joint_vcf(
    vcf: str,
    anonymize_prefix: str = None,
    hap_bed: str = None,
    sexes: list = None,
    outfile: str = None,
    meta_string: str = "anonymize_joint_vcf.py",
):
    """Parse VCF, change sample names, update metadata, and print shuffled VCF"""

    # check if outfile exists
    if outfile:
        if not outfile.endswith(".vcf"):
            raise ValueError("Output file should end in .vcf")
        if os.path.isfile(outfile):
            raise IOError(f"Output file {outfile} already exists.")
        else:
            logger.info(f"Writing anonymized VCF to {outfile}")
            file = open(outfile, "w")
    else:
        logger.info("No output file provided, printing to stdout")
        file = sys.stdout

    # check if bed file exists
    if hap_bed:
        if not os.path.isfile(hap_bed):
            raise IOError(f"BED file {hap_bed} does not exist.")
        hap_regions = []
        logger.info(f"Reading haploid regions from {hap_bed}")
        with open(hap_bed) as file:
            while line := file.readline():
                if line.startswith("#"):
                    continue
                line_list = line.rstrip().split("\t")
                if len(line_list) != 3:
                    raise IOError(f"BED file {hap_bed} is not in the correct format.")
                else:
                    try:
                        hap_regions.append(
                            (line_list[0], int(line_list[1]), int(line_list[2]))
                        )
                        logger.debug(
                            f"Added haploid region \
                                {line_list[0]}:{line_list[1]}-{line_list[2]}"
                        )
                    except ValueError:
                        raise IOError(
                            f"BED file {hap_bed} is not in the correct format."
                        )
            if not hap_regions:
                raise IOError(f"BED file {hap_bed} is empty.")
        if not sexes:
            raise IOError(
                "BED file of haploid regions was provided but no \
                          sample sexes were specified."
            )

    vcf = VCFParser(vcf)

    if sexes:
        if not hap_bed:
            raise IOError(
                "Sample sexes were specified but no BED file of haploid \
                          regions was provided."
            )
        sexes_dict = dict(s.split("+") for s in sexes)
        sexes_sorted = [sexes_dict[sample] for sample in vcf.samples]

    # change sample names
    if anonymize_prefix:
        if (n_samples := len(vcf.samples)) < 2:
            raise IOError(
                "This VCF only contains one sample, so anonymization is not useful."
            )
        else:
            vcf.change_sample_names(
                [f"{anonymize_prefix}_{i}" for i in range(1, n_samples + 1)]
            )

    # add metadata line
    vcf.add_meta("commandline", meta_string)

    # print header
    print(vcf.print_header(), file=file)

    # print each record with samples shuffled
    for variant in vcf:
        if hap_regions and sexes:
            variant.fix_ploidy(sexes=sexes_sorted, regions=hap_regions)
        if anonymize_prefix:
            variant.shuffle_samples()
        print(variant, file=file)

    logger.success("Successfully terminated")


if __name__ == "__main__":
    """This is executed when run from the command line."""

    args = parse_args(sys.argv[1:])

    logger.configure(handlers=[{"sink": sys.stderr, "level": args.loglevel}])

    postprocess_joint_vcf(
        args.vcf,
        args.anonymize_prefix,
        args.hap_bed,
        args.sexes,
        args.outfile,
        " ".join(sys.argv),
    )
