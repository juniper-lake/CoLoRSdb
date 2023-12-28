#!/usr/bin/env python3
"""Fix VCF ploidy and optionally anonymizes multi-sample VCF files by:
(1) randomizing the order of sample-specific data on a per-variant basis and
(2) changing the sample names to an arbitrary counter starting with the provided prefix.

It also adds a metadata line to the VCF file to reflect use of this script.
"""

__version__ = '0.1.0'

import argparse
import os
import sys
from codecs import open
from parser import VCFParser

from loguru import logger


def parse_args(args):
    """Parse command line parameters."""
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('vcf', help='optionally gzipped VCF file to be randomized')
    parser.add_argument(
        '--anonymize_prefix',
        dest='anonymize_prefix',
        help='sample name prefix for header, triggers anonymization',
    )
    parser.add_argument(
        '--non_diploid_regions',
        dest='non_diploid_regions_file',
        type=str,
        help='TSV file with haploid regions in 1-based closed [start, end] coordinates,\n'
        'columns include chrom, start, end, sex, ploidy, required for fixing ploidy',
    )
    parser.add_argument(
        '--sample_sexes',
        dest='sexes',
        type=str,
        nargs='+',
        help='List of sample_name+sex (e.g. sample1+M sample2+female), required for fixing ploidy',
    )
    parser.add_argument('--outfile', dest='outfile', help='output file name')
    parser.add_argument(
        '--loglevel',
        dest='loglevel',
        choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'],
        default='INFO',
        help='set the logging level',
    )
    parser.add_argument(
        '--version',
        action='version',
        version=f'%(prog)s (version {__version__})',
    )

    return parser.parse_args()


def get_non_diploid_regions(non_diploid_regions_file: str = None, sexes: list = None):
    """Read non diploid regions from file and return as list of tuples."""
    # check if hap regions file exists
    non_diploid_regions = []
    if non_diploid_regions_file:
        if not os.path.isfile(non_diploid_regions_file):
            raise OSError(f'File {non_diploid_regions_file} does not exist.')
        logger.info(f'Reading haploid regions from {non_diploid_regions_file}')
        with open(non_diploid_regions_file) as file:
            while line := file.readline():
                if line.startswith('#'):
                    continue
                line_list = line.rstrip().split()
                if len(line_list) != 5:
                    raise OSError(f'Non diploid file {non_diploid_regions_file} is not in the correct format.')
                else:
                    try:
                        non_diploid_regions.append(
                            (line_list[0], int(line_list[1]), int(line_list[2]), line_list[3], int(line_list[4]))
                        )
                        logger.debug(
                            f'Added region {line_list[0]}:{line_list[1]}-{line_list[2], line_list[3], line_list[4]}'
                        )
                    except ValueError:
                        raise OSError(
                            f'File {non_diploid_regions_file} with haploid regions is not in the correct format.'
                        )
            if not non_diploid_regions:
                raise OSError(f'File {non_diploid_regions_file} is empty.')
        if not sexes:
            raise OSError('File of haploid regions was provided but no sample sexes were specified.')
    return non_diploid_regions


def postprocess_joint_vcf(
    vcf: str,
    anonymize_prefix: str = None,
    non_diploid_regions_file: str = None,
    sexes: list = None,
    outfile: str = None,
    meta_string: str = 'postprocess_joint_vcf.py',
):
    """Parse VCF, change sample names, update metadata, and print shuffled VCF."""
    non_diploid_regions = get_non_diploid_regions(non_diploid_regions_file, sexes)

    # check if outfile exists
    if outfile:
        if not outfile.endswith('.vcf'):
            raise ValueError('Output file should end in .vcf')
        if os.path.isfile(outfile):
            raise OSError(f'Output file {outfile} already exists.')
        else:
            logger.info(f'Writing output VCF to {outfile}')
            file = open(outfile, 'w')
    else:
        logger.info('No output file provided, printing to stdout')
        file = sys.stdout

    vcf = VCFParser(vcf)

    if sexes:
        if not non_diploid_regions_file:
            raise OSError('Sample sexes were specified but no file of haploid regions was provided.')
        sexes_dict = dict(s.split('+') for s in sexes)
        sexes_sorted = [sexes_dict[sample] for sample in vcf.samples]

    # change sample names
    if anonymize_prefix:
        if (n_samples := len(vcf.samples)) < 2:
            raise OSError('This VCF only contains one sample, so anonymization is not useful.')
        else:
            vcf.change_sample_names([f'{anonymize_prefix}_{i}' for i in range(1, n_samples + 1)])

    # add metadata line
    vcf.add_meta('commandline', meta_string)

    # print header
    print(vcf.print_metadata(), file=file)

    # print headerline
    print(vcf.print_headerline(), file=file)

    # print each record with samples shuffled
    for variant in vcf:
        if non_diploid_regions_file and sexes:
            variant.fix_ploidy(sexes=sexes_sorted, non_diploid_regions=non_diploid_regions)
        if anonymize_prefix:
            variant.shuffle_samples()
        print(variant, file=file)

    logger.success('Successfully terminated')


if __name__ == '__main__':
    """This is executed when run from the command line."""

    args = parse_args(sys.argv[1:])

    logger.configure(handlers=[{'sink': sys.stderr, 'level': args.loglevel}])

    postprocess_joint_vcf(
        args.vcf,
        args.anonymize_prefix,
        args.non_diploid_regions_file,
        args.sexes,
        args.outfile,
        ' '.join(sys.argv),
    )
