#!/usr/bin/env python3
"""Aggregate multiple TRGT output VCFs into a single
gzipped TSV and optionally anonymizes the data by:
(1) randomizing the order of sample-specific data on a per-variant basis and
(2) changing the sample names to an arbitrary counter starting with the provided prefix.

It also adds a metadata line to the VCF file to reflect use of this script.

TRGT VCFs should be sorted and optionally gzipped.
"""

__version__ = '0.1.0'

import argparse
import os
import random
import sys
from codecs import open
from parser import VCFParser

from loguru import logger


def catch(func, *args, exception, handle, message, **kwargs):
    """Execute a function and raise an exception if it fails."""
    try:
        return func(*args, **kwargs)
    except exception:
        raise handle(message) from None


def release(func, *args, exception, handle, message, **kwargs):
    """Execute a function and raise an exception if it succeeds."""
    try:
        func(*args, **kwargs)
        raise handle(message)
    except exception:
        pass


def flatten(lst):
    """Flatten a list of lists."""
    return [item for sublist in lst for item in sublist]


def parse_args(args):
    """Parse command line parameters."""
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument(
        'vcfs',
        nargs='*',
        help='list of optionally gzipped TRGT VCF files to be aggregated and anonymized',
    )
    parser.add_argument(
        '--bed', dest='bed', help='BED file used to generate VCF, should be sorted the same way as VCFs.'
    )
    parser.add_argument(
        '--anonymize_prefix',
        dest='anonymize_prefix',
        help='sample name prefix for header, triggers anonymization',
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


def get_alleles(samples, vcf_records, bed_trgt_id, bed_start):
    """Get ref and alts, check that TRID matches and POS is as expected."""
    # get ref allele
    if int(vcf_records[0].pos) == (bed_start + 1):
        ref = vcf_records[0].ref
    elif int(vcf_records[0].pos) == bed_start:
        ref = vcf_records[0].ref[1:]

    # get alt alleles
    alts = []
    for sample_idx in range(len(vcf_records)):
        if (vcf_trgt_id := vcf_records[sample_idx].info.split(';')[0].split('=')[1]) != bed_trgt_id:
            raise ValueError(f'{samples[sample_idx]} VCF TRID {vcf_trgt_id} does not match BED TRID {bed_trgt_id}.')
        if int(vcf_records[sample_idx].pos) == (bed_start + 1):
            alts = alts + vcf_records[sample_idx].alts
        # for alleles of length zero (i.e. entire ref range deleted, the POS is the same as the BED start instead of start + 1)
        elif int(vcf_records[sample_idx].pos) == bed_start:
            # update vcf_record alts to remove leading base, replace complete deletions with N
            replacement_alts = [alt[1:] if (len(alt) > 1) else 'N' for alt in vcf_records[sample_idx].alts]
            vcf_records[sample_idx].update_alleles(ref=ref, alts=replacement_alts)
            alts = alts + replacement_alts
        else:
            raise ValueError(
                f'Sample {samples[sample_idx]} VCF position {vcf_records[sample_idx].pos} for TRGT ID {bed_trgt_id} does not match BED position {bed_start}'
            )

    alts = sorted(list(set(alts)))

    return ref, alts


def update_genotypes(ref, alts, vcf_records, bed_trgt_id, samples):
    """Update genotypes to reflect full list of alt alleles."""
    # if there's more than one alt (either real or "."), update genotype numbers to reflect full list of alt alleles
    if len(alts) > 1:
        if '.' in alts:
            alts.remove('.')
        for sample_idx in range(len(vcf_records)):
            logger.debug(f'Updating genotypes for {samples[sample_idx]} at trgt id {bed_trgt_id}')
            vcf_records[sample_idx].update_genotypes(ref=ref, alts=alts)
    sample_data = [vcf_record.sample_data[0] for vcf_record in vcf_records]
    return sample_data


def merge_trgt_vcfs(
    vcf_files: list,
    bed: str,
    anonymize_prefix: str = None,
    outfile: str = None,
    meta_string: str = 'merge_trgt_vcfs.py',
):
    """Combine multiple TRGT VCFs into single VCF and optionally anonymize samples."""
    # check if outfile exists
    if outfile:
        if not outfile.endswith('.vcf'):
            raise ValueError('Output file should end in .vcf')
        if os.path.isfile(outfile):
            raise OSError(f'Output file {outfile} already exists.')
        else:
            file = open(outfile, 'w')
    else:
        logger.info('No output file provided, printing to stdout')
        file = sys.stdout

    if len(vcf_files) < 2:
        raise ValueError('Please provide at least two VCF files to aggregate')

    # parse all VCFs, create iterators for all vcfs
    vcfs = []
    vcf_iterators = []
    for vcf_file in vcf_files:
        vcfs.append(VCFParser(vcf_file))
        vcf_iterators.append(iter(vcfs[-1]))

    # print metadata of first vcf
    vcfs[0].add_meta('commandline', meta_string)
    print(vcfs[0].print_metadata(), file=file)

    # print headerline
    samples = flatten([vcf.samples for vcf in vcfs])
    if anonymize_prefix:
        samples = [f'{anonymize_prefix}_{i}' for i in range(1, len(samples) + 1)]
    headerline = '#' + '\t'.join(vcfs[0].header_columns + samples)
    print(headerline, file=file)

    # iterate through bed, get corresponding records from vcfs
    with open(bed) as bed_file:
        while trgt_bed_line := bed_file.readline():
            bed_trgt_id = trgt_bed_line.split('\t')[3].split(';')[0].split('=')[1]
            bed_chrom = trgt_bed_line.split('\t')[0]
            bed_start = int(trgt_bed_line.split('\t')[1])

            # go to next record for each vcf iterator, if any are empty, raise IndexError
            vcf_records = [
                catch(
                    next,
                    vcf_iterator,
                    exception=StopIteration,
                    handle=IndexError,
                    message='VCFs have different numbers of variants',
                )
                for vcf_iterator in vcf_iterators
            ]

            # get ref and alts
            ref, alts = get_alleles(samples, vcf_records, bed_trgt_id, bed_start)

            logger.info(f'Aggregating VCF records corresponding to TRGT ID {bed_trgt_id}')

            sample_data = update_genotypes(ref, alts, vcf_records, bed_trgt_id, samples)

            if anonymize_prefix:
                sample_data = random.sample(sample_data, len(sample_data))

            # ['CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT',samples]
            aggregated_record = [
                bed_chrom,
                str(bed_start + 1),
                '.',
                ref,
                ','.join(alts),
                '.',
                '.',
                vcf_records[0].info,
                vcf_records[0].format,
            ] + sample_data
            print('\t'.join(aggregated_record), file=file)

    # make sure the first vcf doesn't have fewer variants than the others
    for vcf_iterator in vcf_iterators:
        release(
            next,
            vcf_iterator,
            exception=StopIteration,
            handle=IndexError,
            message='VCFs have different numbers of variants',
        )

    logger.success('Successfully terminated')


if __name__ == '__main__':
    """This is executed when run from the command line."""

    args = parse_args(sys.argv[1:])

    logger.configure(handlers=[{'sink': sys.stderr, 'level': args.loglevel}])

    merge_trgt_vcfs(args.vcfs, args.bed, args.anonymize_prefix, args.outfile, ' '.join(sys.argv))
