#!/usr/bin/env python3
"""Check a BAM file for the presence of specific tags and infer the history of the BAM file."""

__version__ = '0.1.0'

import collections
import json
import math
import statistics
import sys
from argparse import ArgumentParser

import pysam


def get_minimizer(sequence):
    """Determine minimizer of homopolymer-compressed sequence."""
    rev_comp = sequence[::-1].translate(str.maketrans('ATCG', 'TAGC'))
    return min(sequence, rev_comp)


def rq_from_bq(bq):
    """Compute read quality from an array of base qualities."""
    # Compute read quality from an array of base qualities; cap at Q60
    readLen = len(bq)
    expectedErrors = sum([math.pow(10, -0.1 * x) for x in bq])
    return min(60, math.floor(-10 * math.log10(expectedErrors / readLen)))


def error_rate_from_rqv(rqv):
    """Compute error rate from read quality value."""
    # Compute error rate from read quality value
    return math.pow(10, -0.1 * rqv)


def rqv_from_error_rate(error_rate):
    """Compute read quality value from error rate."""
    # Compute read quality value from error rate
    return -10 * math.log10(error_rate)


def check_bam_file(bam_file_path, n_records):
    """Check a BAM file for the presence of specific tags and infer the history of the BAM file."""
    # Load the BAM file
    save = pysam.set_verbosity(0)  # suppress [E::idx_find_and_load]
    bam_file = pysam.AlignmentFile(bam_file_path, 'rb', check_sq=False)
    pysam.set_verbosity(save)  # restore warnings

    # Check if the BAM file is aligned
    aligned = bool(bam_file.nreferences)

    # Initialize a flag to check if the BAM file contains non-CCS reads
    ccs = True
    hifi = True

    # Initialize sets to collect tag names, bq scores, and rq scores
    unique_tags = set()
    unique_bq_scores = set()
    readqv_list = list()
    readlength_list = list()
    first_last_bases = collections.Counter()

    # Initialize the output dictionary
    output = dict()

    # Read the first n records from the BAM file
    for i, record in enumerate(bam_file):
        if i >= n_records:
            break

        if 'ccs' not in record.query_name:
            ccs = False

        # Add the read quality to the list of read qualities
        if record.has_tag('rq'):
            errorrate = 1.0 - record.get_tag('rq')
            readqv = 60 if errorrate == 0 else math.floor(-10 * math.log10(errorrate))
        else:
            readqv = rq_from_bq(record.query_qualities)
        readqv_list.append(readqv)

        if readqv < 20:
            hifi = False

        first_last_bases[get_minimizer(record.query_sequence[:16])] += 1
        first_last_bases[get_minimizer(record.query_sequence[-16:])] += 1

        # Add the read length to the list of read lengths
        readlength_list.append(record.query_length)
        # Add the tag names to the set of unique tags
        unique_tags.update(tag[0] for tag in record.tags)
        # Add the base quality scores to the set of unique scores
        unique_bq_scores.update(record.query_qualities)

    # Check for the presence of specific tags and infer the history of the BAM file
    kinetics_tags = {'fi', 'ri', 'fp', 'rp', 'ip', 'pw'}
    segment_reads_tags = {'di', 'dl', 'dr', 'ds'}
    demultiplexed_tags = {'bc', 'bq', 'bl', 'bt', 'ql', 'qt', 'bx', 'ls'}
    base_modification_tags = {'MM', 'ML', 'Mm', 'Ml'}
    haplotag_tags = {'HP', 'PS', 'PC'}

    output['file'] = bam_file_path
    output['kinetics'] = bool(unique_tags & kinetics_tags)
    output['segment_reads'] = bool(unique_tags & segment_reads_tags)
    output['demultiplexed'] = bool(unique_tags & demultiplexed_tags)
    output['base_modification'] = bool(unique_tags & base_modification_tags)
    output['aligned'] = aligned
    output['haplotagged'] = bool(unique_tags & haplotag_tags)

    read_groups = bam_file.header.get('RG', [])
    output['read_groups'] = [';'.join([':'.join([key, val]) for key, val in rg.items()]) for rg in read_groups]
    output['instrument_models'] = list({rg.get('PM', 'N/A') for rg in read_groups})
    output['libraries'] = list({rg.get('LB', 'N/A') for rg in read_groups})
    output['samples'] = list({rg.get('SM', 'N/A') for rg in read_groups})
    output['movies'] = list({rg.get('PU', 'N/A') for rg in read_groups})

    output['ccs'] = ccs
    output['hifi'] = hifi
    error_rates = [error_rate_from_rqv(rqv) for rqv in readqv_list]
    output['read_quality'] = {
        'mean': round(rqv_from_error_rate(statistics.mean(error_rates)), 2),
        'median': round(rqv_from_error_rate(statistics.median(error_rates)), 2),
        'min': round(rqv_from_error_rate(max(error_rates)), 2),
        'max': round(rqv_from_error_rate(min(error_rates)), 2),
    }
    output['read_length'] = {
        'mean': round(statistics.mean(readlength_list), 2),
        'median': round(statistics.median(readlength_list)),
        'min': round(min(readlength_list)),
        'max': round(max(readlength_list)),
    }
    output['bq_bins'] = len(unique_bq_scores)
    output['potentially_multiplexed'] = first_last_bases.most_common(1)[0][1] / n_records > 0.05

    # Close the BAM file
    bam_file.close()
    return output


def main(arguments):
    """Parse the command-line arguments and check the BAM file."""
    output = check_bam_file(args.bam_file_path, args.n_records)
    file = sys.stdout
    print(json.dumps(output, indent=2), file=file)


if __name__ == '__main__':
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('bam_file_path', help='Path to the BAM file.')
    parser.add_argument(
        '--n_records',
        type=int,
        default=10000,
        help='Number of records to read from the BAM file.  Default: 10000.',
    )
    parser.add_argument('--version', action='version', version=__version__)
    args = parser.parse_args()
    main(args)
