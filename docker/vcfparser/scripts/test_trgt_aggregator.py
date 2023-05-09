#!/usr/bin/env python3
"""
Tests to make sure trgt aggregation behaves as expected.

Usage: pytest -v test_trgt_aggregator.py
"""

__version__ = "0.1.0"

import pytest
from tempfile import NamedTemporaryFile
from parser import VCFParser
import random
import inspect
import sys
from aggregate_trgt_vcfs import aggregate_trgt_vcfs


def get_vcf_file(vcf_lines):
    """
    Take an iterator with vcf lines and prints them to a temporary file.
    
    Arguments:
        vcf_lines (iterator): An iterator with vcf lines
    
    Returns:
        filename (str): The path to the vcf file
    """
    vcf_file = NamedTemporaryFile(mode='w+t', delete=False, suffix='.vcf')
    vcf_file.writelines(vcf_lines)
    vcf_file.seek(0)
    vcf_file.close()
    
    return vcf_file.name


def test_trgt_aggregation():
    vcf_lines1 = [
        '##fileformat=VCFv4.2\n',
        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1\n',
        'chr2\t126162\t.\tREF\tALT1\t0\t.\tTRID=chr2_126161_126197;END=126197;MOTIFS=CTCTCC;STRUC=(CTCTCC)n\tGT:AL:ALLR:SD:MC:MS:AP:AM\t0/1:36,33:35-38,31-33:18,14:6,6:0(0-36),0(0-33):1,1:.,.\n'
    ]

    vcf_lines2 = [
        '##fileformat=VCFv4.2\n',
        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample2\n',
        'chr2\t126162\t.\tREF\tALT2,ALT3\t0\t.\tTRID=chr2_126161_126197;END=126197;MOTIFS=CTCTCC;STRUC=(CTCTCC)n\tGT:AL:ALLR:SD:MC:MS:AP:AM\t1/2:36,33:35-38,31-33:18,14:6,6:0(0-36),0(0-33):1,1:.,.\n'
    ]

    vcf_lines3 = [
        '##fileformat=VCFv4.2\n',
        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample3\n',
        'chr2\t126162\t.\tREF\t.\t0\t.\tTRID=chr2_126161_126197;END=126197;MOTIFS=CTCTCC;STRUC=(CTCTCC)n\tGT:AL:ALLR:SD:MC:MS:AP:AM\t0/0:36,33:35-38,31-33:18,14:6,6:0(0-36),0(0-33):1,1:.,.\n'
    ]

    vcf_lines4 = [
        '##fileformat=VCFv4.2\n',
        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample4\n',
        'chr2\t126162\t.\tREF\t.\t0\t.\tTRID=chr2_126161_126197;END=126197;MOTIFS=CTCTCC;STRUC=(CTCTCC)n\tGT:AL:ALLR:SD:MC:MS:AP:AM\t0/0:36,33:35-38,31-33:18,14:6,6:0(0-36),0(0-33):1,1:.,.\n'
        'chr3\t126162\t.\tREF\t.\t0\t.\tTRID=chr2_126161_126197;END=126197;MOTIFS=CTCTCC;STRUC=(CTCTCC)n\tGT:AL:ALLR:SD:MC:MS:AP:AM\t0/0:36,33:35-38,31-33:18,14:6,6:0(0-36),0(0-33):1,1:.,.\n'
    ]


    vcf_file1 = get_vcf_file(vcf_lines1)
    vcf_file2 = get_vcf_file(vcf_lines2)
    vcf_file3 = get_vcf_file(vcf_lines3)
    vcf_file4 = get_vcf_file(vcf_lines4)

    # good check output
    aggregate_trgt_vcfs([vcf_file1, vcf_file2, vcf_file3])

    # first vcf has fewer variants than at least one of the others
    with pytest.raises(IndexError):
        aggregate_trgt_vcfs([vcf_file1, vcf_file2, vcf_file3, vcf_file4])
    # first vcf has more variants than at least one of the others
    with pytest.raises(IndexError):
        aggregate_trgt_vcfs([vcf_file4,vcf_file1, vcf_file2, vcf_file3])


def run_all_test_functions(mod):
    """Run all functions that don't require inputs"""
    all_functions = inspect.getmembers(mod, inspect.isfunction)
    for key, value in all_functions:
        if str(inspect.signature(value)) == "()":
            value()



if __name__ == '__main__':
    """ This is executed when run from the command line """

    run_all_test_functions(sys.modules[__name__])
