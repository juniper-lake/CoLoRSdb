#!/usr/bin/env python3
"""
Tests to make sure the vcf parser behaves as expected.

Usage: pytest -v test_parser.py
"""

__version__ = "0.1.0"

import pytest
from tempfile import NamedTemporaryFile
from parser import VCFParser
import random
import inspect
import sys


def get_vcf_file(vcf_lines):
    """
    Take an iterator with vcf lines and prints them to a temporary file.

    Arguments:
        vcf_lines (iterator): An iterator with vcf lines

    Returns:
        filename (str): The path to the vcf file
    """
    vcf_file = NamedTemporaryFile(mode="w+t", delete=False, suffix=".vcf")
    vcf_file.writelines(vcf_lines)
    vcf_file.seek(0)
    vcf_file.close()

    return vcf_file.name


def test_add_meta():
    """Test adding metadata line and printing header"""
    vcf_lines = [
        "##fileformat=VCFv4.1\n",
        '##INFO=<ID=MQ,Number=1,Type=Float,Description="RMS Mapping Quality">\n',
        "##contig=<ID=1,length=249250621,assembly=b37>\n",
        "##reference=file:///humgen/gsa-hpprojects/GATK/bundle"
        "/current/b37/human_g1k_v37.fasta\n",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
        "father\tmother\tproband\n",
        "1\t11900\t.\tA\tT,C\t100\tPASS\tMQ=1\tGT:GQ\t0/1:60\t0/2:60\t1/2:60\n",
    ]
    vcf_file = get_vcf_file(vcf_lines)
    vcf = VCFParser(vcf_file)

    test_key, test_value = "commandline", "VCF anonymizer"
    vcf.add_meta(test_key, test_value)
    assert (
        f"##{test_key}={test_value}" in vcf.metadata
    ), "metadata line was not correctly added"

    test_key, test_value = "INFO", "VCF anonymizer"
    with pytest.raises(ValueError):
        vcf.add_meta(test_key, test_value)


def test_change_sample_names():
    """Test adding metadata line and printing header"""
    vcf_lines = [
        "##fileformat=VCFv4.1\n",
        '##INFO=<ID=MQ,Number=1,Type=Float,Description="RMS Mapping Quality">\n',
        "##contig=<ID=1,length=249250621,assembly=b37>\n",
        "##reference=file:///humgen/gsa-hpprojects/GATK/bundle"
        "/current/b37/human_g1k_v37.fasta\n",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
        "father\tmother\tproband\n",
        "1\t11900\t.\tA\tT,C\t100\tPASS\tMQ=1\tGT:GQ\t0/1:60\t0/2:60\t1/2:60\n",
    ]
    vcf_file = get_vcf_file(vcf_lines)
    vcf = VCFParser(vcf_file)

    vcf.change_sample_names(["sample1", "sample2", "sample3"])
    assert "sample1" in vcf.samples, "sample names were not correctly added"

    with pytest.raises(ValueError):
        vcf.change_sample_names(["sample1", "sample2"])

    with pytest.raises(ValueError):
        vcf.change_sample_names(["sample1", "sample 2", "sample3"])


def test_sample_shuffling():
    """Test iterating through variants"""
    vcf_lines = [
        "##fileformat=VCFv4.1\n",
        '##INFO=<ID=MQ,Number=1,Type=Float,Description="RMS Mapping Quality">\n',
        "##contig=<ID=1,length=249250621,assembly=b37>\n",
        "##reference=file:///humgen/gsa-hpprojects/GATK/bundle"
        "/current/b37/human_g1k_v37.fasta\n",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
        "father\tmother\tproband\n",
        "1\t11900\t.\tA\tT,C\t100\tPASS\tMQ=1\tGT:GQ\t0/1:60\t0/2:60\t1/2:60\n",
    ]
    vcf_file = get_vcf_file(vcf_lines)
    vcf = VCFParser(vcf_file)

    random.seed(10)

    for variant in vcf:
        variant.shuffle_samples()
        assert variant.sample_data == [
            "1/2:60",
            "0/1:60",
            "0/2:60",
        ], "samples were not correctly shuffled"


def test_update_genotype():
    """Test updating the genotype based on order of the alt alleles"""
    vcf_lines = [
        "##fileformat=VCFv4.1\n",
        '##INFO=<ID=MQ,Number=1,Type=Float,Description="RMS Mapping Quality">\n',
        "##contig=<ID=1,length=249250621,assembly=b37>\n",
        "##reference=file:///humgen/gsa-hpprojects/GATK/bundle"
        "/current/b37/human_g1k_v37.fasta\n",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
        "father\tmother\tproband\n",
        "1\t11900\t.\tA\tT,C\t100\tPASS\tMQ=1\tGT:GQ\t0/1:60\t0/2:60\t1/2:60\n",
    ]
    vcf_file = get_vcf_file(vcf_lines)
    vcf = VCFParser(vcf_file)

    for variant in vcf:
        variant.update_genotypes(ref="A", alts=["G", "C", "T"])
        assert [sample.split(":")[0] for sample in variant.sample_data] == [
            "0/3",
            "0/2",
            "3/2",
        ], "genotypes were not correctly updated"


def test_fix_ploidy_unsupported_format():
    """Test updating the genotype based on order of the alt alleles"""
    vcf_lines = [
        "##fileformat=VCFv4.1\n",
        '##INFO=<ID=MQ,Number=1,Type=Float,Description="RMS Mapping Quality">\n',
        "##contig=<ID=1,length=249250621,assembly=b37>\n",
        "##reference=file:///humgen/gsa-hpprojects/GATK/bundle"
        "/current/b37/human_g1k_v37.fasta\n",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
        "father\tmother\tproband\n",
        "1\t11900\t.\tA\tT,C\t100\tPASS\tMQ=1\tGT:GQ\t0/1:60\t0/2:60\t1/2:60\n",
    ]
    vcf_file = get_vcf_file(vcf_lines)
    vcf = VCFParser(vcf_file)

    for variant in vcf:
        with pytest.raises(ValueError):
            variant.fix_ploidy(
                sexes=["male", "female", "male"], regions=[("1", 1000, 2000000)]
            )


def test_fix_ploidy():
    """Test updating the genotype based on order of the alt alleles"""
    vcf_lines = [
        "##fileformat=VCFv4.1\n",
        '##INFO=<ID=MQ,Number=1,Type=Float,Description="RMS Mapping Quality">\n',
        "##contig=<ID=1,length=249250621,assembly=b37>\n",
        "##reference=file:///humgen/gsa-hpprojects/GATK/bundle"
        "/current/b37/human_g1k_v37.fasta\n",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
        "father\tmother\tproband\n",
        "1\t11900\t.\tA\tT,C\t100\tPASS\tMQ=1\tGT:AD\t0/1:9,1,3\t0/2:4,4\t1/2:0,0,7\n",
    ]
    vcf_file = get_vcf_file(vcf_lines)
    vcf = VCFParser(vcf_file)

    for variant in vcf:
        variant.fix_ploidy(
            sexes=["male", "female", "male"], regions=[("1", 1000, 2000000)]
        )
        assert variant.sample_data == [
            "0:9,1,3",
            "0/2:4,4",
            "2:0,0,7",
        ], "ploidy was not updated correctly"


def run_all_test_functions(mod):
    """Run all functions that don't require inputs"""
    all_functions = inspect.getmembers(mod, inspect.isfunction)
    for key, value in all_functions:
        if str(inspect.signature(value)) == "()":
            value()


if __name__ == "__main__":
    """This is executed when run from the command line"""

    # this is not exectuted by pytest, it's for development purposes
    run_all_test_functions(sys.modules[__name__])
