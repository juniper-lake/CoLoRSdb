#!/usr/bin/env python3
"""
Tests to make sure the relatedness screener behaves as expected.

Usage: pytest -v test_relatedness_screener.py
"""

__version__ = "0.1.0"

from tempfile import NamedTemporaryFile
import pytest
from relatedness_screener import NoRelation, most_common, flag_related_samples
import inspect
import sys


def get_file(lines):
    """
    Take an iterator with vcf lines and prints them to a temporary file.
    
    Arguments:
        vcf_lines (iterator): An iterator with vcf lines
    
    Returns:
        filename (str): The path to the vcf file
    """
    file = NamedTemporaryFile(mode='w+t', delete=False, suffix='.vcf')
    file.writelines(lines)
    file.seek(0)
    file.close()
    
    return file.name


def test_wrong_header():
    """Test adding metadata line and printing header"""
    lines = [
        '#sample_a\tsample_\trelatedness\tibs0\tibs2\thom_concordance\thets_a\thets_b\thets_ab\t'
        'shared_hets\thom_alts_a\thom_alts_b\tshared_hom_alts\tn\tx_ibs0\tx_ibs2\texpected_relatedness\n',
        'HG002\tHG003\t0.982\t0\t982\t0.977\t293\t280\t336\t165\t480\t475\t464\t988\t0\t85\t-1.0',
        ]
    file = get_file(lines)

    with pytest.raises(IOError):
        NoRelation(file)


def test_relatedness_range():
    """Test that max_relatedness is within range of 0 and 1"""
    lines = [
        '#sample_a\tsample_b\trelatedness\tibs0\tibs2\thom_concordance\thets_a\thets_b\thets_ab\t'
        'shared_hets\thom_alts_a\thom_alts_b\tshared_hom_alts\tn\tx_ibs0\tx_ibs2\texpected_relatedness\n',
        'HG002\tHG003\t0.982\t0\t982\t0.977\t293\t280\t336\t165\t480\t475\t464\t988\t0\t85\t-1.0\n',
        ]
    file = get_file(lines)

    with pytest.raises(ValueError):
        NoRelation(file, max_relatedness=2.0)


def test_most_common():
    """Test that the most common function works as expected"""
    assert most_common([1, 1, 2, 2, 3,3]) == ([1,2,3],2), "all three values occur twice"
    assert most_common(['child', 'child', 'child', 'dad', 'dad', 'mom']) == (['child'],3), "child would get dropped"


def test_iteration():
    """Test that more than one sample is removed"""
    lines = [
        '#sample_a\tsample_b\trelatedness\tibs0\tibs2\thom_concordance\thets_a\thets_b\thets_ab\t'
        'shared_hets\thom_alts_a\thom_alts_b\tshared_hom_alts\tn\tx_ibs0\tx_ibs2\texpected_relatedness\n',
        'child\tdad\t0.475\t0\t982\t0.977\t293\t280\t336\t165\t480\t475\t464\t988\t0\t85\t-1.0\n',
        'child\tmom\t0.512\t0\t982\t0.977\t293\t280\t336\t165\t480\t475\t464\t988\t0\t85\t-1.0\n',
        'child\tbro\t0.497\t0\t982\t0.977\t293\t280\t336\t165\t480\t475\t464\t988\t0\t85\t-1.0\n',
        'bro\tdad\t0.475\t0\t982\t0.977\t293\t280\t336\t165\t480\t475\t464\t988\t0\t85\t-1.0\n',
        'bro\tmom\t0.503\t0\t982\t0.977\t293\t280\t336\t165\t480\t475\t464\t988\t0\t85\t-1.0\n',
        'mom\tdad\t0.05\t0\t982\t0.977\t293\t280\t336\t165\t480\t475\t464\t988\t0\t85\t-1.0\n'
        ]
    file = get_file(lines)

    # with pytest.raises(ValueError):
    cohort = NoRelation(file)
    assert cohort.dropped_samples == ['child','bro'], "only parents should be kept"


def test_sort_by_coverage():
    """Test dropped samples with lowest coverage if there's a tie in number of relations"""
    lines = [
        '#sample_a\tsample_b\trelatedness\tibs0\tibs2\thom_concordance\thets_a\thets_b\thets_ab\t'
        'shared_hets\thom_alts_a\thom_alts_b\tshared_hom_alts\tn\tx_ibs0\tx_ibs2\texpected_relatedness\n',
        'child\tdad\t0.475\t0\t982\t0.977\t293\t280\t336\t165\t480\t475\t464\t988\t0\t85\t-1.0\n'
        ]
    file = get_file(lines)

    cohort = NoRelation(file,sample_order=['child','dad'], coverages=[1,30])
    assert cohort.dropped_samples == ['child'], "child has lower coverage than dad"

    cohort = NoRelation(file,sample_order=['child','dad'], coverages=[30,1])
    assert cohort.dropped_samples == ['dad'], "dad has lower coverage than child"

    cohort = NoRelation(file,sample_order=['child','dad'])
    assert cohort.dropped_samples == ['child'], "child is first in file"


def test_wrong_sample_names():
    """Test if sample names set on CLI are in the file"""
    lines = [
        '#sample_a\tsample_b\trelatedness\tibs0\tibs2\thom_concordance\thets_a\thets_b\thets_ab\t'
        'shared_hets\thom_alts_a\thom_alts_b\tshared_hom_alts\tn\tx_ibs0\tx_ibs2\texpected_relatedness\n',
        'child\tdad\t0.475\t0\t982\t0.977\t293\t280\t336\t165\t480\t475\t464\t988\t0\t85\t-1.0\n',
        'child\tmom\t0.512\t0\t982\t0.977\t293\t280\t336\t165\t480\t475\t464\t988\t0\t85\t-1.0\n',
        'mom\tdad\t0.05\t0\t982\t0.977\t293\t280\t336\t165\t480\t475\t464\t988\t0\t85\t-1.0\n'
        ]
    file = get_file(lines)

    with pytest.raises(ValueError):
        NoRelation(file, sample_order=['not_child','dad','mom'])


def test_wrong_coverage_type():
    """Test if coverages are integers or not"""
    lines = [
        '#sample_a\tsample_b\trelatedness\tibs0\tibs2\thom_concordance\thets_a\thets_b\thets_ab\t'
        'shared_hets\thom_alts_a\thom_alts_b\tshared_hom_alts\tn\tx_ibs0\tx_ibs2\texpected_relatedness\n',
        'child\tdad\t0.475\t0\t982\t0.977\t293\t280\t336\t165\t480\t475\t464\t988\t0\t85\t-1.0\n',
        'child\tmom\t0.512\t0\t982\t0.977\t293\t280\t336\t165\t480\t475\t464\t988\t0\t85\t-1.0\n',
        'mom\tdad\t0.05\t0\t982\t0.977\t293\t280\t336\t165\t480\t475\t464\t988\t0\t85\t-1.0\n'
        ]
    file = get_file(lines)

    with pytest.raises(ValueError):
        NoRelation(file, sample_order=['child','dad','mom'], coverages=['10X','20X','30X'])


def test_sample_order():
    """Test sort order of output since this affects WDL functionality"""
    lines = [
        '#sample_a\tsample_b\trelatedness\tibs0\tibs2\thom_concordance\thets_a\thets_b\thets_ab\t'
        'shared_hets\thom_alts_a\thom_alts_b\tshared_hom_alts\tn\tx_ibs0\tx_ibs2\texpected_relatedness\n',
        'child\tdad\t0.475\t0\t982\t0.977\t293\t280\t336\t165\t480\t475\t464\t988\t0\t85\t-1.0\n',
        'child\tmom\t0.512\t0\t982\t0.977\t293\t280\t336\t165\t480\t475\t464\t988\t0\t85\t-1.0\n',
        'child\tbro\t0.497\t0\t982\t0.977\t293\t280\t336\t165\t480\t475\t464\t988\t0\t85\t-1.0\n',
        'bro\tdad\t0.475\t0\t982\t0.977\t293\t280\t336\t165\t480\t475\t464\t988\t0\t85\t-1.0\n',
        'bro\tmom\t0.503\t0\t982\t0.977\t293\t280\t336\t165\t480\t475\t464\t988\t0\t85\t-1.0\n',
        'mom\tdad\t0.05\t0\t982\t0.977\t293\t280\t336\t165\t480\t475\t464\t988\t0\t85\t-1.0\n'
        ]
  
    file = get_file(lines)

    cohort = NoRelation(file)
    assert cohort.samples == ['bro', 'child', 'dad', 'mom'], "sample name order determined from pairs file"
    assert cohort.sorted_keep_drop == ['drop', 'drop', 'keep', 'keep'], "sample name order determined from pairs file"

    cohort = NoRelation(file, sample_order=['child','dad','mom','bro'])
    assert cohort.samples == ['child','dad','mom','bro'], "sample order explicitly set"
    assert cohort.sorted_keep_drop == ['drop', 'keep', 'keep', 'drop'], "sample order explicitly set"
    

def test_flag_related_samples(capfd):
    """Test output of main function to flag related samples"""
    lines = [
        '#sample_a\tsample_b\trelatedness\tibs0\tibs2\thom_concordance\thets_a\thets_b\thets_ab\t'
        'shared_hets\thom_alts_a\thom_alts_b\tshared_hom_alts\tn\tx_ibs0\tx_ibs2\texpected_relatedness\n',
        'child\tdad\t0.475\t0\t982\t0.977\t293\t280\t336\t165\t480\t475\t464\t988\t0\t85\t-1.0\n',
        'child\tmom\t0.512\t0\t982\t0.977\t293\t280\t336\t165\t480\t475\t464\t988\t0\t85\t-1.0\n',
        'child\tbro\t0.497\t0\t982\t0.977\t293\t280\t336\t165\t480\t475\t464\t988\t0\t85\t-1.0\n',
        'bro\tdad\t0.475\t0\t982\t0.977\t293\t280\t336\t165\t480\t475\t464\t988\t0\t85\t-1.0\n',
        'bro\tmom\t0.503\t0\t982\t0.977\t293\t280\t336\t165\t480\t475\t464\t988\t0\t85\t-1.0\n',
        'mom\tdad\t0.05\t0\t982\t0.977\t293\t280\t336\t165\t480\t475\t464\t988\t0\t85\t-1.0\n'
        'dad_bro\tchild\t0.26\t0\t982\t0.977\t293\t280\t336\t165\t480\t475\t464\t988\t0\t85\t-1.0\n'
        'dad_bro\tbro\t0.24\t0\t982\t0.977\t293\t280\t336\t165\t480\t475\t464\t988\t0\t85\t-1.0\n'
        'dad_bro\tdad\t0.50\t0\t982\t0.977\t293\t280\t336\t165\t480\t475\t464\t988\t0\t85\t-1.0\n'
        'dad_bro\tmom\t0.03\t0\t982\t0.977\t293\t280\t336\t165\t480\t475\t464\t988\t0\t85\t-1.0\n'
        ]
  
    file = get_file(lines)
    flag_related_samples(file, max_relatedness=0.125, sample_order=['dad','dad_bro', 'mom','bro','child'], coverages=[30, 40, 20, 50, 60])
    out, err = capfd.readouterr()
    assert out == "dad\tdad_bro\tmom\tbro\tchild\ndrop\tkeep\tkeep\tdrop\tdrop\n", "related samples flagged correctly"


def run_all_test_functions(mod):
    """Run all functions that don't require inputs"""
    all_functions = inspect.getmembers(mod, inspect.isfunction)
    for key, value in all_functions:
        if str(inspect.signature(value)) == "()":
            value()


if __name__ == '__main__':
    """ This is executed when run from the command line """

    run_all_test_functions(sys.modules[__name__])
