#!/usr/bin/env python3
"""
This tool anonymizes multi-sample vcf files by:
(1) randomizing the order of sample-specific data on a per-variant basis and
(2) changing the sample names to an arbitrary counter starting with the provided prefix.

Reference: https://github.com/moonso/vcf_parser/blob/develop/vcf_parser/parser.py
"""

__version__ = "0.1.0"

import argparse
import gzip
import logging
import os
from codecs import open, getreader
import re
import random


class VCFParser(object):
    """docstring for VCFParser"""
    def __init__(self, infile, source=None):
        super().__init__()
        self.logger = logging.getLogger(__name__)
        self.vcf = None
        self.beginning = True
        self.metadata = []
        self.header_columns=['CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT']
        self.samples = []
        self.sample_pattern = re.compile('^[a-zA-Z0-9_]+$')
        self.shuffle_samples = False

        self.logger.info(f"Reading vcf from file {infile}")
        file_name, file_extension = os.path.splitext(infile)
        if file_extension == '.gz':
            self.logger.debug("VCF is zipped")
            self.vcf = getreader('utf-8')(gzip.open(infile), errors='strict')
        elif file_extension == '.vcf':
            self.vcf = open(infile, mode='r', encoding='utf-8', errors='strict')
        else:
            raise IOError("File is not in a supported format!\n"
                                " Or use correct ending(.vcf or .vcf.gz)")
        
        # Parse the metadata lines
        self.logger.debug("Reading first line.")
        self.next_line = self.vcf.readline().rstrip()
        # First line is always a metadata line
        if not self.next_line.startswith('#'):
            raise IOError("VCF files allways have to start with a metadata line.")
        while self.next_line.startswith('#'):
            if self.next_line.startswith('##'):
                self.metadata.append(self.next_line)
            elif self.next_line.startswith('#'):
                self.parse_header_line(self.next_line)
            self.next_line = self.vcf.readline().rstrip()


    def parse_header_line(self, line):
        """Get samples names from header line"""
        self.header = line[1:].rstrip().split('\t')
        if self.header[:9] != self.header_columns:
            raise IOError("VCF header does not contain the correct fields.")
        self.samples = self.header[9:]
        if len(self.samples) < 1:
            raise IOError("This VCF file does not contain any samples.")
        

    def add_meta(self, key:str, value:str):
        """Add metadata line to the metadata list"""
        if key.lower() not in ['info','filter','format','contig','fileformat']:
            if (meta_line := f'##{key}={value}') not in self.metadata:
                self.metadata.append(meta_line)
        else:
            raise ValueError(f"Metadata key {key} is reserved. Please use another key.")
    

    def change_sample_names(self, sample_names:list):
        """Change sample names"""
        if len(sample_names) != len(self.samples):
            raise ValueError(f"Number of sample names ({len(sample_names)}) does not match number of samples ({len(self.samples)})")
        elif not all(self.sample_pattern.match(sample_name) for sample_name in sample_names):
            raise ValueError(f"Sample names can only contain letters, numbers and underscores.")
        else:
            self.logger.info(f"Changing sample names to: {sample_names}")
            self.samples = sample_names
            self.header = self.header[:9] + self.samples


    
    def print_header(self):
        """Print the metadata and header lines"""
        return '\n'.join(self.metadata) + '\n' + '#' + '\t'.join(self.header) + '\n'
    

    def __iter__(self):
        # We need to treat the first case as an exception because we read it in the init
        if self.shuffle_samples:
            self.logger.info(f"Sample data has been shuffled")
            
        if self.beginning:
            if self.next_line:
                variant_line = self.next_line.split('\t')

                self.logger.debug("Checking if variant line is malformed")
                if len(self.header) != len(variant_line):
                    raise SyntaxError("One of the variant lines is malformed: {0}".format(
                        line
                    ))
                
                variant = VariantRecord(variant_line, self.shuffle_samples)

                yield variant
                
                self.beginning = False

        for line in self.vcf:
            line = line.rstrip()
            variant_line = line.split('\t')
 
            if not line.startswith('#'):
                self.logger.debug("Checking if variant line is malformed")
                if len(self.header) != len(variant_line):
                    raise SyntaxError("One of the variant lines is malformed: {0}".format(
                        line
                    ))
                
                variant = VariantRecord(variant_line, self.shuffle_samples)

                yield variant

    def __call__(self, shuffle_samples):
        self.shuffle_samples = shuffle_samples
        return self


class VariantRecord(object):
    """docstring for VariantRecords"""
    def __init__(self, variant_line, shuffle_samples=False):
        super().__init__()
        self.logger = logging.getLogger(__name__)
        self.logger.debug("Creating VariantRecords object")

        self.variant_data = variant_line[:9]
        self.original_sample_data = variant_line[9:]                
        self.sample_data = self.original_sample_data.copy()
        
        if shuffle_samples:
            self.sample_data = random.sample(self.sample_data, len(self.sample_data))


    def __str__(self):
        return '\t'.join(self.variant_data + self.sample_data)
    

def cli(args):
    """   docstring for cli """
    log_format = "%(asctime)s::%(levelname)s::%(name)s::"\
                "%(filename)s::%(lineno)d::%(message)s"
    logging.basicConfig(level=args.loglevel, format=log_format)
    
    if args.outfile:
        if os.path.isfile(args.outfile):
            raise IOError(f"Output file {args.outfile} already exists.")
        
    vcf = VCFParser(args.vcf)

    # change sample names
    if (n_samples := len(vcf.samples)) < 2:
       raise IOError("This VCF file only contains one sample, so anonymization is not useful.")
    else:
        vcf.change_sample_names([f'{args.prefix}_{i}' for i in range(1, n_samples+1)])

    if args.outfile:
        with open(args.outfile, mode='w', encoding='utf-8', errors='strict') as out:
            # print metadata and header line
            out.write(vcf.print_header())

            # print each record with samples shuffled
            for variant in vcf(shuffle_samples=True):
                out.write(str(variant) + '\n')
    else:
        # print header
        print(vcf.print_header())

        # print each record with samples shuffled
        for variant in vcf(shuffle_samples=True):
            print(variant)

    
if __name__ == '__main__':
    """ This is executed when run from the command line """
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("vcf", help="optionally gzipped VCF file to be randomized")
    parser.add_argument("prefix", help="prefix of output VCF file")
    parser.add_argument(
        "--outfile",
        dest="outfile",
        help="output file name"
    )
    parser.add_argument(
        "--loglevel", 
        dest="loglevel", 
        choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'], 
        default='INFO', 
        help="set the logging level"
    )
    parser.add_argument(
        "--version",
        action="version",
        version="%(prog)s (version {version})".format(version=__version__)
    )

    args = parser.parse_args()

    cli(args)
