#!/usr/bin/env python3
"""
This tool anonymizes multi-sample vcf files by:
(1) randomizing the order of sample-specific data on a per-variant basis and
(2) changing the sample names to an arbitrary counter starting with the provided prefix.

It also adds a metadata line to the vcf file to indicate that the file has been anonymized.
"""

__version__ = "0.1.0"

from loguru import logger
import sys
import argparse
import gzip
import os
from codecs import open, getreader
import re
import random

# reference: https://github.com/moonso/vcf_parser/blob/develop/vcf_parser/parser.py
class VCFParser(object):
    """Object representing a VCF file including iterator that returns VariantRecord objects."""
    def __init__(self, infile:str):
        super().__init__()
        self.vcf = None
        self.beginning = True
        self.metadata = []
        self.header_columns=['CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT']
        self.samples = []
        self.sample_pattern = re.compile('^[a-zA-Z0-9_]+$')
        self.shuffle_samples = False

        logger.info(f"Reading vcf from file {infile}")
        file_name, file_extension = os.path.splitext(infile)
        if file_extension == '.gz':
            logger.debug("VCF is zipped")
            self.vcf = getreader('utf-8')(gzip.open(infile), errors='strict')
        elif file_extension == '.vcf':
            self.vcf = open(infile, mode='r', encoding='utf-8', errors='strict')
        else:
            raise IOError("File is not in a supported format!\n"
                                " Or use correct ending(.vcf or .vcf.gz)")
        
        # Parse the metadata lines
        logger.debug("Reading first line.")
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


    def parse_header_line(self, line:str):
        """Make sure header line matches expectations and get samples names from header line."""
        self.header = line[1:].rstrip().split('\t')
        if self.header[:9] != self.header_columns:
            raise IOError("VCF header does not contain the correct fields.")
        self.samples = self.header[9:]
        if len(self.samples) < 1:
            raise IOError("This VCF file does not contain any samples.")
        

    def add_meta(self, key:str, value:str):
        """Add metadata line in form of "key=value" to end of metadata list."""
        if key.lower() not in ['info','filter','format','contig','fileformat']:
            if (meta_line := f'##{key}={value}') not in self.metadata:
                self.metadata.append(meta_line)
        else:
            raise ValueError(f"Metadata key {key} is reserved. Please use another key.")
    

    def change_sample_names(self, sample_names:list):
        """Change sample names and update header line with the provided list of names."""
        if len(sample_names) != len(self.samples):
            raise ValueError(f"Number of sample names ({len(sample_names)}) does not match number of samples ({len(self.samples)})")
        elif not all(self.sample_pattern.match(sample_name) for sample_name in sample_names):
            raise ValueError(f"Sample names can only contain letters, numbers and underscores.")
        else:
            logger.info(f"Changing sample names to: {sample_names}")
            self.samples = sample_names
            self.header = self.header[:9] + self.samples


    
    def print_header(self):
        """Returns a string with the metadata and header line."""
        return '\n'.join(self.metadata) + '\n' + '#' + '\t'.join(self.header) + '\n'
    

    def __iter__(self):
        """Iterate over the variants in the VCF file."""
        # We need to treat the first case as an exception because we read it in the init
        if self.shuffle_samples:
            logger.info(f"Sample data is being shuffled.")
            
        if self.beginning:
            if self.next_line:
                variant_line = self.next_line.split('\t')

                logger.debug("Checking if variant line is malformed")
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
                logger.debug("Checking if variant line is malformed")
                if len(self.header) != len(variant_line):
                    raise SyntaxError("One of the variant lines is malformed: {0}".format(
                        line
                    ))
                
                variant = VariantRecord(variant_line, self.shuffle_samples)

                yield variant


    def __call__(self, shuffle_samples:bool):
        """Return an iterator over the variants in the VCF file, allowing for optional shuffling of samples."""
        self.shuffle_samples = shuffle_samples
        return self


class VariantRecord(object):
    """Object representing a single variant record in a VCF file."""
    def __init__(self, variant_line, shuffle_samples=False):
        super().__init__()
        logger.debug("Creating VariantRecords object")

        self.variant_data = variant_line[:9]
        self.original_sample_data = variant_line[9:]                
        self.sample_data = self.original_sample_data.copy()
        
        if shuffle_samples:
            self.sample_data = random.sample(self.sample_data, len(self.sample_data))


    def __str__(self):
        """Return a string representation of the variant record."""
        return '\t'.join(self.variant_data + self.sample_data)
    

def cli(args):
    """Parse VCF, change sample names, update metadata, and print shuffled VCF to stdout or file."""
    
    if args.outfile:
        if os.path.isfile(args.outfile):
            raise IOError(f"Output file {args.outfile} already exists.")
        
    vcf = VCFParser(args.vcf)

    # change sample names
    if (n_samples := len(vcf.samples)) < 2:
       raise IOError("This VCF file only contains one sample, so anonymization is not useful.")
    else:
        vcf.change_sample_names([f'{args.prefix}_{i}' for i in range(1, n_samples+1)])

    # add metadata line
    vcf.add_meta('commandline', 'vcf_anonymizer.py was used to anonymize this file by randomizing the order of \
                 sample-specific data (i.e. genotypes) on a per-variant basis and changing the sample names')

    if args.outfile:
        logger.info(f"Writing anonymized VCF to {args.outfile}")
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

    logger.success("Successfully terminated")
    
    
if __name__ == '__main__':
    """This is executed when run from the command line."""
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

    logger.configure(handlers=[{"sink": sys.stderr, "level": args.loglevel}])

    cli(args)
