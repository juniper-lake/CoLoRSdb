#!/usr/bin/env python3
"""
The VCFParser class was built to facilitate anonymizing VCF files.
"""

__version__ = "0.1.0"

from loguru import logger
import gzip
import os
from codecs import open, getreader
import re
import random


class VCFParser(object):
    """Object representing a VCF file with iterator of variant records."""

    def __init__(self, infile: str):
        super().__init__()
        self.vcf = None
        self.beginning = True
        self.metadata = []
        self.header_columns = [
            "CHROM",
            "POS",
            "ID",
            "REF",
            "ALT",
            "QUAL",
            "FILTER",
            "INFO",
            "FORMAT",
        ]
        self.samples = []
        self.sample_pattern = re.compile("^[a-zA-Z0-9_]+$")
        self.shuffle_samples = False

        logger.info(f"Reading vcf from file {infile}")
        file_name, file_extension = os.path.splitext(infile)
        if file_extension == ".gz":
            logger.debug("VCF is zipped")
            self.vcf = getreader("utf-8")(gzip.open(infile), errors="strict")
        elif file_extension == ".vcf":
            self.vcf = open(infile, mode="r", encoding="utf-8", errors="strict")
        else:
            raise IOError(
                "File is not in a supported format!\n"
                " Or use correct ending(.vcf or .vcf.gz)"
            )

        # Parse the metadata lines
        logger.debug("Reading first line.")
        self.next_line = self.vcf.readline().rstrip()
        # First line is always a metadata line
        if not self.next_line.startswith("#"):
            raise IOError("VCF files allways have to start with a metadata line.")
        while self.next_line.startswith("#"):
            if self.next_line.startswith("##"):
                self.metadata.append(self.next_line)
            elif self.next_line.startswith("#"):
                self.parse_header_line(self.next_line)
            self.next_line = self.vcf.readline().rstrip()

    def parse_header_line(self, line: str):
        """Make sure header line columns are correct and get samples names."""
        self.header = line[1:].rstrip().split("\t")
        if self.header[:9] != self.header_columns:
            raise IOError("VCF header does not contain the correct fields.")
        self.samples = self.header[9:]
        if len(self.samples) < 1:
            raise IOError("This VCF file does not contain any samples.")

    def add_meta(self, key: str, value: str):
        """Add metadata line in form of "key=value" to end of metadata list."""
        if key.lower() not in ["info", "filter", "format", "contig", "fileformat"]:
            if (meta_line := f"##{key}={value}") not in self.metadata:
                self.metadata.append(meta_line)
        else:
            raise ValueError(f"Metadata key {key} is reserved. Please use another key.")

    def change_sample_names(self, sample_names: list):
        """Change sample names and update header line with provided list of names."""
        if len(sample_names) != len(self.samples):
            raise ValueError(
                f"Number of sample names ({len(sample_names)}) does \
                             not match number of samples ({len(self.samples)})"
            )
        elif not all(
            self.sample_pattern.match(sample_name) for sample_name in sample_names
        ):
            raise ValueError(
                "Sample names can only contain letters, numbers \
                             and underscores."
            )
        else:
            logger.info(f"Changing sample names to: {sample_names}")
            self.samples = sample_names
            self.header = self.header[:9] + self.samples

    def print_metadata(self):
        """Returns a string with the metadata and header line."""
        return "\n".join(self.metadata)

    def print_headerline(self):
        """Returns a string with the metadata and header line."""
        return "#" + "\t".join(self.header)

    def __iter__(self):
        """Iterate over the variants in the VCF file."""
        # We need to treat the first case as an exception because we read it in the init
        if self.shuffle_samples:
            logger.info("Sample data is being shuffled.")

        if self.beginning:
            if self.next_line:
                variant_line = self.next_line.split("\t")

                logger.debug("Checking if variant line is malformed")
                if len(self.header) != len(variant_line):
                    raise SyntaxError(
                        f"One of the variant lines is malformed: \
                                      {self.next_line}"
                    )

                variant = VariantRecord(variant_line, self.shuffle_samples)

                yield variant

                self.beginning = False

        for line in self.vcf:
            line = line.rstrip()
            variant_line = line.split("\t")

            if not line.startswith("#"):
                logger.debug("Checking if variant line is malformed")
                if len(self.header) != len(variant_line):
                    raise SyntaxError(f"One of the variant lines is malformed: {line}")

                variant = VariantRecord(variant_line, self.shuffle_samples)

                yield variant

    def __call__(self, shuffle_samples: bool):
        """Return an iterator over variants in VCF, optionally shuffling samples."""
        self.shuffle_samples = shuffle_samples
        return self


class VariantRecord(object):
    """Object representing a single variant record in a VCF file."""

    def __init__(self, variant_line, shuffle_samples=False):
        super().__init__()
        logger.debug("Creating VariantRecords object")
        self.line = variant_line
        self.chrom = variant_line[0]
        self.pos = variant_line[1]
        self.id = variant_line[2]
        self.ref = variant_line[3]
        self.alts = variant_line[4].split(",")
        self.alleles = [self.ref] + self.alts
        self.info = variant_line[7]
        self.format = variant_line[8]
        self.variant_data = variant_line[:9]
        self.sample_data = variant_line[9:]
        self.genotypes = [
            map(int, sample.split(":")[0].split("/")) for sample in self.sample_data
        ]
        self.full_genotypes = [
            [self.alleles[hap] for hap in genotype] for genotype in self.genotypes
        ]

        if shuffle_samples:
            self.sample_data = random.sample(self.sample_data, len(self.sample_data))
            self.line = self.variant_data + self.sample_data

    def __str__(self):
        """Return a string representation of the variant record."""
        return "\t".join(self.line)

    def update_genotypes(self, ref: str, alts: list):
        """Update genotype and haplotypes based on provided alleles."""
        alleles = [ref] + alts
        # get new genotype for each sample
        self.genotypes = [
            [alleles.index(hap) for hap in full_genotype]
            for full_genotype in self.full_genotypes
        ]
        # use new genotype to update sample data string for each sample
        self.sample_data = [
            ":".join(
                ["/".join(map(str, self.genotypes[i]))]
                + self.sample_data[i].split(":")[1:]
            )
            for i in range(len(self.genotypes))
        ]
        # update full variant line with new sample data
        self.line = self.variant_data + self.sample_data
