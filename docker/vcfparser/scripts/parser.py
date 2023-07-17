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
import math


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

        if self.beginning:
            if self.next_line:
                variant_line = self.next_line.split("\t")

                logger.debug("Checking if variant line is malformed")
                if len(self.header) != len(variant_line):
                    raise SyntaxError(
                        f"One of the variant lines is malformed: \
                                      {self.next_line}"
                    )

                variant = VariantRecord(variant_line)

                yield variant

                self.beginning = False

        for line in self.vcf:
            line = line.rstrip()
            variant_line = line.split("\t")

            if not line.startswith("#"):
                logger.debug("Checking if variant line is malformed")
                if len(self.header) != len(variant_line):
                    raise SyntaxError(f"One of the variant lines is malformed: {line}")

                variant = VariantRecord(variant_line)

                yield variant


class VariantRecord(object):
    """Object representing a single variant record in a VCF file."""

    def __init__(self, variant_line):
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
        self.original_sample_data = variant_line[9:]
        self.sample_data = variant_line[9:]

    def __str__(self):
        """Return a string representation of the variant record."""

        return "\t".join(self.line)

    def update_genotypes(self, ref: str, alts: list):
        """Update genotype and haplotypes based on provided alleles."""
        old_alleles = self.alleles + ["."]
        new_alleles = [ref] + alts + ["."]
        unknown_index = len(old_alleles) - 1
        
        # get new genotype as number for each sample, replace "." with unknown index
        genotypes = [
            map(lambda s: int(s.replace(".", str(unknown_index))), 
                sample.split(":")[0].split("/")) 
            for sample in self.sample_data
        ]

        full_genotypes = [
            [old_alleles[hap] for hap in genotype] for genotype in genotypes
        ]
        new_genotypes = [
            [new_alleles.index(hap) if (hap != ".") else hap for hap in full_genotype]
            for full_genotype in full_genotypes
        ]

        # use new genotype to update sample data string for each sample
        self.sample_data = [
            ":".join(
                ["/".join(map(str, new_genotypes[i]))]
                + self.sample_data[i].split(":")[1:]
            )
            for i in range(len(genotypes))
        ]
        # update full variant line with new sample data
        self.line = self.variant_data + self.sample_data

    def shuffle_samples(self):
        """Shuffle sample data."""
        self.sample_data = random.sample(self.sample_data, len(self.sample_data))
        self.line = self.variant_data + self.sample_data

    def check_ploidy(self, chrom: str, pos: int, 
                        regions: list[tuple[str, int, int, str, int]]
                        ):
        """Check if variant is in any of the provided regions."""
        male = ["male", "m", "xy"]
        female = ["female", "f", "xx"]
        valid_sex = male + female
        ploidy_dict = {}
        for region in regions:
            (region_chrom, region_start, region_end, sex, ploidy) = region
            sex = sex.lower()
            if chrom == region_chrom and (region_start < pos <= region_end):
                if sex not in valid_sex:
                    raise ValueError(f"Sex can only be: {valid_sex}")
                elif sex in male:
                    target_sex = "male"
                else:
                    target_sex = "female"
                if ploidy not in [0,1]:
                    raise ValueError("Ploidy in non-diploid regions should be less \
                                     than 0 or 1.")
                ploidy_dict |= {target_sex: ploidy}
        return ploidy_dict
    
    def convert_to_missing(self, sample_idx: int):
        """Convert diploid genotype to missing."""
        format = self.format.split(":")
        sample = self.sample_data[sample_idx].split(":")
        sample_dict = dict(zip(format, sample))
        sample_dict |= {"GT": "./."}
        new_sample = [sample_dict[key] for key in format]
        return new_sample

    def convert_to_haploid(self, sample_idx: int):
        """Convert diploid genotype to haploid."""
        format = self.format.split(":")
        sample = self.sample_data[sample_idx].split(":")
        sample_dict = dict(zip(format, sample))
        # deepvariant
        if "PL" in format:
            logger.debug(f"Using PL to fix ploidy at {self.chrom}:{self.pos}")
            pls = sample_dict["PL"].split(",")
            pls = [int(pl) for pl in pls]
            if ("|" not in sample_dict["GT"]) and ("/" not in sample_dict["GT"]):
                    logger.critical("Genotype is already haploid.")
            else:
                # update PL, GQ, and GT
                
                num_alleles = len(self.alts) + 1
                un_normalized_probs = [10 ** (pl / -10) for pl in pls]
                homozygous_un_normalized_probs = []
                for i in range(num_alleles):
                    for j in range(i, num_alleles):
                        if i == j:
                            pl_index = int(i * (i + 1) / 2 + j)
                            homozygous_un_normalized_probs.append(
                                un_normalized_probs[pl_index]
                            )
                haploid_probs = tuple(
                    [
                        p / sum(homozygous_un_normalized_probs)
                        for p in homozygous_un_normalized_probs
                    ]
                )
                haploid_pls = [int(-10 * math.log10(p)) for p in haploid_probs]
                min_pl = min(haploid_pls)
                gq = 10000
                haploid_pls = [pl - min_pl for pl in haploid_pls]
                for i, pl in enumerate(haploid_pls):
                    if pl == 0:
                        # maintain a no call
                        gt = "." if sample_dict["GT"] == "./." else i
                    elif pl < gq:
                        gq = pl
                haploid_pls = [str(pl) for pl in haploid_pls]
                sample_dict |= {"PL": ",".join(haploid_pls), "GT": str(gt), "GQ": str(gq)}
        # pbsv
        elif "AD" in format:
            logger.debug(f"Using AD to fix ploidy at {self.chrom}:{self.pos}")
            ads = [int(ad) for ad in sample_dict["AD"].split(",")]
            max_ad = max(ads)
            # genotype unknown if AD is tied, otherwise highest AD
            if ads.count(max_ad) != 1:
                gt = "."
            else:
                gt = ads.index(max_ad)
                sample_dict |= {"GT": str(gt)}
        # sniffles
        elif ("DR" in format) and ("DV" in format):
            logger.debug(f"Using DR+DV to fix ploidy at {self.chrom}:{self.pos}")
            ads = [int(sample_dict["DR"]), int(sample_dict["DV"])]
            max_ad = max(ads)
            if ads.count(max_ad) == 1:
                gt = "."
            else:
                gt = ads.index(max_ad)
                sample_dict |= {"GT": str(gt)}
        else:
            raise ValueError(f"Format {format} is not supported for fixing ploidy.")
        new_sample = [sample_dict[key] for key in format]
        return new_sample

    def fix_ploidy(self, sexes: list[str], 
                   non_diploid_regions: list[tuple[str, int, int]], 
                   ):
        """Fix ploidy for hemizygous loci."""
        if (n_sexes := len(sexes)) != (n_samples := len(self.sample_data)):
            raise ValueError(
                f"Number of sexes {n_sexes} does not match number of \
                             samples {n_samples}"
            )

        if self.original_sample_data != self.sample_data:
            raise ValueError(
                "Sample data has already been modified. Fix \
                             ploidy before shuffling samples or updating genotypes."
            )
        
        if (ploidy_dict := self.check_ploidy(self.chrom, int(self.pos), non_diploid_regions)):
            for sample_idx in range(len(self.sample_data)):
                sex = sexes[sample_idx].lower()
                if "male" in ploidy_dict:
                    if sex in ["m", "male", "xy"]:
                        if ploidy_dict["male"] == 1:
                            new_sample = self.convert_to_haploid(sample_idx)
                            self.sample_data[sample_idx] = ":".join(new_sample)
                        else:
                            logger.warning("Are you sure there should be regions \
                                           with 0 ploidy specified for males?")
                            new_sample = self.convert_to_missing(sample_idx)
                            self.sample_data[sample_idx] = ":".join(new_sample)
                if "female" in ploidy_dict:
                    if sex in ["f", "female", "xx"]:
                        if ploidy_dict["female"] == 1:
                            logger.warning("Are you sure there should be haploid \
                                           regions specified for females?")
                            new_sample = self.convert_to_haploid(sample_idx)
                            self.sample_data[sample_idx] = ":".join(new_sample)
                        else:
                            new_sample = self.convert_to_missing(sample_idx)
                            self.sample_data[sample_idx] = ":".join(new_sample)
                elif sex not in ["m", "male", "xy","f", "female", "xx"]:
                    raise ValueError(
                        f"The specified sex '{sex}' is not valid. Please \
                                     use m, male, or xy for males and f, female, \
                                     or xx for females."
                    )
                self.line = self.variant_data + self.sample_data
