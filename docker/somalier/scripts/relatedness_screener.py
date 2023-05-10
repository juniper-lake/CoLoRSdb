#!/usr/bin/env python3
"""
This tool parses the pairs.tsv output of somalier relate to determine which samples to 
drop from a cohort in order to keep relatedness below some threshold.
"""

__version__ = "0.1.0"

from loguru import logger
import pandas as pd
import sys
import argparse
from itertools import chain
from collections import Counter
import os


class NoRelation(object):
    """Object representing a group of samples with relatedness below some threshold."""

    def __init__(
        self,
        infile: str,
        max_relatedness: float = 0.125,
        sample_order: list = None,
        coverages: list = None,
    ):
        super().__init__()
        self.header = [
            "#sample_a",
            "sample_b",
            "relatedness",
            "ibs0",
            "ibs2",
            "hom_concordance",
            "hets_a",
            "hets_b",
            "hets_ab",
            "shared_hets",
            "hom_alts_a",
            "hom_alts_b",
            "shared_hom_alts",
            "n",
            "x_ibs0",
            "x_ibs2",
            "expected_relatedness",
        ]
        self.samples = sample_order
        self.max_relatedness = max_relatedness
        self.pairs = None
        self.dropped_samples = []
        self.coverages = coverages
        self.sorted_keep_drop = None

        if self.max_relatedness < 0.0 or self.max_relatedness > 1.0:
            raise ValueError(
                f"max_relatedness must be between 0.0 and 1.0, but \
                    was {self.max_relatedness}"
            )

        # read the somalier pairs file
        logger.info(f"Reading somalier pairwise relatedness values from file {infile}")
        self.pairs = pd.read_csv(infile, sep="\t", header=0)

        # check that the header is as expected
        if self.pairs.columns.tolist() != self.header:
            raise IOError(
                f"File {infile} does not have the expected header: {self.header}"
            )

        # get sample names from pairs file if order is not provided
        if self.samples is None:
            if self.coverages is not None:
                logger.warning(
                    "Provided coverages will be ignored because sample order is \
                        not provided."
                )
                self.coverages = None
            logger.warning(
                "Sample order is not explicitly set, so it will be determined by the \
                    order of samples in the pairs file."
            )
            self.samples = sorted(
                list(
                    set(
                        self.pairs["#sample_a"].tolist()
                        + self.pairs["sample_b"].tolist()
                    )
                )
            )

        # check that the number of samples matches the number of coverages
        if self.coverages is not None:
            if len(self.coverages) != len(self.samples):
                raise ValueError(
                    f"Number of coverages ({len(self.coverages)}) does not match \
                        number of samples ({len(self.samples)})"
                )
            else:
                try:
                    self.coverages = [float(x) for x in self.coverages]
                except ValueError:
                    raise ValueError(
                        f"Coverages must be integers, but provided coverages \
                            were {self.coverages}"
                    )
                self.samples_sorted_by_coverage = [
                    x
                    for _, x in sorted(zip(self.coverages, self.samples), reverse=False)
                ]

        # check that the samples in the pairs file match the provided samples
        if set(self.samples) != set(
            self.pairs["#sample_a"].tolist() + self.pairs["sample_b"].tolist()
        ):
            raise ValueError(
                f"Provided samples \
                    [{set(self.samples) - set(self.pairs['#sample_a'].tolist() + self.pairs['sample_b'].tolist())}] \
                    do not match samples in pairs file "
            )

        # check if any samples are related
        if self.pairs["relatedness"].max() <= self.max_relatedness:
            logger.info(
                f"All samples are below the relatedness threshold of \
                    {self.max_relatedness}, so no samples will be removed."
            )
        else:
            logger.info(f"Removing samples with relatedness > {self.max_relatedness}")
            self.remove_related_samples()
        self.sorted_keep_drop = [
            "drop" if sample in self.dropped_samples else "keep"
            for sample in self.samples
        ]

    def remove_related_samples(self):
        while (relatedness := self.pairs.relatedness.max()) > self.max_relatedness:
            logger.info(
                f"Maximum pairwise relatedness among samples {relatedness} is above \
                    the threshold {self.max_relatedness}"
            )
            # list of samples with relatedness above threshold, where frequency
            # of each sample name reflects number of relations above threshold
            related_samples = list(
                chain(
                    *(
                        self.pairs.loc[
                            self.pairs.loc[:, "relatedness"] > self.max_relatedness,
                            ["#sample_a", "sample_b"],
                        ].values
                    )
                )
            )
            # find samples with highest number of relations above threshold
            most_common_samples, count = most_common(related_samples)

            if self.coverages is None:
                # if coverages not provided, drop first index of most common samples
                self.dropped_samples.append(most_common_samples[0])
            else:
                # if coverages provided, drop sample with lowest coverage
                self.dropped_samples.append(
                    sorted(
                        most_common_samples, key=self.samples_sorted_by_coverage.index
                    )[0]
                )
            logger.info(
                f"Removing sample [{self.dropped_samples[-1]}], related to {count} \
                    other samples"
            )
            self.pairs = self.pairs.loc[
                ~(self.pairs["#sample_a"] == self.dropped_samples[-1])
                & ~(self.pairs["sample_b"] == self.dropped_samples[-1])
            ]


def most_common(lst):
    """Returns list of most common element in a list and the count of those elements"""
    counts = Counter(lst)
    max_count = counts.most_common(1)[0][1]
    items = [value for value, count in counts.most_common() if count == max_count]
    return items, max_count


def relatedness_float(x, min=0.0, max=1.0):
    """Restrict the range of a float to [0.0, 1.0]"""
    try:
        x = float(x)
    except ValueError:
        raise argparse.ArgumentTypeError(f"{x} is not a floating-point literal")

    if x < min or x > max:
        raise argparse.ArgumentTypeError(f"{x} not in range [{min}, {max}]")
    return x


def coverage_float(x, min=0.0):
    """Restrict the range of a float to above 0"""
    try:
        x = float(x)
    except ValueError:
        raise argparse.ArgumentTypeError(f"{x} is not a floating-point literal")

    if x < min:
        raise argparse.ArgumentTypeError(f"{x} not in range [{min}, inf]")
    return x


def parse_args(args):
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        "somalier_pairs_tsv", help="pairs tsv output from somalier relate"
    )
    parser.add_argument(
        "--max_relatedness",
        dest="max_relatedness",
        type=relatedness_float,
        default=0.125,
        help="maximum relatedness to consider unrelated",
    )
    parser.add_argument(
        "--sample_order",
        dest="sample_order",
        type=str,
        nargs="+",
        help="list of samples names sorted in preferred order of output, \
            should be same order as coverages",
    )
    parser.add_argument(
        "--coverages",
        dest="coverages",
        type=coverage_float,
        nargs="+",
        help="list of coverages used to determine which sample to drop in case of tie, \
            should be same order as sample_order",
    )
    parser.add_argument(
        "--outfile", dest="outfile", help="output file name, should be TSV"
    )
    parser.add_argument(
        "--loglevel",
        dest="loglevel",
        choices=["TRACE", "DEBUG", "INFO", "SUCCESS", "WARNING", "ERROR", "CRITICAL"],
        default="INFO",
        help="set the logging level",
    )
    parser.add_argument(
        "--version",
        action="version",
        version="%(prog)s (version {version})".format(version=__version__),
    )

    return parser.parse_args()


def flag_related_samples(
    somalier_pairs: str,
    max_relatedness: float,
    sample_order: list = None,
    coverages: list = None,
    outfile: str = None,
):
    """Identify related samples to remove, print samples names and keep/drop boolean"""

    # check if outfile exists
    if outfile:
        if not outfile.endswith(".vcf"):
            raise ValueError("Output file should end in .vcf")
        if os.path.isfile(outfile):
            raise IOError(f"Output file {outfile} already exists.")
        else:
            file = open(outfile, "w")
    else:
        logger.info("No output file provided, printing to stdout")
        file = sys.stdout

    cohort = NoRelation(
        somalier_pairs,
        max_relatedness=max_relatedness,
        sample_order=sample_order,
        coverages=coverages,
    )

    # print sample names in provided order
    print("\t".join(cohort.samples), file=file)

    # print keep/drop for each sample
    print("\t".join(cohort.sorted_keep_drop), file=file)

    logger.success("Successfully terminated")


if __name__ == "__main__":
    """This is executed when run from the command line."""

    args = parse_args(sys.argv[1:])

    logger.configure(handlers=[{"sink": sys.stderr, "level": args.loglevel}])

    flag_related_samples(
        args.somalier_pairs_tsv,
        args.max_relatedness,
        args.sample_order,
        args.coverages,
        args.outfile,
    )
