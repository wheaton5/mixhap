#!/usr/bin/env python

import argparse

parser = argparse.ArgumentParser(description = "haplotype phasing via sparse mixture model clustering")

parser.add_argument("-v","--vcf", required=True, help = "vcf containing variants called on the bam")
parser.add_argument("-b","--bam", required=True, help = "bam of this data")
parser.add_argument("-f","--reference", required=True, help = "reference fasta")
parser.add_argument("--mode", required=True, help = "10x, Hi-C, or longread")

