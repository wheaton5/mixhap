

import argparse
import subprocess
import os

parser = argparse.ArgumentParser(description = "find heterozygous kmers from short read data")

parser.add_argument("--txg_r1s", required=False, nargs='+', help = "10x genomics read1s (takes multiple)")
parser.add_argument("--txg_r2s", required=False, nargs='+', help = "10x genomics read2s (takes multiple)")
parser.add_argument("--hic_r1s", required=False, nargs='+', help = "Hi-C read1s (takes multiple)")
parser.add_argument("--hic_r2s", required=False, nargs='+', help = "Hi-C read2s (takes multiple)")
parser.add_argument("--long_reads", required=False, nargs='+', help= "long reads (pacbio or oxford, takes multiple)")
parser.add_argument("-o", "--output", required=True, help = "prefix output")
parser.add_argument("-t", "--threads", required=True, help="threads")
parser.add_argument("--min_coverage", required=True, type= int, help="min coverage for het kmer")
parser.add_argument("--max_coverage", required=True, type=int, help= "max coverage for het kmer")
parser.add_argument("--max_total_coverage", required=True, type=int, help ="max total coverage for het kmer pair")
parser.add_argument("--max_error", required=False, default=1, type=int, help="max coverage of 3rd kmer with same outer bases")
parser.add_argument("--whitelist", required=True, help="10x barcode whitelist")


parser.add_argument("-m", "--memory", required=False, type=int, default = 24, help="memory in GB, default 24")
args = parser.parse_args()

mypath = os.path.dirname(os.path.realpath(__file__))
cmd = [mypath+"/het_snp_kmers/target/release/het_snp_kmers", "--kmer_counts", args.output+"/kmer_counts.tsv"]
cmd.extend(["--max_coverage", args.max_coverage, "--min_coverage", args.min_coverage, 
            "--max_error", args.max_error, "--max_total_coverage",args.max_total_coverage])
with open(args.output+"/het_kmers.tsv", 'w') as out:
    subprocess.check_call([cmd],stdout=out)

cmd = [mypath+"/molecule_kmers/target/release/molecule_kmers"]
cmd.extend(["--txg_barcodes", args.whitelist, "--txg_r1s"])
cmd.extend(args.txg_r1s)
cmd.extend(['--txg_r1s'])
cmd.extend(args.txg_r2s)
if args.hic_r1s:
    cmd.extend(["--hic_r1s"])
    cmd.extend(args.hic_r1s)
if args.hic_r2s:
    cmd.extend(["--hic_r2s"])
    cmd.extend(args.hic_r2s)

if args.long_reads:
    cmd.extend(["--long_reads"])
    cmd.extend(args.long_reads)
cmd.extend(["--paired_kmers", args.output+"/het_kmers.tsv"])

with open(args.output+"/molecule_kmers.custom_binary", 'wb') as out:
    subprocess.check_call(cmd, stdout = out)


cmd = [mypath+"/mixhap/target/release/mixhap"]
cmd.extend(["--ploidy", "1", "--variants", args.output + "/molecule_kmers.custom_binary"])
cmd.extend(["--barcode_whitelist", args.whitelist])

with open(args.output+"/something_to_be_determined.out", 'w') as out:
    subprocess.check_call(cmd, stdout = out)

