#!/usr/bin/env python3

import os
import sys
import argparse
import csv

def parse_args(args=None):
    parser = argparse.ArgumentParser(description="QC QUANTTB report for a sample ",
                                     epilog="Example usage: quanttb_qc.py -i sampleName.quant.txt -o sampleName.quant.qc.csv -t 0.80 -hdr false -n SAMPLE1234")


    parser.add_argument(dest="quanttb_stats_files",
                        help="list of statsfiles",
                        nargs='+')
    return parser.parse_args()

def main(args=None):
    args = parse_args(args)
    accepted = [["SAMPLE","REFNAME","TOTSCORE","RELABUNDANCE","RELABUNDANCE_THRESHOLD_MET","DEPTH","DERIVED_NAME"]]
    rejected = [["SAMPLE","REFNAME","TOTSCORE","RELABUNDANCE","RELABUNDANCE_THRESHOLD_MET","DEPTH","DERIVED_NAME"]]
    for filename in args.quanttb_stats_files:
        with open(filename, 'r') as quanttb_stats_handle:
            quanttb_stats = csv.reader(quanttb_stats_handle)
            quanttb_stats = list(quanttb_stats)[0]
            if float(quanttb_stats[4]) == 0:
                rejected.append(quanttb_stats)
            elif float(quanttb_stats[4]) == 1:
                accepted.append(quanttb_stats)
            else:
                print('FU')
                print(quanttb_stats)
    print(accepted)
    print(rejected)
    with open('approved_samples.quanttb_cohort_stats.tsv', 'w') as csvfile:
        csvwriter = csv.writer(csvfile, delimiter='\t', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        csvwriter.writerows(accepted)
    with open('rejected_samples.quanttb_cohort_stats.tsv', 'w') as csvfile:
        csvwriter = csv.writer(csvfile, delimiter='\t', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        csvwriter.writerows(rejected)

if __name__ == "__main__":
    sys.exit(main())

