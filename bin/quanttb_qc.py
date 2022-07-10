#!/usr/bin/env python3

import os
import sys
import argparse
import csv


#####################################################
# Arg processing
#####################################################


def validate_file(f):
    if not os.path.exists(f):
        raise argparse.ArgumentTypeError("\n\nERROR: The file {0} does not exist.".format(f))
    return f


def parse_args(args=None):
    parser = argparse.ArgumentParser(description="QC QUANTTB report for a sample ",
                                     epilog="Example usage: quanttb_qc.py -i sampleName.quant.txt -o sampleName.quant.qc.csv -t 0.80 -hdr false -n SAMPLE1234")


    parser.add_argument("-n",
                        "--derived_sample_name",
                        dest="derived_sample_name",
                        required=True,
                        help="the unique derived name of the sample")


    parser.add_argument("-i",
                        "--input",
                        dest="input_quanttb_stats_file",
                        required=True,
                        type=validate_file,
                        help="input file",
                        metavar="FILE")

    parser.add_argument("-o",
                        "--output",
                        dest="output_quanttb_qc_file",
                        required=True,
                        help="output file",
                        metavar="FILE")


    parser.add_argument("-t",
                        "--relative_abundance_threshold",
                        dest="relative_abundance_threshold",
                        required=True,
                        help="relative abundance threshold value",
                        type=float)

    parser.add_argument("-hdr",
                        "--write_header",
                        dest="write_header",
                        required=True,
                        help="write header to the output file",
                        choices=("true", "false"))



    args = parser.parse_args()
    return args


#####################################################
# QC QuantTB report
#####################################################

def read_compute_write_qc_report(input_filename, output_filename, relabundance_threshold, write_header, derived_name):

    with open(input_filename) as input_quanttb_stats_file:

        #NOTE  Expected default headers have extra spaces in the beginning
        # ['sample', ' refname', ' totscore', ' relabundance', ' depth']

        input_quanttb_stats_orddict = csv.DictReader(input_quanttb_stats_file, delimiter=',')

        relabundance_threshold_met = 0

        max_relabundance = ['NA', 'NA', 'NA', 'NA', 0, 'NA', derived_name] # Initialize default values in cases where QUANTTB returns an empty file
        for row in input_quanttb_stats_orddict:
            sample = row['sample']
            refname = row[' refname']
            totscore = row[' totscore']
            relabundance = float(row[' relabundance'])
            depth = float(row[' depth'])

            if max_relabundance[3] == 'NA' or relabundance > max_relabundance[3]: # If there is no relabundance yet, or the current relabundance is higher than the current max, update the max_relabundance parameter
                max_relabundance = [sample, refname, totscore, relabundance, 0, depth, derived_name]

        if max_relabundance[3] != 'NA' and max_relabundance[3] >= relabundance_threshold:
            max_relabundance[4] = 1
        print("Sample name: ", max_relabundance[0])
        print("Derived sample name: ", max_relabundance[6])

        print("Relabundance threshold: ", max_relabundance[3])
        print("Threshold met: ", max_relabundance[4])

        #NOTE: Expected output field names, this is then used in the nextflow layer.
        output_fieldnames = ['sample', 'refname', 'totscore', 'relabundance', 'relabundance_threshold_met', 'depth', 'derived_name']

        with open(output_filename, "w") as output_quanttb_qc_csv:
            writer = csv.writer(output_quanttb_qc_csv, delimiter=',')

            if write_header == "true":
                writer.writerow(output_fieldnames)
            writer.writerow(max_relabundance)


#####################################################
# Main
#####################################################


def main(args=None):
    args = parse_args(args)
    read_compute_write_qc_report(args.input_quanttb_stats_file,
                                 args.output_quanttb_qc_file,
                                 args.relative_abundance_threshold,
                                 args.write_header,
                                 args.derived_sample_name)


#####################################################
#####################################################

if __name__ == "__main__":
    sys.exit(main())

#####################################################
#####################################################
