#!/usr/bin/env python3

import argparse
import vcf
from scipy.stats import binom_test


def filter_vcf_file(vcf_input, vcf_output, pval=0.5):

    # Open the VCF file for reading

    vcf_reader = vcf.Reader(open(vcf_input, 'r'))

    # Create a list to store the filtered records

    filtered_records = []

    # Iterate over each record (variant) in the VCF file

    for record in vcf_reader:

        # Extract the DP4 field from the INFO dictionary

        dp4_values = record.INFO.get('DP4')

        # Perform binomial test on the last two values

        if dp4_values:

            p_value = binom_test([dp4_values[-2], dp4_values[-1]], n=sum(dp4_values[-2:]), p=pval)

            # Check if p-value is above or equal to 0.05
            if p_value >= pval:
                filtered_records.append(record)


    # Write the filtered records to a new VCF file
    vcf_writer = vcf.Writer(open(vcf_output, 'w'), vcf_reader)

    for record in filtered_records:
        vcf_writer.write_record(record)

    # Close the VCF writer
    vcf_writer.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Reduce strand bias in VCF file')
    parser.add_argument('input', metavar='input_vcf', type=str, help='The initial VCF file')
    parser.add_argument('output', metavar='output_vcf', type=str, help='The output VCF filed')
    args = vars(parser.parse_args())

    input_vcf_file = args['input']
    output_vcf_file = args['output']

    filter_vcf_file(input_vcf_file, output_vcf_file)
