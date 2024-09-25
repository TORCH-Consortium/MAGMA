#!/usr/bin/env python3
#
# Copyright (c) 2021-2024 MAGMA pipeline authors, see https://doi.org/10.1371/journal.pcbi.1011648
#
# This file is part of MAGMA pipeline, see https://github.com/TORCH-Consortium/MAGMA
#
# For quick overview of GPL-3 license, please refer
# https://www.tldrlegal.com/license/gnu-general-public-license-v3-gpl-3
#
# - You MUST keep this license with original authors in your copy
# - You MUST acknowledge the original source of this software
# - You MUST state significant changes made to the original software
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program . If not, see <http://www.gnu.org/licenses/>.
#


import argparse
import vcf
from scipy.stats import binom_test, binomtest

import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)

def filter_vcf_file(vcf_input, vcf_output, pval):
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

            # The binom_test tests the null hypothesis that ALT variants are fairly distributed across the read directions, the hypothetical success probability of which is 0.5, p=0.5 reflects this hypothetical probability not the significance level.
            p_value = binom_test([dp4_values[-2], dp4_values[-1]], n=sum(dp4_values[-2:]), p=0.5)

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
    parser.add_argument('pval', metavar='pval', type=str, help='The pval to be used for filtering')
    parser.add_argument('input', metavar='input_vcf', type=str, help='The initial VCF file')
    parser.add_argument('output', metavar='output_vcf', type=str, help='The output VCF filed')
    args = vars(parser.parse_args())

    input_pval = float(args['pval'])
    input_vcf_file = args['input']
    output_vcf_file = args['output']

    filter_vcf_file(input_vcf_file, output_vcf_file, input_pval)
