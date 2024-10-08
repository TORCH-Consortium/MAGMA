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

import ast
import argparse

import pandas as pd

header_formats = """##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=MIN_DP,Number=1,Type=Integer,Description="Minimum DP observed within the GVCF block">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification">
"""

def read_vcf(filename):
    with open(filename, 'r') as vcf:
        header = ''
        for line in vcf:
            if line[:2] != '##':
                break
            else:
                header += line
        try:
            df = pd.read_csv(vcf, header=None, sep='\t')
            df.columns = line[:-1].split('\t')
            not_empty = True
        except pd.errors.EmptyDataError as e:
            df = pd.DataFrame(columns=line[:-1].split('\t'))
            not_empty = False
    return df, header, not_empty

def write_vcf(filename, df, header):
    with open(filename, 'w') as vcf:
        vcf.write(header)
        vcf.write(header_formats)
        df.to_csv(vcf, sep='\t', index=False)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Analyse resistance output from the MAGMA pipeline')
    parser.add_argument('lofreq_vcf_file', metavar='lofreq_vcf_file', type=str, help='The input lofreq vcf file')
    parser.add_argument('lofreq_sample_name', metavar='lofreq_sample_name', type=str, help='The sample name')
    parser.add_argument('outfile', metavar='outfile', type=str, help='The name of the output VCF file')

    args = vars(parser.parse_args())

    vcf, header, not_empty = read_vcf(args['lofreq_vcf_file'])
    header = '\n'.join([i for i in header.split('\n') if 'lofreq' not in i])
    if not_empty:
        vcf['FORMAT'] = 'GT:AD:DP:GQ:PL'

        for idx, row in vcf.iterrows():
            info = [ast.literal_eval(i.split('=')[1]) for i in row['INFO'].split(';')[:4]]
            ref_dp = sum(info[3][:2])
            alt_dp = sum(info[3][2:])
            GT = 1
            AD = '{},{}'.format(ref_dp, alt_dp)
            DP = sum(info[3])
            GQ = 99
            PL = '1800,0'
            vcf.loc[idx, args['lofreq_sample_name']] = '{}:{}:{}:{}:{}'.format(GT,AD,DP,GQ,PL)

    write_vcf(args['outfile'], vcf, header)
