#! /usr/bin/env python3
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

import sys
import argparse

def main(args):
    table = []
    with open(args.table, 'r') as table_file:
        table.append(table_file.readline().strip().split('\t')) # Get the headerline without modifying
        # Process the actual variants
        for idx, l in enumerate(table_file):
            l = l.strip().split('\t')
            l = [i.replace('*', '-').replace('.', '-') for i in l]
            if l.count('-')/len(l) < (1-args.cutoff_site_representation):
                table.append(l)
            else:
                pass
    with open(args.output_fasta, 'w') as fasta_file:
        for l in list(map(list, zip(*table))):
            fasta_file.write('>{}\n{}\n'.format(l[0].replace('.GT', ''), ''.join(l[1:])))

    

parser = argparse.ArgumentParser(description='tbprofiler script',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('table', type=str, help='The input table to convert (stdin if empty)')
parser.add_argument('output_fasta', type=str, help='The output fasta file')
parser.add_argument('cutoff_site_representation', type=float, help='Minimum fraction of samples that need to have a call at a site before it is considered')
parser.set_defaults(func=main)
args = parser.parse_args()
args.func(args)
