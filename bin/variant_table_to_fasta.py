#!/usr/bin/env python3

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
            if l.count('-')/len(l) < (1-args.site_representation_cutoff):
                table.append(l)
            else:
                pass
    with open(args.output_fasta, 'w') as fasta_file:
        for l in list(map(list, zip(*table))):
            fasta_file.write('>{}\n{}\n'.format(l[0].replace('.GT', ''), ''.join(l[1:])))

    

parser = argparse.ArgumentParser(description='tbprofiler script',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('table', type=str, help='The input table to convert (stdin if empty)')
parser.add_argument('output_fasta', type=str, help='The output fasta file')
parser.add_argument('site_representation_cutoff', type=float, help='Minimum fraction of samples that need to have a call at a site before it is considered')
parser.set_defaults(func=main)
args = parser.parse_args()
args.func(args)
