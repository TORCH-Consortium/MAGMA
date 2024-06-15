#!/usr/bin/env python3

import re
import argparse

import pandas as pd

from sys import exit

parser = argparse.ArgumentParser(description='Run the MAGMA pipeline samplesheet validation')
parser.add_argument('input_file', metavar='input_file', type=str, help='The input sample file')
args = vars(parser.parse_args())

name_re = re.compile('^[a-zA-Z0-9\-\_]*$')
ss = pd.read_csv(args['input_file'])

fail = False
for idx, row in ss.iterrows():
    if not name_re.match(row['Study']):
        print('Row {}: {}  Illegal character in study id'.format(idx, row['Study']))
        fail = True
    if not name_re.match(row['Sample']):
        print('Row {}: {} - Illegal character in sample id'.format(idx, row['Sample']))
        fail = True
    if row['R1'] == row['R2']:
        print('Row {}: {}, {} - Same fastq file specified twice'.format(idx, row['R1'], row['R2']))
        fail = True

if not fail:
    ss.to_csv('samplesheet.valid.csv', index=False)
    print('Samplesheet validation checks passed')
    exit(0)
else:
    print('Samplesheet validation checks failed')
    exit(1)
