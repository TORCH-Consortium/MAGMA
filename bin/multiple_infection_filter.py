#!/usr/bin/env python3

import os
import json
import csv
import glob

import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Analyse resistance output from MAGMA pipeline')
    parser.add_argument('indir', metavar='indir', type=str, help='The directory containing the LoFreq TBProfiler output')
    parser.add_argument('relative_abundance_threshold', metavar='relative_abundance_threshold', type=float, help='Minimum relative abundance of the majority strain required to process the sample')

    args = vars(parser.parse_args())
    samples = []
    json_files = glob.glob("*/**.json")
    for file_name in json_files:
        with open(file_name) as json_file:
            samples.append(json.load(json_file))
    accepted = [['SAMPLE', 'LINEAGES', 'FREQUENCIES', 'RELABUNDANCE_THRESHOLD_MET']]
    rejected = [['SAMPLE', 'LINEAGES', 'FREQUENCIES', 'RELABUNDANCE_THRESHOLD_MET']]
    for sample in samples:
        sublins = [i for i in sample['lineage'] if i['lineage'] in sample['sub_lineage'].split(';')]
        lins = [i['lineage'] for i in sublins]
        fracs = [i['fraction'] for i in sublins]
        if not lins:
            rejected.append([sample['id'], 'None', 'None', 0])
        elif max(fracs) < args['relative_abundance_threshold']:
            rejected.append([sample['id'], ';'.join(lins), ';'.join(['{:.0%}'.format(i) for i in fracs]), 0])
        else:
            accepted.append([sample['id'], ';'.join(lins), ';'.join(['{:.0%}'.format(i) for i in fracs]), 1])

    with open('approved_samples.relabundance.tsv', 'w') as csvfile:
        csvwriter = csv.writer(csvfile, delimiter='\t', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        csvwriter.writerows(accepted)

    with open('rejected_samples.relabundance.tsv', 'w') as csvfile:
        csvwriter = csv.writer(csvfile, delimiter='\t', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        csvwriter.writerows(rejected)
