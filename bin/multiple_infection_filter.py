import os
import json
import csv

import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Analyse resistance output from XBS Pipeline')
    parser.add_argument('indir', metavar='indir', type=str, help='The directory containing the LoFreq TBProfiler output')
    parser.add_argument('relative_abundance_threshold', metavar='relative_abundance_threshold', type=float, help='Minimum relative abundance of the majority strain required to process the sample')

    args = vars(parser.parse_args())
    samples = []
    for file_name in os.listdir(args['indir']):
        with open(os.path.join(args['indir'], file_name)) as json_file:
            samples.append(json.load(json_file))
    accepted = [['SAMPLE', 'LINEAGES', 'FREQUENCIES', 'RELABUNDANCE_THRESHOLD_MET']]
    rejected = [['SAMPLE', 'LINEAGES', 'FREQUENCIES', 'RELABUNDANCE_THRESHOLD_MET']]
    for sample in samples:
        sublins = [i for i in sample['lineage'] if i['lin'] in sample['sublin'].split(';')]
        lins = [i['lin'] for i in sublins]
        fracs = [i['frac'] for i in sublins]
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
