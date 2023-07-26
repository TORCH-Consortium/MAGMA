#!/usr/bin/env python3

import ast
import argparse
import re

import pandas as pd

re_mapped_p = re.compile(r'\d* mapped \((.*)%\)')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Process the sample stats')
    parser.add_argument('--sample_name', dest='sample_name', required=True, metavar='sample_name', type=str, help='The sample name')
    parser.add_argument('--flagstat_file', dest='flagstat_file', required=True, metavar='flagstat_file', type=str, help='The flag stats file')
    parser.add_argument('--samtoolsstats_file', dest='samtoolsstats_file', required=True, metavar='samtoolsstats_file', type=str, help='The samtools stats file')
    parser.add_argument('--wgsmetrics_file', dest='wgsmetrics_file', required=True, metavar='wgsmetrics_file', type=str, help='The WGS metrics file')
    parser.add_argument('--ntmfraction_file', dest='ntmfraction_file', required=True, metavar='ntmfraction_file', type=str, help='The NTM fraction file')

    parser.add_argument('--median_coverage_cutoff', metavar='median_coverage_cutoff', default=10, type=float, help='The median coverage cutoff threshold')
    parser.add_argument('--breadth_of_coverage_cutoff', metavar='breadth_of_coverage_cutoff', default=0.9, type=float, help='The breadth of coverage cutoff threshold')
    parser.add_argument('--rel_abundance_cutoff', metavar='rel_abundance_cutoff', default=0.8, type=float, help='The relative abundance cutoff threshold')
    parser.add_argument('--ntm_fraction_cutoff', metavar='ntm_fraction_cutoff', default=0.2, type=float, help='The NTM fraction cutoff threshold')

    args = vars(parser.parse_args())

    with open(args['wgsmetrics_file']) as f:
        for line in f:
            if '## METRICS CLASS' in line:
                rows = [f.readline().strip(), f.readline().strip()]
                wgsmetrics = pd.DataFrame([rows[1].split('\t')], columns=rows[0].split('\t'))
    with open(args['ntmfraction_file']) as f:
        ntm_fraction = float(f.read().strip())
    with open(args['samtoolsstats_file']) as f:
        for line in f:
            if 'insert size average' in line:
                ins_size = float(line.strip().split('\t')[2])
            if 'raw total sequences' in line:
                total_seqs = int(line.strip().split('\t')[2])
            if 'average quality' in line:
                avg_qual = float(line.strip().split('\t')[2])
    with open(args['flagstat_file']) as f:
        for line in f:
            m = re_mapped_p.match(line)
            if m:
                mapped_p = float(m[1])

    if int(wgsmetrics.loc[0, 'MEDIAN_COVERAGE']) >= args['median_coverage_cutoff']:
        coverage_threshold_met = 1
    else:
        coverage_threshold_met = 0

    if float(wgsmetrics.loc[0, 'PCT_1X']) >= args['breadth_of_coverage_cutoff']:
        breadth_of_coverage_threshold_met = 1
    else:
        breadth_of_coverage_threshold_met = 0

    if ntm_fraction <= args['ntm_fraction_cutoff']:
        ntm_fraction_threshold_met = 1
    else:
        ntm_fraction_threshold_met = 0

    if coverage_threshold_met and breadth_of_coverage_threshold_met and ntm_fraction_threshold_met:
        all_thresholds_met = 1
    else:
        all_thresholds_met = 0

    with open('{}.stats.tsv'.format(args['sample_name']), 'w') as f:
        f.write('\t'.join([str(i) for i in [args['sample_name'], ins_size, mapped_p, total_seqs, avg_qual] + list(wgsmetrics.loc[0, ['MEAN_COVERAGE', 'SD_COVERAGE', 'MEDIAN_COVERAGE', 'MAD_COVERAGE', 'PCT_EXC_ADAPTER', 'PCT_EXC_MAPQ', 'PCT_EXC_DUPE', 'PCT_EXC_UNPAIRED', 'PCT_EXC_BASEQ', 'PCT_EXC_OVERLAP', 'PCT_EXC_CAPPED', 'PCT_EXC_TOTAL', 'PCT_1X', 'PCT_5X', 'PCT_10X', 'PCT_30X', 'PCT_50X', 'PCT_100X']]) + [ntm_fraction, ntm_fraction_threshold_met, coverage_threshold_met, breadth_of_coverage_threshold_met, all_thresholds_met]]))
