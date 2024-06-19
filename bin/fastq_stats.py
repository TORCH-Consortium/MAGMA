#!/usr/bin/env python3

import ast
import argparse
import re

import pandas as pd

re_mapped_p = re.compile(r'\d* mapped \((.*)%\)')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Process the sample stats')
    parser.add_argument('--sample_name', dest='sample_name', required=True, metavar='sample_name', type=str, help='The sample name')
    parser.add_argument('--seqkit_stats_file', dest='seqkit_stats_file', required=True, metavar='seqkit_stats_file', type=str, help='The seqkit stats file')
    parser.add_argument('--du_stats_file', dest='du_stats_file', required=True, metavar='du_stats_file', type=str, help='The du stats file')
    parser.add_argument('--md5sum_stats_file', dest='md5sum_stats_file', required=True, metavar='md5sum_stats_file', type=str, help='The md5sum metrics file')
    args = vars(parser.parse_args())

    with open(args['wgsmetrics_file']) as f:
        for line in f:
            if '## METRICS CLASS' in line:
                rows = [f.readline().strip(), f.readline().strip()]
                wgsmetrics = pd.DataFrame([rows[1].split('\t')], columns=rows[0].split('\t'))

    with open(args['ntmfraction_file']) as f:
        ntm_fraction = float(f.read().strip())


    with open(args['flagstat_file']) as f:
        for line in f:
            m = re_mapped_p.match(line)
            if m:
                mapped_p = float(m[1])

    if int(wgsmetrics.loc[0, 'MEDIAN_COVERAGE']) >= args['median_coverage_cutoff']:
        coverage_threshold_met = 1
    else:
        coverage_threshold_met = 0


    with open('{}.fastq_stats.tsv'.format(args['sample_name']), 'w') as f:
        f.write('\t'.join([str(i) for i in [args['sample_name'], ins_size, mapped_p, total_seqs, avg_qual] + list(wgsmetrics.loc[0, ['MEAN_COVERAGE', 'SD_COVERAGE', 'MEDIAN_COVERAGE', 'MAD_COVERAGE', 'PCT_EXC_ADAPTER', 'PCT_EXC_MAPQ', 'PCT_EXC_DUPE', 'PCT_EXC_UNPAIRED', 'PCT_EXC_BASEQ', 'PCT_EXC_OVERLAP', 'PCT_EXC_CAPPED', 'PCT_EXC_TOTAL', 'PCT_1X', 'PCT_5X', 'PCT_10X', 'PCT_30X', 'PCT_50X', 'PCT_100X']]) + [ntm_fraction, ntm_fraction_threshold_met, coverage_threshold_met, breadth_of_coverage_threshold_met, all_thresholds_met]]))
        f.write('\n')
