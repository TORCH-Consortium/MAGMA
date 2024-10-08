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

    parser.add_argument('--cutoff_median_coverage', metavar='cutoff_median_coverage', default=10, type=float, help='The median coverage cutoff threshold')
    parser.add_argument('--cutoff_breadth_of_coverage', metavar='cutoff_breadth_of_coverage', default=0.9, type=float, help='The breadth of coverage cutoff threshold')
    parser.add_argument('--cutoff_ntm_fraction', metavar='cutoff_ntm_fraction', default=0.2, type=float, help='The NTM fraction cutoff threshold')

## NOTE: This is computed by the multiple_infection_filter script
#    parser.add_argument('--cutoff_rel_abundance', metavar='cutoff_rel_abundance', default=0.8, type=float, help='The relative abundance cutoff threshold')

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

    if int(wgsmetrics.loc[0, 'MEDIAN_COVERAGE']) >= args['cutoff_median_coverage']:
        coverage_threshold_met = 1
    else:
        coverage_threshold_met = 0

    if float(wgsmetrics.loc[0, 'PCT_1X']) >= args['cutoff_breadth_of_coverage']:
        breadth_of_coverage_threshold_met = 1
    else:
        breadth_of_coverage_threshold_met = 0

    if ntm_fraction <= args['cutoff_ntm_fraction']:
        ntm_fraction_threshold_met = 1
    else:
        ntm_fraction_threshold_met = 0

    if coverage_threshold_met and breadth_of_coverage_threshold_met and ntm_fraction_threshold_met:
        all_thresholds_met = 1
    else:
        all_thresholds_met = 0

    with open('{}.stats.tsv'.format(args['sample_name']), 'w') as f:
        f.write('\t'.join([str(i) for i in [args['sample_name'], ins_size, mapped_p, total_seqs, avg_qual] + list(wgsmetrics.loc[0, ['MEAN_COVERAGE', 'SD_COVERAGE', 'MEDIAN_COVERAGE', 'MAD_COVERAGE', 'PCT_EXC_ADAPTER', 'PCT_EXC_MAPQ', 'PCT_EXC_DUPE', 'PCT_EXC_UNPAIRED', 'PCT_EXC_BASEQ', 'PCT_EXC_OVERLAP', 'PCT_EXC_CAPPED', 'PCT_EXC_TOTAL', 'PCT_1X', 'PCT_5X', 'PCT_10X', 'PCT_30X', 'PCT_50X', 'PCT_100X']]) + [ntm_fraction, ntm_fraction_threshold_met, coverage_threshold_met, breadth_of_coverage_threshold_met, all_thresholds_met]]))
        f.write('\n')
