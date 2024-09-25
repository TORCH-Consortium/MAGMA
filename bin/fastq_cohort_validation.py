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

import glob
import argparse
import csv
import json

import pandas as pd

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Summarize the input FASTQ validation report')

    parser.add_argument('magma_samplesheet', metavar='magma_samplesheet',  type=str, help='')

    parser.add_argument('merged_fastq_reports', metavar='merged_fastq_reports',  type=str, help='')

    parser.add_argument('magma_analysis', metavar='magma_analysis',  type=str, help='')


    args = vars(parser.parse_args())

    # ============================================
    # Parse the validation reports for exact sample names which passed/failed
    # ============================================

    # Load the JSON file into a dictionary
    with open(args['merged_fastq_reports'], 'r') as f:
        fastq_report_dict = json.load(f)

    with open(args['magma_samplesheet'], 'r') as f:
        magma_samplesheet_json = json.load(f)

    fastq_report_keys_list = list(fastq_report_dict.keys())

    magma_analysis_json = {}

    for elem in magma_samplesheet_json:
        elem["fastq_report"] = {}

        if elem['R1'] is not None:
            fastq_1_name = elem['R1'].split("/")[-1]
            elem["fastqs_approved"] = True
            if fastq_1_name in fastq_report_keys_list:
                elem["fastq_report"][fastq_1_name] = {"file": fastq_report_dict[fastq_1_name]}
            else:
                elem["fastq_report"][fastq_1_name] = {"fastq_utils_check": "failed"}
                elem["fastqs_approved"] = False

        if elem['R2'] != "" :
            fastq_2_name = elem['R2'].split("/")[-1]
            if fastq_2_name in fastq_report_keys_list:
                elem["fastq_report"][fastq_2_name] = {"file": fastq_report_dict[fastq_2_name]}
            else:
                elem["fastq_report"][fastq_2_name] = {"fastq_utils_check": "failed"}
                elem["fastqs_approved"] = False

        magma_analysis_json[elem["magma_sample_name"]] = elem

    with open(args['magma_analysis'], 'w') as f:
        json.dump(magma_analysis_json, f, indent=4, ensure_ascii= False)
