#!/usr/bin/env python3

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
