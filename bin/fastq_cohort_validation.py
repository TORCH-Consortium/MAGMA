#!/usr/bin/env python3

import glob
import argparse
import csv
import json

import pandas as pd

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Summarize the input FASTQ validation report')

    parser.add_argument('joint_vcf_name', metavar='joint_vcf_name', default="joint", type=str,
                        help='The prefix name of joint analysis used for used in the pipeline')

    args = vars(parser.parse_args())

    vcf_name = args['joint_vcf_name']

    # ============================================
    # Parse the validation reports for exact sample names which passed/failed
    # ============================================

    # Load the JSON file into a dictionary
    with open('merged_fastq_reports.json', 'r') as f:
        fastq_report_dict = json.load(f)

    with open('samplesheet.json', 'r') as f:
        magma_analysis_dict = json.load(f)

    fastq_report_keys_list = list(fastq_report_dict.keys())

    for k in magma_analysis_dict.keys():
        magma_analysis_dict[k]["fastq_report"] = {}

        if magma_analysis_dict[k]['R1'] is not None:
            fastq_1_name = magma_analysis_dict[k]['R1'].split("/")[-1]
            magma_analysis_dict[k]["fastqs_approved"] = True
            if fastq_1_name in fastq_report_keys_list:
                magma_analysis_dict[k]["fastq_report"][fastq_1_name] = {"file": fastq_report_dict[fastq_1_name]}
            else:
                magma_analysis_dict[k]["fastq_report"][fastq_1_name] = {"fastq_utils_check": "failed"}
                magma_analysis_dict[k]["fastqs_approved"] = False

        if magma_analysis_dict[k]['R2'] is not None:
            fastq_2_name = magma_analysis_dict[k]['R2'].split("/")[-1]
            if fastq_2_name in fastq_report_keys_list:
                magma_analysis_dict[k]["fastq_report"][fastq_2_name] = {"file": fastq_report_dict[fastq_2_name]}
            else:
                magma_analysis_dict[k]["fastq_report"][fastq_2_name] = {"fastq_utils_check": "failed"}
                magma_analysis_dict[k]["fastqs_approved"] = False

    with open('magma_analysis.json', 'w') as f:
        json.dump(magma_analysis_dict, f, indent=4)

# ============================================
# Parse the validation reports for exact sample names which passed/failed
# ============================================

    # Filter the dictionary for samples with fastqs_approved == True and False
    approved_samples = {k for k, v in magma_analysis_dict.items() if v["fastqs_approved"] == True}
    # Write approved_samples to a txt file with newline
    with open("approved_samples.txt", "w") as f:
        for sample in approved_samples:
            f.write(sample + "\n")

    rejected_samples = {k for k, v in magma_analysis_dict.items() if v["fastqs_approved"] == False}
    # Write approved_samples to a txt file with newline
    with open("rejected_samples.txt", "w") as f:
        for sample in rejected_samples:
            f.write(sample + "\n")
