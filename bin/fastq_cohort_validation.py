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

    # # Check for files matching *check.passed*
    # passed_files = glob.glob("fastq_validation/*check.passed*")
    # passed_data = []  # Store contents of passed files
    # if passed_files:
    #     for fname in passed_files:
    #         with open(fname) as infile:
    #             passed_data.append(infile.read())
    #     with open(f"{vcf_name}.fastqs.passed.tsv", "w") as outfile:
    #         for data in passed_data:
    #             outfile.write(data)
    # else:
    #     print("No samples passed!")
    #
    # # Create the failed file anyhow, since this is an optional output
    # open(f"{vcf_name}.fastqs.failed.tsv", 'a').close()
    #
    # # Check for files matching *check.failed*
    # failed_files = glob.glob("fastq_validation/*check.failed*")
    # failed_data = []  # Store contents of failed files
    # if failed_files:
    #     for fname in failed_files:
    #         with open(fname) as infile:
    #             failed_data.append(infile.read())
    #     with open(f"{vcf_name}.fastqs.failed.tsv", "w") as outfile:
    #         for data in failed_data:
    #             outfile.write(data)
    # else:
    #     print("No samples failed!")

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
            if fastq_1_name in fastq_report_keys_list:
                magma_analysis_dict[k]["fastq_report"][fastq_1_name] = {"file": fastq_report_dict[fastq_1_name]}
            else:
                magma_analysis_dict[k]["fastq_report"][fastq_1_name] = {"fastq_utils_check": "failed"}

        if magma_analysis_dict[k]['R2'] is not None:
            fastq_2_name = magma_analysis_dict[k]['R2'].split("/")[-1]
            if fastq_2_name in fastq_report_keys_list:
                magma_analysis_dict[k]["fastq_report"][fastq_2_name] = {"file": fastq_report_dict[fastq_2_name]}
            else:
                magma_analysis_dict[k]["fastq_report"][fastq_2_name] = {"fastq_utils_check": "failed"}

    with open('magma_analysis.json', 'w') as f:
        json.dump(magma_analysis_dict, f, indent=4)

# ============================================
# Parse the validation reports for exact sample names which passed/failed
# ============================================
#
# validation_and_stats_dict = {}
# validate_passed_samples = []
#
# for row in passed_data:
#     row_split = row.split("\t")
#     derived_magma_name = row_split[0]
#     sample_name = derived_magma_name.split(".")[1]
#     validate_passed_samples.append(sample_name)
#     validation_and_stats_dict[sample_name] = {"magma_name": derived_magma_name,
#                                               "fastq_validation_status": "passed"}
#
# validate_failed_samples = []
# for row in failed_data:
#     row_split = row.split("\t")
#     derived_magma_name = row_split[0]
#     sample_name = derived_magma_name.split(".")[1]
#     validate_passed_samples.append(sample_name)
#     validation_and_stats_dict[sample_name] = {"magma_name": derived_magma_name,
#                                               "fastq_validation_status": "failed"}
#

# print(validation_and_stats_dict)
