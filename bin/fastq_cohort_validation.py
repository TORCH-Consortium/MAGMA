#!/usr/bin/env python3

import glob
import argparse
import csv

import pandas as pd

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Summarize the input FASTQ validation report')

    parser.add_argument('joint_vcf_name', metavar='joint_vcf_name', default="joint", type=str,
                        help='The prefix name of joint analysis used for used in the pipeline')

    args = vars(parser.parse_args())

    vcf_name = args['joint_vcf_name']

    # Check for files matching *check.passed*
    passed_files = glob.glob("fastq_validation/*check.passed*")
    passed_data = []  # Store contents of passed files
    if passed_files:
        for fname in passed_files:
            with open(fname) as infile:
                passed_data.append(infile.read())
        with open(f"{vcf_name}.fastqs.passed.tsv", "w") as outfile:
            for data in passed_data:
                outfile.write(data)
    else:
        print("No samples passed!")

    # Create the failed file anyhow, since this is an optional output
    open(f"{vcf_name}.fastqs.failed.tsv", 'a').close()

    # Check for files matching *check.failed*
    failed_files = glob.glob("fastq_validation/*check.failed*")
    failed_data = []  # Store contents of failed files
    if failed_files:
        for fname in failed_files:
            with open(fname) as infile:
                failed_data.append(infile.read())
        with open(f"{vcf_name}.fastqs.failed.tsv", "w") as outfile:
            for data in failed_data:
                outfile.write(data)
    else:
        print("No samples failed!")

    # ============================================
    # Parse the validation reports for exact sample names which passed/failed
    # ============================================

    # Example usage:
    csv_file = 'initial_samplesheet.csv'  # Replace with your CSV file path
    samplesheet_df = pd.read_csv(csv_file)
    # print(samplesheet_df)

    # Create another column by adding Sample and Attempt columns
    samplesheet_df['MagmaSampleName'] = samplesheet_df['Study'].astype(str) + \
                                        "." + samplesheet_df['Sample'].astype(str) + \
                                        ".L" + samplesheet_df['Library'].astype(str) + \
                                        ".A" + samplesheet_df['Attempt'].astype(str) + \
                                        "." + samplesheet_df['Flowcell'].astype(str) + \
                                        "." + samplesheet_df['Lane'].astype(str) + \
                                        "." + samplesheet_df['Index Sequence'].astype(str)

    magma_dict = samplesheet_df.set_index('MagmaSampleName').to_dict(orient='index')

    print(magma_dict)

# ============================================
# Parse the validation reports for exact sample names which passed/failed
# ============================================

validation_and_stats_dict = {}
validate_passed_samples = []

for row in passed_data:
    row_split = row.split("\t")
    derived_magma_name = row_split[0]
    sample_name = derived_magma_name.split(".")[1]
    validate_passed_samples.append(sample_name)
    validation_and_stats_dict[sample_name] = {"magma_name": derived_magma_name,
                                              "fastq_validation_status": "passed"}

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
