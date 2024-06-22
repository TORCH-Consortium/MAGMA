#!/usr/bin/env python3

import glob
import os
import argparse


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Summarize the input FASTQ validation report')

    parser.add_argument('joint_vcf_name', metavar='joint_vcf_name', default="joint", type=str, help='The prefix name of joint analysis used for used in the pipeline')

    args = vars(parser.parse_args())



    # Replace with your actual params.vcf_name
    vcf_name = args['joint_vcf_name']

    # Check for files matching *check.passed*
    passed_files = glob.glob("fastq_validation/*check.passed*")
    if passed_files:
        with open(f"{vcf_name}.fastqs.passed.tsv", "w") as outfile:
            for fname in passed_files:
                with open(fname) as infile:
                    outfile.write(infile.read())
    else:
        print("No samples passed!")

    # Create the failed file anyhow, since this is an optional output
    open(f"{vcf_name}.fastqs.failed.tsv", 'a').close()

    # Check for files matching *check.failed*
    failed_files = glob.glob("fastq_validation/*check.failed*")
    if failed_files:
        with open(f"{vcf_name}.fastqs.failed.tsv", "w") as outfile:
            for fname in failed_files:
                with open(fname) as infile:
                    outfile.write(infile.read())
    else:
        print("No samples failed!")
