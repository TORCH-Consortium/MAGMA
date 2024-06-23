#!/usr/bin/env python3

import re
import argparse
import csv

from sys import exit

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run the MAGMA pipeline samplesheet validation')
    parser.add_argument('input_file', metavar='input_file', type=str, help='The input sample file')
    parser.add_argument('output_file', metavar='output_file', type=str, help='The validate output sample file')
    args = vars(parser.parse_args())

    name_re = re.compile('^[a-zA-Z0-9\-_]*$')

    # Read the CSV file using the csv library
    with open(args['input_file'], 'r') as infile, open(args['output_file'], 'w', newline='') as outfile:
        reader = csv.DictReader(infile)
        writer = csv.DictWriter(outfile, fieldnames=reader.fieldnames + ['magma_sample_name', 'magma_bam_rg_string'])
        writer.writeheader()

        fail = False
        for row in reader:

            # Perform validation checks
            if not name_re.match(row['Study']):
                print(f'Row {reader.line_num}: {row["Study"]}  Illegal character in STUDY id')
                fail = True
            if not name_re.match(row['Sample']):
                print(f'Row {reader.line_num}: {row["Sample"]} - Illegal character in SAMPLE id')
                fail = True
            if row['R1'] == row['R2']:
                print(f'Row {reader.line_num}: {row["R1"]}, {row["R2"]} - DUPLICATED fastq file specified')
                fail = True


            # Create the magma_sample_name column
            row['magma_sample_name'] = f"{row['Study']}.{row['Sample']}.L{row['Library']}.A{row['Attempt']}.{row['Flowcell']}.{row['Lane']}.{row['Index Sequence']}"

            # Create the magma_bwa_rg_string column
            row['magma_bam_rg_string'] = f"@RG\\tID:{row['Flowcell']}.{row['Lane']}\\tSM:{row['Study']}.{row['Sample']}\\tPL:illumina\\tLB:lib{row['Library']}\\tPU:{row['Flowcell']}.{row['Lane']}.{row['Index Sequence']}"


            # Write the validated row to the output file
            writer.writerow(row)

        if not fail:
            print('Samplesheet format validation checks PASSED')
            exit(0)
        else:
            print('Samplesheet format validation checks FAILED')
            exit(1)