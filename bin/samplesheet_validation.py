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


import re
import argparse
import csv
import json

from sys import exit

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run the MAGMA pipeline samplesheet validation')
    parser.add_argument('input_file', metavar='input_file', type=str, help='The input sample file')
    parser.add_argument('output_file', metavar='output_file', type=str, help='The validate output sample file')
    args = vars(parser.parse_args())

    name_re = re.compile('^[a-zA-Z0-9\-_]*$')

    standard_magma_fields = ['Study', 'Sample', 'Library', 'Attempt', 'R1', 'R2', 'Flowcell', 'Lane', 'Index Sequence']

    # Read the CSV file using the csv library
    with open(args['input_file'], 'r') as infile, open(args['output_file'], 'w', newline='') as outfile:
        reader = csv.DictReader(infile)

        writer = csv.DictWriter(outfile,
                                fieldnames=standard_magma_fields + ['magma_sample_name', 'magma_bam_rg_string'],
                                quoting=csv.QUOTE_NONE,  # Disable quoting
                                escapechar='\\')  # Use backslash as escape character
        writer.writeheader()

        samplesheet_data = []
        fail = False
        for row in reader:

            if 'Library' not in reader.fieldnames or row['Library'] == "" or row['Library'] == None:
                row['Library'] = "1"

            if 'Attempt' not in reader.fieldnames or row['Attempt'] == "" or row['Attempt'] == None:
                row['Attempt'] = "1"

            if 'Flowcell' not in reader.fieldnames or row['Flowcell'] == "" or row['Flowcell'] == None:
                row['Flowcell'] = "1"

            if 'Lane' not in reader.fieldnames or row['Lane'] == "" or row['Lane'] == None:
                row['Lane'] = "1"

            # if 'IndexSequence' not in reader.fieldnames or row['IndexSequence'] == "":
            #     row['IndexSequence'] = "1"

            if 'Index Sequence' not in reader.fieldnames or row['Index Sequence'] == "" or row['Index Sequence'] == None:
                row['Index Sequence'] = "1"

            if 'Study' not in reader.fieldnames or row['Study'] == "":
                row['Study'] = "MAGMA"

            # Perform validation checks
            if 'Study' in reader.fieldnames and not name_re.match(row['Study']):
                print(f'Row {reader.line_num}: {row["Study"]}  Illegal character in STUDY id')
                fail = True
            if not name_re.match(row['Sample']):
                print(f'Row {reader.line_num}: {row["Sample"]} - Illegal character in SAMPLE id')
                fail = True
            if row['R1'] == row['R2']:
                print(f'Row {reader.line_num}: {row["R1"]}, {row["R2"]} - DUPLICATED fastq file specified')
                fail = True

            # Create the magma_sample_name column
            row[
                'magma_sample_name'] = f"{row['Study']}.{row['Sample']}.L{row['Library']}.A{row['Attempt']}.{row['Flowcell']}.{row['Lane']}.{row['Index Sequence']}"

            # Create the magma_bwa_rg_string column
            row[
                'magma_bam_rg_string'] = f"@RG\\tID:{row['Flowcell']}.{row['Lane']}\\tSM:{row['Study']}.{row['Sample']}\\tPL:illumina\\tLB:lib{row['Library']}\\tPU:{row['Flowcell']}.{row['Lane']}.{row['Index Sequence']}"

            # Write the validated row to the output file
            writer.writerow(row)

            # Append the row to the samplesheet_data list
            samplesheet_data.append(row)

        # Convert the samplesheet_data list to JSON and write to a file
        with open('magma_samplesheet.json', 'w') as f:
            json.dump(samplesheet_data, f, indent=4, ensure_ascii=False)

        if not fail:
            print('Samplesheet format validation checks PASSED')
            exit(0)
        else:
            print('Samplesheet format validation checks FAILED')
            exit(1)
