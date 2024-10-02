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

def validate_study(study, line_num):
    name_re = re.compile('^[a-zA-Z0-9\-_]*$')
    if not name_re.match(study):
        print(f'Row {line_num}: {study} - Illegal character in STUDY id')
        return False
    return True

def validate_sample(sample, line_num):
    name_re = re.compile('^[a-zA-Z0-9\-_]*$')
    if not name_re.match(sample):
        print(f'Row {line_num}: {sample} - Illegal character in SAMPLE id')
        return False
    return True

def validate_fastq_files(r1, r2, line_num):
    if r1 == r2:
        print(f'Row {line_num}: {r1}, {r2} - DUPLICATED fastq file specified')
        return False
    return True

def set_default_value(field, default_value):
    return field if field and field.strip() else default_value

def process_row(row, line_num, reader):
    row['Library'] = set_default_value(row.get('Library'), "1")
    row['Attempt'] = set_default_value(row.get('Attempt'), "1")
    row['Flowcell'] = set_default_value(row.get('Flowcell'), "1")
    row['Lane'] = set_default_value(row.get('Lane'), "1")
    row['Index Sequence'] = set_default_value(row.get('Index Sequence'), "1")
    row['Study'] = set_default_value(row.get('Study'), "MAGMA")

    valid = True
    valid &= validate_study(row['Study'], line_num)
    valid &= validate_sample(row['Sample'], line_num)
    valid &= validate_fastq_files(row['R1'], row['R2'], line_num)

    # Create the magma_sample_name column
    row['magma_sample_name'] = f"{row['Study']}.{row['Sample']}.L{row['Library']}.A{row['Attempt']}.{row['Flowcell']}.{row['Lane']}.{row['Index Sequence']}"

    # Create the magma_bwa_rg_string column
    row['magma_bam_rg_string'] = f"@RG\\tID:{row['Flowcell']}.{row['Lane']}\\tSM:{row['Study']}.{row['Sample']}\\tPL:illumina\\tLB:lib{row['Library']}\\tPU:{row['Flowcell']}.{row['Lane']}.{row['Index Sequence']}"

    return row, valid

def main():
    parser = argparse.ArgumentParser(description='Run the MAGMA pipeline samplesheet validation')
    parser.add_argument('input_file', metavar='input_file', type=str, help='The input sample file')
    parser.add_argument('output_file', metavar='output_file', type=str, help='The validate output sample file')
    args = vars(parser.parse_args())

    standard_magma_fields = ['Study', 'Sample', 'Library', 'Attempt', 'R1', 'R2', 'Flowcell', 'Lane', 'Index Sequence']

    with open(args['input_file'], 'r') as infile, open(args['output_file'], 'w', newline='') as outfile:
        reader = csv.DictReader(infile)
        writer = csv.DictWriter(outfile,
                                fieldnames=standard_magma_fields + ['magma_sample_name', 'magma_bam_rg_string'],
                                quoting=csv.QUOTE_NONE,  # Disable quoting
                                escapechar='\\')  # Use backslash as escape character
        writer.writeheader()

        fail = False
        samplesheet_data = []
        for row in reader:
            row, valid = process_row(row, reader.line_num, reader)
            if not valid:
                fail = True
            writer.writerow(row)
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

if __name__ == '__main__':
    main()