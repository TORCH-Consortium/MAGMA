#!/usr/bin/env python3

import argparse
import csv




def prepare_cohort_stats(merged_cohort_stats_file,output_name)

    cohort_stats_output = []
    open(merged_cohort_stats_file) as tsv_file:
        tsv_reader = csv.reader(tsv_file, delimiter='\t')
        cohort_stats_output.append(("\t").join(tsv_reader[0]))
        for row in tsv_reader[1::]:
            # replace colums 25 - 30
            for column in range(25,30):
                # replace 1 by APPROVED in "*_MET" and 0 by REJECTED
                row[column].replace("1","APPROVED").replace("0","REJECTED")
            cohort_stats_output.append("\t".join(row))
    with open(output_name, 'w') as f:
        f.write("\n".join(cohort_stats_output))



joint.merged_cohort_stats.tsv

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Prepare TSV files for multiqc report generation.')
    parser.add_argument('--merged_cohort_stats', required=True, help='Path to the merged cohort status TSV file.')
    parser.add_argument('--skip_merge_analysis', action="store_false", help='If used, specify that merge analysis skipped.')


    args = parser.parse_args()
    prepare_cohort_stats(args.merged_cohort_stats)

    if skip_merge analysis:
        print("merge was not skipped")

