#!/usr/bin/env python3

import argparse
import csv

def prepare_cohort_stats(merged_cohort_stats_file,output_name):

    cohort_stats_output = []
    with open(merged_cohort_stats_file) as tsv_file:
        tsv_reader = csv.reader(tsv_file, delimiter='\t')
        cohort_stats_output.append("\t".join(next(tsv_reader, None)))
        for row in tsv_reader:
            # replace colums 25 - 30
            for column in range(26,31):
                # replace 1 by APPROVED in "*_MET" and 0 by REJECTED
                row[column] = row[column].replace("1","APPROVED").replace("0","REJECTED")
            cohort_stats_output.append("\t".join(row))
    with open(output_name, 'w') as f:
        f.write("\n".join(cohort_stats_output))


def prepare_distance_matrix(matrix_file, output_file):
    samples_to_remove =["MTb_L1.SAMN10185847",
                        "MTb_L2.SAMEA1877219",
                        "MTb_L3.SAMEA1877181",
                        "MTb_L410.ERR216945",
                        "MTb_L43.ERR1193883",
                        "MTb_L5.SAMEA1877169",
                        "MTb_L6.SAMEA1877150",
                        "MTb_L7.ERR1971849",
                        "MTb_L8.SRR10828835",
                        "MTb_L9.ERR4162024",
                        "MTB_L10.ERR2707158",
                        "Mcanettii.ERR5104570"]
    with open(matrix_file, "r") as f:
        lines = [line.strip().split("\t") for line in f.readlines()]
        header = lines[0]
        filtered_indices = [i for i, sample in enumerate(header) if sample not in samples_to_remove]
        filtered_matrix = [[header[i] for i in filtered_indices]]

        for row in lines[1:]:
            sample_name = row[0]
        if sample_name not in samples_to_remove:
            filtered_matrix.append([row[i] for i in filtered_indices])

    with open(output_file, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerows(filtered_matrix)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Prepare TSV files for multiqc report generation.')
    parser.add_argument('--merged_cohort_stats', required=True, help='Path to the merged cohort status TSV file.')
    
    parser.add_argument('--output_prefix',required=False,default="prepared", type=str) 
    parser.add_argument('--skip_merge_analysis', action="store_false", help='If used, specify that merge analysis skipped.')
    parser.add_argument('--distance_matrix', required=False, help='Path to the distance matrix TSV file')



    args = parser.parse_args()
    prepare_cohort_output_name = args.output_prefix + "_" + args.merged_cohort_stats 
    prepare_cohort_stats(args.merged_cohort_stats,prepare_cohort_output_name)

    if args.skip_merge_analysis:
        print("merge was not skipped")       
        prepare_distance_matrix_output_name = args.output_prefix + "_" + args.distance_matrix
        prepare_distance_matrix(args.distance_matrix, prepare_distance_matrix_output_name)
        

