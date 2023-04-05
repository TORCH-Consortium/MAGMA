#!/usr/bin/env python3

import argparse
import pandas as pd

    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Summarize the CALL_WF and MINOR_VARIANTS_ANALYSIS_WF analysis results')
    parser.add_argument('--relabundance_approved_tsv', default="approved_samples.relabundance.tsv", metavar='relabundance_approved_tsv', type=str, help='File enlisting the approved samples from MINOR_VARIANTS_ANALYSIS_WF')
    parser.add_argument('--relabundance_rejected_tsv', default="rejected_samples.relabundance.tsv", metavar='relabundance_rejected_tsv', type=str, help='File enlisting the rejected samples from MINOR_VARIANTS_ANALYSIS_WF')
    parser.add_argument('--call_wf_cohort_stats_tsv', default="joint.cohort_stats.tsv", metavar='call_wf_cohort_stats_tsv', type=str, help='File enlisting the cohort results of CALL_WF')
    parser.add_argument('--output_file', default="joint.merged_cohort_stats.tsv", metavar='output_file', type=str, help='Name of the output file merged cohort statistics')
    args = vars(parser.parse_args())
    
    # Read all CSV files
    df_cohort_stats = pd.read_csv(args['call_wf_cohort_stats_tsv'], sep="\t", index_col="SAMPLE")
    df_approved_relabundance_stats =  pd.read_csv(args['relabundance_approved_tsv'], sep="\t", index_col="SAMPLE")
    df_rejected_relabundance_stats =  pd.read_csv(args['relabundance_rejected_tsv'], sep="\t", index_col="SAMPLE")
    
    # Join the datasets
    df_relabundance_stats_concat = pd.concat([df_approved_relabundance_stats, df_rejected_relabundance_stats])
    df_joint_cohort_stats =  df_cohort_stats.join(df_relabundance_stats_concat, on="SAMPLE", how="left")
  
    # Reorder the columns
    new_cols = ['AVG_INSERT_SIZE', 'MAPPED_PERCENTAGE', 'RAW_TOTAL_SEQS', 'AVERAGE_BASE_QUALITY', 'MEAN_COVERAGE', 'SD_COVERAGE', 'MEDIAN_COVERAGE', 'MAD_COVERAGE', 'PCT_EXC_ADAPTER', 'PCT_EXC_MAPQ', 'PCT_EXC_DUPE', 'PCT_EXC_UNPAIRED', 'PCT_EXC_BASEQ', 'PCT_EXC_OVERLAP', 'PCT_EXC_CAPPED', 'PCT_EXC_TOTAL', 'PCT_1X', 'PCT_5X', 'PCT_10X', 'PCT_30X', 'PCT_50X', 'PCT_100X', 'LINEAGES', 'FREQUENCIES', 'MAPPED_NTM_FRACTION_16S', 'MAPPED_NTM_FRACTION_16S_THRESHOLD_MET', 'COVERAGE_THRESHOLD_MET', 'BREADTH_OF_COVERAGE_THRESHOLD_MET', 'RELABUNDANCE_THRESHOLD_MET', 'ALL_THRESHOLDS_MET']
    df_final_cohort_stats= df_joint_cohort_stats[new_cols]
    df_final_cohort_stats.to_csv(args['output_file'], sep="\t")
