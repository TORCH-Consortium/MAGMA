#!/usr/bin/env python3

import argparse
import pandas as pd


def has_semicolon(value) -> bool:
    if pd.isna(value):
        return False
    return ";" in str(value)


parser = argparse.ArgumentParser()
parser.add_argument("--merged-cohort-stats", required=True)
parser.add_argument("--lofreq-cutoff", required=True, type=float)
parser.add_argument("--output", default="warnings.txt")
args = parser.parse_args()

df = pd.read_csv(args.merged_cohort_stats, sep="\t")

warnings = []

for _, row in df.iterrows():
    sample = str(row["SAMPLE"])

    if float(row["MEDIAN_COVERAGE"]) < args.lofreq_cutoff:
        warnings.append(
            f"WARNING: Cannot reliably report low frequency variants in this {sample} "
            f"due to median coverage less than {args.lofreq_cutoff:g}x"
        )

    if float(row.get("MAPPED_NTM_FRACTION_16S", 0)) > 0:
        warnings.append(
            f"WARNING: Presence of NTM detected in {sample}, refer to the QC stats for additional information"
        )

    if has_semicolon(row.get("LINEAGES")) or has_semicolon(row.get("FREQUENCIES")):
        warnings.append(
            f"WARNING: {sample} may contain multiple Mtb lineages, refer to the QC stats for additional information"
        )

with open(args.output, "w") as out:
    if warnings:
        out.write("\n".join(warnings) + "\n")
