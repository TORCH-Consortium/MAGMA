#!/usr/bin/env python3

import argparse
import csv
import re
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(
        description="Convert samtools bedcov output into one-row DR gene coverage TSV."
    )
    parser.add_argument("--sample-name", required=True)
    parser.add_argument("--bedcov", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument(
        "--gene-name-column",
        type=int,
        default=5,
        help=(
            "1-based column index in the BED/bedcov file containing the gene name. "
            "For tbprofiler_whov2plus_genes.bed this is column 5."
        ),
    )
    return parser.parse_args()


def clean_column_name(value: str) -> str:
    value = value.strip()
    value = re.sub(r"[^A-Za-z0-9_.-]+", "_", value)
    value = value.strip("_")
    return value or "unknown_gene"


def main():
    args = parse_args()

    gene_idx = args.gene_name_column - 1

    header = ["sample"]
    values = [args.sample_name]
    seen = {}

    with args.bedcov.open() as handle:
        reader = csv.reader(handle, delimiter="\t")

        for row in reader:
            if not row:
                continue

            if len(row) <= gene_idx:
                raise ValueError(
                    f"Expected gene name column {args.gene_name_column}, "
                    f"but row only has {len(row)} columns: {row}"
                )

            bed_start = int(row[1])
            bed_stop = int(row[2])
            gene_name = row[4]
            summed_depth = float(row[-1])

            region_size = bed_stop - bed_start
            mean_depth = summed_depth / region_size if region_size > 0 else "NA"

            # Guard against duplicated gene names. Your uploaded BED has unique gene names,
            # but this makes the script safe if that changes later.
            seen[gene_name] = seen.get(gene_name, 0) + 1
            if seen[gene_name] == 1:
                column_name = f"dr_gene_{gene_name}_mean_depth"
            else:
                column_name = f"dr_gene_{gene_name}_{seen[gene_name]}_mean_depth"

            header.append(f"dr_gene_{gene_name}_mean_depth")
            values.append(mean_depth)

    with args.output.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(header)
        writer.writerow(values)


if __name__ == "__main__":
    main()
