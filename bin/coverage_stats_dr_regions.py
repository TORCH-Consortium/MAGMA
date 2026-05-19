#!/usr/bin/env python3

import argparse
import csv
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(
        description="Convert samtools bedcov output into one-row DR-region coverage TSV."
    )
    parser.add_argument("--sample-name", required=True)
    parser.add_argument("--bedcov", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args()


def main():
    args = parse_args()

    header = ["sample"]
    values = [args.sample_name]

    with args.bedcov.open() as handle:
        reader = csv.reader(handle, delimiter="\t")

        for row in reader:
            if not row:
                continue

            chrom = row[0]
            bed_start = int(row[1])
            bed_stop = int(row[2])
            summed_depth = float(row[-1])

            region_size = bed_stop - bed_start
            display_start = bed_start + 1

            region_name = f"{chrom}_{display_start}_{bed_stop}"
            mean_depth = summed_depth / region_size if region_size > 0 else "NA"

            header.append(f"dr_region_{region_name}_mean_depth")
            values.append(mean_depth)

    with args.output.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(header)
        writer.writerow(values)


if __name__ == "__main__":
    main()
