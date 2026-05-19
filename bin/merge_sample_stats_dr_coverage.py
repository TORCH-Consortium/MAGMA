#!/usr/bin/env python3

import argparse
from pathlib import Path

import pandas as pd


SAMPLE_ID_COLUMNS = [
    "sample",
    "sampleName",
    "sample_name",
    "sample_id",
    "Sample",
    "SampleName",
]


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Merge a one-row MAGMA sample stats TSV with a one-row "
            "DR-region coverage TSV."
        )
    )

    parser.add_argument(
        "--sample-name",
        required=True,
        help="Sample name used for validation and output naming.",
    )
    parser.add_argument(
        "--sample-stats",
        required=True,
        type=Path,
        help="Input sample stats TSV.",
    )
    parser.add_argument(
        "--dr-coverage",
        required=True,
        type=Path,
        help="Input DR-region coverage TSV.",
    )
    parser.add_argument(
        "--output",
        required=True,
        type=Path,
        help="Output merged sample stats TSV.",
    )

    return parser.parse_args()


def read_one_row_tsv(path: Path, label: str) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"{label} file does not exist: {path}")

    df = pd.read_csv(path, sep="\t")

    if len(df.index) != 1:
        raise ValueError(
            f"Expected {label} file to contain exactly one row, "
            f"but found {len(df.index)} rows: {path}"
        )

    return df


def main():
    args = parse_args()

    sample_stats = pd.read_csv(args.sample_stats, sep="\t", header=None)
    dr_coverage = pd.read_csv(args.dr_coverage, sep="\t")

    # The DR coverage file may include a sample-identifying first column.
    # Drop it before column-wise concatenation to avoid duplicate IDs.
    dr_coverage = dr_coverage.drop(
        columns=[col for col in SAMPLE_ID_COLUMNS if col in dr_coverage.columns]
    )

    overlapping_columns = set(sample_stats.columns).intersection(dr_coverage.columns)
    if overlapping_columns:
        raise ValueError(
            "Refusing to merge because these DR coverage columns already exist "
            "in the sample stats file: "
            + ", ".join(sorted(overlapping_columns))
        )

    merged = pd.concat(
        [
            sample_stats.reset_index(drop=True),
            dr_coverage.reset_index(drop=True),
        ],
        axis=1,
    )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()
