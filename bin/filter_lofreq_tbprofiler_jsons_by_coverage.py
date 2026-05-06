#!/usr/bin/env python3

import argparse
import shutil
from pathlib import Path

import pandas as pd


def sample_from_tbprofiler_json(path: Path) -> str:
    name = path.name

    for suffix in [".results.json", ".json"]:
        if name.endswith(suffix):
            return name[: -len(suffix)]

    return path.stem


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--resistance-jsons", nargs="+", required=True)
    parser.add_argument("--cohort-stats", required=True)
    parser.add_argument("--cutoff", type=float, required=True)
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    stats = pd.read_csv(args.cohort_stats, sep="\t")

    keep_samples = set(
        stats.loc[
            stats["MEDIAN_COVERAGE"].astype(float) >= args.cutoff,
            "SAMPLE"
        ].astype(str)
    )

    for json_path in map(Path, args.resistance_jsons):
        sample = sample_from_tbprofiler_json(json_path)

        if sample in keep_samples:
            shutil.copy2(json_path, outdir / json_path.name)


if __name__ == "__main__":
    main()
