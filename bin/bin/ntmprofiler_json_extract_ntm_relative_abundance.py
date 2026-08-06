import argparse
import json
import math

parser = argparse.ArgumentParser()
parser.add_argument("input_json")
parser.add_argument("output_file")
args = parser.parse_args()

with open(args.input_json) as handle:
    result = json.load(handle)

taxa = result.get("taxa")
if not isinstance(taxa, list):
    raise ValueError("NTMProfiler JSON does not contain a valid 'taxa' list")

ntm_percentage = 0.0

for taxon in taxa:
    species = taxon.get("species")
    abundance = taxon.get("relative_abundance")

    if abundance is None:
        raise ValueError(f"Missing relative_abundance for {species!r}")

    abundance = float(abundance)

    if not math.isfinite(abundance) or abundance < 0 or abundance > 100:
        raise ValueError(
            f"Invalid relative_abundance for {species!r}: {abundance}"
        )

    if species != "Mycobacterium tuberculosis":
        ntm_percentage += abundance

if ntm_percentage > 100:
    raise ValueError(
        f"Summed non-tuberculosis relative abundance exceeds 100: "
        f"{ntm_percentage}"
    )

with open(args.output_file, "w") as handle:
    handle.write(f"{ntm_percentage / 100.0:.10g}\n")
