#!/usr/bin/env python3

"""
Cluster samples from a snp-dists distance matrix and generates annotation for phylogeny.

Reads a snp-dists tab-separated distance matrix (*.tsv), identifies all samples connected
by a distance <= threshold, computes connected clusters, and writes a two-column
TSV file: Cluster Sample. A treefile can converted to a NEXUS format file with 
unique colours for each cluster, best displayed in FigTree

Usage:
    python cluster_samples.py input_matrix.tsv output_clusters.tsv

Options:
    --threshold FLOAT    Distance threshold (default: 5)
    --tree               The input IQtree treefile 
    --tree-out           Filename for the output NEXUS tree

"""

import argparse
import csv
import re
import os
import colorsys
from collections import defaultdict


def read_distance_matrix(filename):
    """Read a square distance matrix from a TSV file."""

    with open(filename, newline="") as f:
        reader = csv.reader(f, delimiter="\t")
        rows = list(reader)

    # First row contains column sample names
    header = rows[0][1:]
    samples = header

    matrix = {}

    for row in rows[1:]:
        row_sample = row[0]
        values = row[1:]

        if len(values) != len(samples):
            raise ValueError(
                f"Row '{row_sample}' has {len(values)} values, expected {len(samples)}."
            )

        matrix[row_sample] = {}

        for col_sample, value in zip(samples, values):
            try:
                matrix[row_sample][col_sample] = float(value)
            except ValueError:
                raise ValueError(
                    f"Invalid distance value '{value}' between "
                    f"{row_sample} and {col_sample}"
                )

    return samples, matrix


def build_graph(samples, matrix, threshold):
    """Build adjacency graph using the distance threshold."""
    graph = defaultdict(set)

    for s in samples:
        graph[s]  # ensure node exists

    for i, s1 in enumerate(samples):
        for s2 in samples[i + 1:]:
            if matrix[s1][s2] <= threshold:
                graph[s1].add(s2)
                graph[s2].add(s1)

    return graph


def connected_components(graph):
    """Return connected components using DFS."""
    visited = set()
    clusters = []

    for node in sorted(graph):
        if node in visited:
            continue

        stack = [node]
        component = []

        while stack:
            current = stack.pop()
            if current in visited:
                continue

            visited.add(current)
            component.append(current)

            for neighbor in graph[current]:
                if neighbor not in visited:
                    stack.append(neighbor)

        clusters.append(sorted(component))

    return clusters


def write_clusters(outfile, clusters):
    """Write clusters to TSV, excluding singleton clusters."""
    with open(outfile, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")

        writer.writerow(["Cluster", "Sample"])

        cluster_id = 1
        for members in clusters:
            if len(members) < 2:
                continue

            for sample in members:
                writer.writerow([cluster_id, sample])

            cluster_id += 1
def write_sample_cluster_files(outfile, clusters, query_samples, threshold):
    output_dir = os.path.dirname(os.path.abspath(outfile))

    cluster_lookup = {
        sample: members
        for members in clusters
        for sample in members
    }

    for sample in query_samples:
        members = cluster_lookup.get(sample, [])

        cluster_mates = [
            other_sample
            for other_sample in members
            if other_sample != sample
        ]

        output_file = os.path.join(
            output_dir,
            f"{sample}.{threshold}SNPcluster.txt",
        )

        with open(output_file, "w") as handle:
            for cluster_mate in cluster_mates:
                handle.write(f"{cluster_mate}\n")
				
def write_coloured_tree(treefile, outfile, clusters, samples):
    """Write a Nexus tree with taxon labels coloured by cluster."""

    # Sample -> cluster
    sample_to_cluster = {}
    for cid, members in enumerate(clusters, start=1):
        if len(members) > 1:
            for sample in members:
                sample_to_cluster[sample] = cid

    # Generate one colour per cluster
    cluster_colour = {}

    n = max(len(clusters), 1)

    for cid in range(1, n + 1):
        h = (cid - 1) / n
        r, g, b = colorsys.hsv_to_rgb(h, 0.75, 0.95)
        cluster_colour[cid] = "#{:02X}{:02X}{:02X}".format(
            int(r * 255),
            int(g * 255),
            int(b * 255)
        )

    # Generate taxon labels
    taxlabels = []
    for sample in sorted(samples):
        if sample in sample_to_cluster:
            colour = cluster_colour[sample_to_cluster[sample]]
        else:
            colour = "#000000"   # singleton = black

        taxlabels.append(f"    {sample}[&!color={colour}]")

    # Read tree
    with open(treefile) as f:
        tree = f.read().strip()

    # Replace each taxon label
    for sample in samples:
        if sample in sample_to_cluster:
            colour = cluster_colour[sample_to_cluster[sample]]
        else:
            colour = "#000000"

        tree = re.sub(
            rf'(?<![\w.-]){re.escape(sample)}(?=[:),;])',
            f'{sample}[&!color={colour}]',
            tree
        )

    # Write Nexus
    with open(outfile, "w") as out:
        out.write("#NEXUS\n\n")

        # Taxa block (colours tip labels)
        out.write("Begin taxa;\n")
        out.write(f"    Dimensions ntax={len(samples)};\n")
        out.write("    Taxlabels\n")
        out.write("\n".join(taxlabels))
        out.write("\n    ;\n")
        out.write("End;\n\n")

        # Trees block (colours branches)
        out.write("Begin trees;\n")
        out.write(f"    Tree TREE1 = {tree}\n")
        out.write("End;\n")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input_tsv", help="Input distance matrix TSV")
    parser.add_argument("output_tsv", help="Output cluster TSV")
    parser.add_argument(
        "--threshold",
        type=float,
        default=5.0,
        help="Maximum distance for connecting samples (default: 5)",
    )
    parser.add_argument(
        "--tree",
        help="Input Newick tree to colour by cluster."
    )
    parser.add_argument(
        "--tree-out",
        help="Output coloured Nexus tree."
    )
    parser.add_argument(
        "--query-samples",
        required=True,
        help=(
            "Comma-separated sample IDs for which .cluster.txt "
            "files should be written."
        ),
    )

    args = parser.parse_args()

    query_samples = [
        sample.strip()
        for sample in args.query_samples.split(",")
        if sample.strip()
    ]

    samples, matrix = read_distance_matrix(args.input_tsv)

	print(f"Distance matrix file: {args.input_tsv}")
	print(f"Distance matrix samples: {len(samples)}")
	print(f"First samples: {samples[:5]}")

    graph = build_graph(samples, matrix, args.threshold)
    clusters = connected_components(graph)

    exported_clusters = [c for c in clusters if len(c) > 1]

    write_clusters(args.output_tsv, exported_clusters)
    write_sample_cluster_files(
        args.output_tsv,
        clusters,
        query_samples,
		args.threshold,
    )
	
    if args.tree:
        if not args.tree_out:
            raise ValueError("--tree-out must be supplied when using --tree")

        write_coloured_tree(
            args.tree,
            args.tree_out,
            clusters,
            samples
        )

    print(f"Found {len(clusters)} connected samples.")
    print(f"Exported {len(exported_clusters)} clusters containing "
	      f"{sum(len(c) for c in exported_clusters)} samples.")
    print(f"Results written to {args.output_tsv}")


if __name__ == "__main__":
    main()
