#!/usr/bin/env python3

"""
Cluster samples from a snp-dists distance matrix and generate annotations
for a phylogeny.

Reads a snp-dists tab-separated distance matrix (*.tsv), identifies all
samples connected by a distance less than or equal to the threshold,
computes connected clusters, and writes a two-column TSV file containing
Cluster and Sample.

A treefile can also be converted to NEXUS format with unique colours for
each cluster, best displayed in FigTree.

Usage:
    python catnip.py input_matrix.tsv output_clusters.tsv

Options:
    --threshold FLOAT
        Distance threshold. Default: 5.

    --tree FILE
        Input IQ-TREE treefile.

    --tree-out FILE
        Filename for the output NEXUS tree.

    --query-samples IDS
        Comma-separated sample IDs for which per-sample cluster files
        should be written.
"""

import argparse
import colorsys
import csv
import os
import re
from collections import defaultdict


def read_distance_matrix(filename):
    """Read a square distance matrix from a TSV file."""
    with open(filename, newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        rows = list(reader)

    if not rows:
        raise ValueError(
            f"Distance matrix '{filename}' is empty."
        )

    if len(rows[0]) < 2:
        raise ValueError(
            f"Distance matrix '{filename}' has an invalid header: "
            f"{rows[0]!r}"
        )

    # The first cell is the row-name header. Remaining cells are samples.
    samples = rows[0][1:]

    if not samples:
        raise ValueError(
            f"No sample IDs found in distance matrix '{filename}'."
        )

    matrix = {}

    for row_number, row in enumerate(rows[1:], start=2):
        if not row:
            continue

        row_sample = row[0]
        values = row[1:]

        if not row_sample:
            raise ValueError(
                f"Missing row sample name on line {row_number} "
                f"of '{filename}'."
            )

        if len(values) != len(samples):
            raise ValueError(
                f"Row '{row_sample}' has {len(values)} values, "
                f"expected {len(samples)}."
            )

        matrix[row_sample] = {}

        for col_sample, value in zip(samples, values):
            try:
                matrix[row_sample][col_sample] = float(value)
            except ValueError as error:
                raise ValueError(
                    f"Invalid distance value '{value}' between "
                    f"{row_sample} and {col_sample}."
                ) from error

    missing_rows = [
        sample
        for sample in samples
        if sample not in matrix
    ]

    if missing_rows:
        raise ValueError(
            "Distance matrix is missing rows for: "
            + ", ".join(missing_rows)
        )

    return samples, matrix


def build_graph(samples, matrix, threshold):
    """Build an adjacency graph using the distance threshold."""
    graph = defaultdict(set)

    for sample in samples:
        # Ensure singleton samples are represented in the graph.
        graph[sample]

    for index, sample_1 in enumerate(samples):
        for sample_2 in samples[index + 1:]:
            if matrix[sample_1][sample_2] <= threshold:
                graph[sample_1].add(sample_2)
                graph[sample_2].add(sample_1)

    return graph


def connected_components(graph):
    """Return connected components using depth-first search."""
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

            for neighbour in graph[current]:
                if neighbour not in visited:
                    stack.append(neighbour)

        clusters.append(sorted(component))

    return clusters


def write_clusters(outfile, clusters):
    """Write non-singleton clusters to a TSV file."""
    with open(outfile, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["Cluster", "Sample"])

        cluster_id = 1

        for members in clusters:
            if len(members) < 2:
                continue

            for sample in members:
                writer.writerow([cluster_id, sample])

            cluster_id += 1


def write_sample_cluster_files(
    outfile,
    clusters,
    query_samples,
    threshold,
):
    """
    Write one cluster-mate file for each clustered query sample.

    The full distance matrix is used to determine clusters, but files are
    only written for samples supplied through --query-samples.

    No file is written when the requested sample is absent from the matrix
    or has no cluster mates.
    """
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

        if not cluster_mates:
            continue

        output_file = os.path.join(
            output_dir,
            f"{sample}.{threshold:g}SNPcluster.txt",
        )

        with open(output_file, "w") as handle:
            for cluster_mate in cluster_mates:
                handle.write(f"{cluster_mate}\n")


def write_coloured_tree(
    treefile,
    outfile,
    clusters,
    samples,
):
    """Write a NEXUS tree with taxon labels coloured by cluster."""
    sample_to_cluster = {}

    for cluster_id, members in enumerate(clusters, start=1):
        if len(members) > 1:
            for sample in members:
                sample_to_cluster[sample] = cluster_id

    clustered_ids = sorted(set(sample_to_cluster.values()))
    cluster_colour = {}

    number_of_clusters = max(len(clustered_ids), 1)

    for colour_index, cluster_id in enumerate(clustered_ids):
        hue = colour_index / number_of_clusters
        red, green, blue = colorsys.hsv_to_rgb(
            hue,
            0.75,
            0.95,
        )

        cluster_colour[cluster_id] = (
            "#{:02X}{:02X}{:02X}".format(
                int(red * 255),
                int(green * 255),
                int(blue * 255),
            )
        )

    taxlabels = []

    for sample in sorted(samples):
        if sample in sample_to_cluster:
            colour = cluster_colour[
                sample_to_cluster[sample]
            ]
        else:
            colour = "#000000"

        taxlabels.append(
            f"    {sample}[&!color={colour}]"
        )

    with open(treefile) as handle:
        tree = handle.read().strip()

    if not tree:
        raise ValueError(
            f"Tree file '{treefile}' is empty."
        )

    for sample in samples:
        if sample in sample_to_cluster:
            colour = cluster_colour[
                sample_to_cluster[sample]
            ]
        else:
            colour = "#000000"

        tree = re.sub(
            rf"(?<![\w.-]){re.escape(sample)}(?=[:),;])",
            f"{sample}[&!color={colour}]",
            tree,
        )

    with open(outfile, "w") as handle:
        handle.write("#NEXUS\n\n")

        handle.write("Begin taxa;\n")
        handle.write(
            f"    Dimensions ntax={len(samples)};\n"
        )
        handle.write("    Taxlabels\n")
        handle.write("\n".join(taxlabels))
        handle.write("\n    ;\n")
        handle.write("End;\n\n")

        handle.write("Begin trees;\n")
        handle.write(f"    Tree TREE1 = {tree}\n")
        handle.write("End;\n")


def parse_query_samples(value):
    """Parse a comma-separated query-sample argument."""
    query_samples = [
        sample.strip()
        for sample in value.split(",")
        if sample.strip()
    ]

    if not query_samples:
        raise ValueError(
            "--query-samples did not contain any sample IDs."
        )

    return query_samples


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "input_tsv",
        help="Input distance matrix TSV.",
    )

    parser.add_argument(
        "output_tsv",
        help="Output cluster TSV.",
    )

    parser.add_argument(
        "--threshold",
        type=float,
        default=5.0,
        help=(
            "Maximum distance for connecting samples "
            "(default: 5)."
        ),
    )

    parser.add_argument(
        "--tree",
        help="Input Newick tree to colour by cluster.",
    )

    parser.add_argument(
        "--tree-out",
        help="Output coloured NEXUS tree.",
    )

    parser.add_argument(
        "--query-samples",
        required=True,
        help=(
            "Comma-separated sample IDs for which "
            ".SNPcluster.txt files should be written."
        ),
    )

    args = parser.parse_args()

    if args.threshold < 0:
        raise ValueError(
            "--threshold must be zero or greater."
        )

    query_samples = parse_query_samples(
        args.query_samples
    )

    samples, matrix = read_distance_matrix(
        args.input_tsv
    )

    print(
        f"Distance matrix file: {args.input_tsv}"
    )
    print(
        f"Distance matrix samples: {len(samples)}"
    )
    print(
        f"First samples: {samples[:5]}"
    )

    graph = build_graph(
        samples,
        matrix,
        args.threshold,
    )

    clusters = connected_components(graph)

    exported_clusters = [
        cluster
        for cluster in clusters
        if len(cluster) > 1
    ]

    write_clusters(
        args.output_tsv,
        exported_clusters,
    )

    write_sample_cluster_files(
        args.output_tsv,
        clusters,
        query_samples,
        args.threshold,
    )

    if args.tree:
        if not args.tree_out:
            raise ValueError(
                "--tree-out must be supplied when using --tree."
            )

        write_coloured_tree(
            args.tree,
            args.tree_out,
            clusters,
            samples,
        )
    elif args.tree_out:
        raise ValueError(
            "--tree must be supplied when using --tree-out."
        )

    clustered_sample_count = sum(
        len(cluster)
        for cluster in exported_clusters
    )

    print(
        f"Found {len(clusters)} connected components."
    )
    print(
        f"Exported {len(exported_clusters)} clusters "
        f"containing {clustered_sample_count} samples."
    )
    print(
        f"Results written to {args.output_tsv}"
    )


if __name__ == "__main__":
    main()
