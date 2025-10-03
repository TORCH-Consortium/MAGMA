#!/usr/bin/env python3
#
# Copyright (c) 2021-2024 MAGMA pipeline authors, see https://doi.org/10.1371/journal.pcbi.1011648
#
# This file is part of MAGMA pipeline, see https://github.com/TORCH-Consortium/MAGMA
#
# For quick overview of GPL-3 license, please refer
# https://www.tldrlegal.com/license/gnu-general-public-license-v3-gpl-3
#
# - You MUST keep this license with original authors in your copy
# - You MUST acknowledge the original source of this software
# - You MUST state significant changes made to the original software
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program . If not, see <http://www.gnu.org/licenses/>.
#
import json
import glob
import pandas as pd
import argparse
import os

# Arguments

parser = argparse.ArgumentParser(description="Global summary of variants from different methods.")
parser.add_argument("--major", required=True, help="Folder containing MajorVariantsXBS")
parser.add_argument("--minor", required=True, help="Folder containing MinorVariantsLoFreq")
parser.add_argument("--delly", required=True, help="Folder containing Delly")
parser.add_argument("--output", "-o", default="cohort_resistance_variants_summary.xlsx", help="Name of the output file (Excel)")

args = parser.parse_args()

folders = {
    "Major": args.major,
    "Minor": args.minor,
    "Delly": args.delly
}

output_file = args.output

# Configuration

categories = [
    "Assoc w R",
    "Assoc w R - Interim",
    "Uncertain significance",
    "Not assoc w R - Interim",
    "Not assoc w R",
    "NA?"
]

drugs_list = [
    "rifampicin", "isoniazid", "ethambutol", "pyrazinamide",
    "moxifloxacin", "levofloxacin", "bedaquiline", "delamanid",
    "pretomanid", "linezolid", "streptomycin", "amikacin",
    "kanamycin", "capreomycin", "clofazimine", "ethionamide"
]

# Reference samples to be ignored
ignored_samples = {
    "Mcanettii.ERR5104570.results",
    "MTB_L10.ERR2707158.results",
    "MTb_L1.SAMN10185847.results",
    "MTb_L2.SAMEA1877219.results",
    "MTb_L3.SAMEA1877181.results",
    "MTb_L410.ERR216945.results",
    "MTb_L43.ERR1193883.results",
    "MTb_L5.SAMEA1877169.results",
    "MTb_L6.SAMEA1877150.results",
    "MTb_L7.ERR1971849.results",
    "MTb_L8.SRR10828835.results",
    "MTb_L9.ERR4162024.results"
}

def make_variant_id(var):
    return (
        var.get("gene_name", ""),
        var.get("nucleotide_change", ""),
        var.get("protein_change", ""),
        var.get("chrom", ""),
        str(var.get("pos", "")),
        var.get("alt", "")
    )

def add_variant(drug, confidence, var, sample_id, method_name,
                summary_occurrences, unique_variants, variant_samples, variant_methods):
    if drug not in drugs_list:
        return
    if confidence not in categories:
        confidence = "NA?"

    summary_occurrences.loc[drug, confidence] += 1
    vid = make_variant_id(var)
    unique_variants[(drug, confidence)].add(vid)

    key = (drug, vid, confidence)
    if key not in variant_samples:
        variant_samples[key] = set()
        variant_methods[key] = set()
    variant_samples[key].add(sample_id)
    variant_methods[key].add(method_name)

def process_folder(folder_path, method_name):
    if not os.path.isdir(folder_path):
        raise NotADirectoryError(f"ERROR: Specified folder does not exist {folder_path}")

    summary_occurrences = pd.DataFrame(0, index=drugs_list, columns=categories)
    summary_unique = pd.DataFrame(0, index=drugs_list, columns=categories)
    unique_variants = { (drug, cat): set() for drug in drugs_list for cat in categories }
    variant_samples = {}
    variant_methods = {}

    json_files = glob.glob(os.path.join(folder_path, "*.json"))
    if not json_files:
        print(f"ERROR: no json found on {folder_path}")
        return summary_occurrences, summary_unique, pd.DataFrame(), 0

    for file in json_files:
        sample_id = os.path.basename(file).replace(".json", "")
        if sample_id in ignored_samples:
            continue 

        with open(file) as f:
            data = json.load(f)

        for var in data.get("dr_variants", []):
            for d in var.get("drugs", []):
                drug = d["drug"]
                conf = d.get("confidence", "NA?")
                add_variant(drug, conf, var, sample_id, method_name,
                            summary_occurrences, unique_variants,
                            variant_samples, variant_methods)

        for var in data.get("other_variants", []):
            for ann in var.get("annotation", []):
                drug = ann.get("drug")
                conf = ann.get("confidence", "NA?")
                if drug:
                    add_variant(drug, conf, var, sample_id, method_name,
                                summary_occurrences, unique_variants,
                                variant_samples, variant_methods)

    for drug in drugs_list:
        for cat in categories:
            summary_unique.loc[drug, cat] = len(unique_variants[(drug, cat)])

    summary_occurrences["Total of Variants Detected"] = summary_occurrences.sum(axis=1)
    summary_unique["Total of Unique Variants Detected"] = summary_unique.sum(axis=1)

    total_samples = len([f for f in json_files if os.path.basename(f).replace(".json", "") not in ignored_samples])
    final_rows = []
    for (drug, vid, conf), samples_set in variant_samples.items():
        gene, nuc_change, prot_change, chrom, pos, alt = vid
        hgvs = f"{gene}_{prot_change}" if prot_change else ""
        mutation_str = "|".join([x for x in vid if x])
        num_samples = len(samples_set)
        percent = round((num_samples / total_samples) * 100, 2) if total_samples > 0 else 0
        methods = ",".join(sorted(variant_methods[(drug, vid, conf)]))
        samples_cleaned = [s.replace(".results", "") for s in sorted(samples_set)]
        sample_list = "[" + ",".join(samples_cleaned) + "]"
        final_rows.append([drug, gene, chrom, pos, mutation_str, hgvs, conf,
                           num_samples, percent, methods, sample_list])

    df_final = pd.DataFrame(final_rows, columns=[
        "Drug", "Gene", "Chromosome", "Position", "Mutation", "HGVS",
        "Classification", "Number of samples", "%", "Method", "Samples"
    ])

    return summary_occurrences, summary_unique, df_final, total_samples

def combine_results(results):
    # Combina ocorrências e variantes únicas
    occ_combined = pd.DataFrame(0, index=drugs_list, columns=categories)
    uniq_combined = pd.DataFrame(0, index=drugs_list, columns=categories)
    df_list = []
    total_samples = 0

    for method_name in results:
        occ, uniq, df_final, n_samples = results[method_name]
        occ_combined += occ
        uniq_combined += uniq
        df_list.append(df_final)
        total_samples += n_samples

    df_final_combined = pd.concat(df_list, ignore_index=True)
    # Opcional: remover duplicatas se necessário
    df_final_combined = df_final_combined.drop_duplicates(subset=["Drug", "Mutation", "Classification", "Samples"])

    occ_combined["Total of Variants Detected"] = occ_combined.sum(axis=1)
    uniq_combined["Total of Unique Variants Detected"] = uniq_combined.sum(axis=1)

    return occ_combined, uniq_combined, df_final_combined

# Process the samples
results = {}
for method_name, folder in folders.items():
    occ, uniq, df_final, n_samples = process_folder(folder, method_name)
    results[method_name] = (occ, uniq, df_final, n_samples)

# Relatório combinado
occ_combined, uniq_combined, df_final_combined = combine_results(results)

# Produce the output
with pd.ExcelWriter(output_file) as writer:
    for method_name in ["Major", "Minor", "Delly"]:
        results[method_name][0].to_excel(writer, sheet_name=f"{method_name}_total_occurrences")
        results[method_name][1].to_excel(writer, sheet_name=f"{method_name}_unique_occurrences")
        results[method_name][2].to_excel(writer, sheet_name=f"{method_name}_variants_summary", index=False)
    # Tabelas combinadas
    occ_combined.to_excel(writer, sheet_name="Combined_total_occurrences")
    uniq_combined.to_excel(writer, sheet_name="Combined_unique_occurrences")
    df_final_combined.to_excel(writer, sheet_name="Combined_variants_summary", index=False)

print(f"Summary Created: {output_file}")