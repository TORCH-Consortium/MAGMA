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

import os
import sys
import argparse
import json
from pathlib import Path


#####################################################
# Arg processing
#####################################################


def validate_path(f):
    if not os.path.exists(f):
        raise argparse.ArgumentTypeError(
            "\n\nERROR: The path {0} does not exist.".format(f)
        )
    return f


def parse_args(args=None):
    parser = argparse.ArgumentParser(
        description="Select the best set of annotations based on the minimal minVQSLod score",
        epilog="Example usage: select_best_annotations.py --input_directory annotations_and_tranches_json_files --output_directory best_annotation_files --output_json_file_name all_annotations_data.json",
    )

    parser.add_argument(
        "-i",
        "--input_directory",
        dest="input_folder_with_json_files",
        required=True,
        type=validate_path,
        help="input directory",
        metavar="FOLDER",
    )

    parser.add_argument(
        "-o",
        "--output_directory",
        dest="output_folder_with_best_annotations",
        required=True,
        type=validate_path,
        help="output directory",
        metavar="FOLDER",
    )

    parser.add_argument(
        "-j",
        "--output_json_file_name",
        dest="output_json_file_name",
        required=True,
        help="output directory",
        metavar="FOLDER",
    )

    args = parser.parse_args()
    return args


#####################################################
# Find and move files related to the best annotations
#####################################################


def move_files(annotations_dict, output_directory):

    best_annotations_count = annotations_dict["annotationsCount"]

    p = Path(".")
    grep_pattern = "*" + best_annotations_count + "*"

    best_ann_tranches_files = list(p.glob("./tranches_files/" + grep_pattern))
    print("BEST ANNOTATIONS TRANCHES FILE: ", best_ann_tranches_files)

    for tranches_file in best_ann_tranches_files:
        tranches_file.rename("./" + output_directory + "/" + tranches_file.name)

    best_ann_recalibrated_vcf_files = list(p.glob("./recal_vcf_files/" + grep_pattern))
    print("BEST ANNOTATIONS RECALIBRATED VCF FILES: ", best_ann_recalibrated_vcf_files)

    for recal_vcf_file in best_ann_recalibrated_vcf_files:
        recal_vcf_file.rename("./" + output_directory + "/" + recal_vcf_file.name)


def read_annotations_json_file(annotation_data_json_file_path):

    with open(annotation_data_json_file_path, "r") as annotation_json_data_file:
        annotation_data_dict = json.load(annotation_json_data_file)

    return annotation_data_dict


def collect_annotations_data(base_folder, list_of_files):

    annotations_dict_list = []

    for annotation_json_file in list_of_files:
        annotation_data_json_file_path = base_folder + "/" + annotation_json_file

        annotations_dict_list.append(
            read_annotations_json_file(annotation_data_json_file_path)
        )

    return annotations_dict_list


def find_best_annotations(input_folder, json_file_name):

    # NOTE: Explore the use of path.iterdir
    list_of_annotations_files = os.listdir(input_folder)

    all_annotations_dict_list = collect_annotations_data(
        input_folder, list_of_annotations_files
    )

    print("ALL ANNOTATIONS DATA: ", all_annotations_dict_list)

    with open(json_file_name, "w") as fp:
        json.dump(all_annotations_dict_list, fp, indent=4)

    tentative_best_annotations = all_annotations_dict_list[0]
    tentative_max_minvqslod_score = float(tentative_best_annotations["minVQSLod"])

    for annotations_dict in all_annotations_dict_list[1:]:

        candidate_minvqslod = float(annotations_dict["minVQSLod"])

        tentative_max_minvqslod_score = max(
            tentative_max_minvqslod_score, candidate_minvqslod
        )

    best_annotations_dict = list(
        filter(
            lambda a_dict: float(a_dict["minVQSLod"]) == tentative_max_minvqslod_score,
            all_annotations_dict_list,
        )
    )[0]

    return best_annotations_dict


#####################################################
# Main
#####################################################


def main(args=None):
    args = parse_args(args)

    best_annotations_dict = find_best_annotations(
        args.input_folder_with_json_files, args.output_json_file_name
    )

    print("BEST ANNOTATIONS: ", best_annotations_dict)

    move_files(best_annotations_dict, args.output_folder_with_best_annotations)


#####################################################
#####################################################

if __name__ == "__main__":
    sys.exit(main())

#####################################################
#####################################################
