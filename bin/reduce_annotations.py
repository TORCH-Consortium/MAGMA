#!/usr/bin/env python3

import os
import sys
import argparse


#####################################################
# Arg processing
#####################################################


def validate_file(f):
    if not os.path.exists(f):
        raise argparse.ArgumentTypeError("\n\nERROR: The file {0} does not exist.".format(f))
    return f


def parse_args(args=None):
    parser = argparse.ArgumentParser(description="Extract and reshape the annotations for variant recalibration",
                                     epilog="Example usage: reduce_annotations.py -i annotations.log -o new.annotations.txt")

    parser.add_argument("-i",
                        "--input",
                        dest="input_annotations_file",
                        required=True,
                        type=validate_file,
                        help="input file",
                        metavar="FILE")

    parser.add_argument("-o",
                        "--output",
                        dest="output_annotations_file",
                        required=True,
                        help="output file",
                        metavar="FILE")

    parser.add_argument("-t",
                        "--tranches",
                        dest="input_tranches_file",
                        required=True,
                        type=validate_file,
                        help="input tranches file",
                        metavar="FILE")


    parser.add_argument("-d",
                        "--tranches-and-annotations",
                        dest="output_annotations_and_tranches_data",
                        required=True,
                        help="output annotations and tranches data",
                        metavar="FILE")


    args = parser.parse_args()
    return args


#####################################################
# Extract annotations
#####################################################

def extract_reduce_write_annotations(input_filename, output_filename):

    annotations_order_marker = "VariantDataManager - Annotation order is"

    with open(input_filename) as input_annotations_file:
        lines = input_annotations_file.readlines()
        annotations_order_string = list(filter(lambda a_line: a_line.find(annotations_order_marker) != -1, lines))[0]
        annotations_order_list = annotations_order_string.strip().split(":")[-1].replace('[','').replace(']','').split(',')
        print("ANNOTATIONS_ANALYSIS: Initial ordered annotations list: ", annotations_order_list)

        eliminated_annotation_string = '-an ' + ' -an '.join(annotations_order_list[:-1]).strip()
        print("ANNOTATIONS_ANALYSIS: Dropping : " + annotations_order_list[-1])
        print("ANNOTATIONS_ANALYSIS: Final annotation string : ", eliminated_annotation_string)

    with open(output_filename, "w") as output:
        output.write(eliminated_annotation_string)

    return eliminated_annotation_string


#####################################################
# Tranches and annotations file
#####################################################

def tranches_and_annotations(input_annotations_file, input_tranches_filename, output_tranches_file):

    annotations_order_marker = "VariantDataManager - Annotation order is"

    with open(input_annotations_file) as input_annotations_file:
        lines = input_annotations_file.readlines()
        annotations_order_string = list(filter(lambda a_line: a_line.find(annotations_order_marker) != -1, lines))[0]
        annotations_order_list = annotations_order_string.strip().split(":")[-1].replace('[','').replace(']','').split(',')
        print("ANNOTATIONS_ANALYSIS: Ordered annotations list: ", annotations_order_list)


    with open(input_tranches_filename) as input_tranches_file:
        lines = input_tranches_file.readlines()
        minVQSLod_score = lines[-2].split(",")[5]
        print("TRANCHES_ANALYSIS: minVQSLod score corresponding to targetTruthSensitivity: 99.90 is => ", minVQSLod_score)


    with open(output_tranches_file, "w") as output:
        output.writelines("\n".join([minVQSLod_score, annotations_order_list]))

#####################################################
# Main
#####################################################


def main(args=None):
    args = parse_args(args)

    extract_reduce_write_annotations(args.input_annotations_file,
                                     args.output_annotations_file)

    tranches_and_annotations(args.input_annotations_file,
                             args.input_tranches_file,
                             args.output_annotations_and_tranches_data)

#####################################################
#####################################################

if __name__ == "__main__":
    sys.exit(main())

#####################################################
#####################################################
