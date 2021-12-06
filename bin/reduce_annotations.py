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
                                     epilog="Example usage: reduce_annotations.py -i annotations.log")

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
        print("Initial ordered annotations list: ", annotations_order_list)

        eliminated_annotation_string = '-an ' + ' -an '.join(annotations_order_list[:-1]).strip()
        print("Dropping : " + annotations_order_list[-1])
        print("Final annotation string : ", eliminated_annotation_string)

    with open(output_filename, "w") as output:
        output.write(eliminated_annotation_string)

    return eliminated_annotation_string



#####################################################
# Main
#####################################################


def main(args=None):
    args = parse_args(args)

    extract_reduce_write_annotations(args.input_annotations_file, args.output_annotations_file)


#####################################################
#####################################################

if __name__ == "__main__":
    sys.exit(main())

#####################################################
#####################################################
