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
                                     epilog="Example usage: reshape_annotations.py annotations.log")

    parser.add_argument("-i",
                        "--input",
                        dest="annotations_file",
                        required=True,
                        type=validate_file,
                        help="input file",
                        metavar="FILE")

    args = parser.parse_args()
    return args


#####################################################
# Extract annotations
#####################################################

def extract_annotations(filename):
    annotations_order_marker = "VariantDataManager - Annotation order is"
    with open(filename) as annotations_file:
        lines = annotations_file.readlines()
        annotations_order_string = list(filter(lambda a_line: a_line.find(annotations_order_marker) != -1, lines))[0]
        print(annotations_order_string.strip().split(":")[-1].replace())
        return annotations_order_string



#####################################################
# Main
#####################################################


def main(args=None):
    args = parse_args(args)

    extract_annotations(args.annotations_file)


#####################################################
#####################################################

if __name__ == "__main__":
    sys.exit(main())

#####################################################
#####################################################
