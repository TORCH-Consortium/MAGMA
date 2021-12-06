#!/usr/bin/env python3

import os
import sys
import json
import logging
import argparse

logger = logging.getLogger()

def parse_args(args=None):
    Description = "Extract and reshape the annotations for variant recalibration"
    Epilog = "Example usage: reshape_annotations.py annotations.log"

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("annotations_file", type=str, help="The path to the annotations file")
    return parser.parse_args(args)

def main(args=None):
    args = parse_args(args)
    logging.basicConfig(level=args.log_level, format="[%(levelname)s] %(message)s")


    logger.error(f"Annotations file '{annotations_file}' does not exist. Please check the path provided.")

#####################################################
#####################################################

if __name__ == "__main__":
    sys.exit(main())

#####################################################
#####################################################
