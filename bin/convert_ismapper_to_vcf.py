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

import argparse
import csv
import glob
import os
from Bio import SeqIO

# Define the VCF header
vcf_header = """##fileformat=VCFv4.2
##source=ISMapper
##FILTER=<ID=PASS,Description="All filters passed">
##FILTER=<ID=LowQual,Description="Poor quality and insufficient number of PEs and SRs.">
##ALT=<ID=DEL,Description="Deletion">
##ALT=<ID=DUP,Description="Duplication">
##ALT=<ID=INV,Description="Inversion">
##ALT=<ID=BND,Description="Breakend">
##ALT=<ID=INS,Description="Insertion">
##INFO=<ID=CIEND,Number=2,Type=Integer,Description="PE confidence interval around END">
##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="PE confidence interval around POS">
##INFO=<ID=CHR2,Number=1,Type=String,Description="Chromosome for POS2 coordinate in case of an inter-chromosomal translocation">
##INFO=<ID=POS2,Number=1,Type=Integer,Description="Genomic position for CHR2 in case of an inter-chromosomal translocation">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the structural variant">
##INFO=<ID=PE,Number=1,Type=Integer,Description="Paired-end support of the structural variant">
##INFO=<ID=MAPQ,Number=1,Type=Integer,Description="Median mapping quality of paired-ends">
##INFO=<ID=SRMAPQ,Number=1,Type=Integer,Description="Median mapping quality of split-reads">
##INFO=<ID=SR,Number=1,Type=Integer,Description="Split-read support">
##INFO=<ID=SRQ,Number=1,Type=Float,Description="Split-read consensus alignment quality">
##INFO=<ID=CONSENSUS,Number=1,Type=String,Description="Split-read consensus sequence">
##INFO=<ID=CE,Number=1,Type=Float,Description="Consensus sequence entropy">
##INFO=<ID=CT,Number=1,Type=String,Description="Paired-end signature induced connection type">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Insertion length for SVTYPE=INS">
##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">
##INFO=<ID=PRECISE,Number=0,Type=Flag,Description="Precise structural variation">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=SVMETHOD,Number=1,Type=String,Description="Type of approach used to detect SV">
##INFO=<ID=INSLEN,Number=1,Type=Integer,Description="Predicted length of the insertion">
##INFO=<ID=HOMLEN,Number=1,Type=Integer,Description="Predicted microhomology length using a max. edit distance of 2">
##INFO=<ID=Orientation,Number=1,Type=String,Description="Orientation of the transposable element">
##INFO=<ID=Gap,Number=1,Type=Integer,Description="Gap between sequences">
##INFO=<ID=Call,Number=1,Type=String,Description="Call information">
##INFO=<ID=Percent_ID,Number=1,Type=Float,Description="Percent identity">
##INFO=<ID=Percent_cov,Number=1,Type=Float,Description="Percent coverage">
##INFO=<ID=Left_gene,Number=1,Type=String,Description="Left gene">
##INFO=<ID=Left_description,Number=1,Type=String,Description="Description of the left gene">
##INFO=<ID=Left_strand,Number=1,Type=String,Description="Strand of the left gene">
##INFO=<ID=Left_distance,Number=1,Type=Integer,Description="Distance to the left gene">
##INFO=<ID=Right_gene,Number=1,Type=String,Description="Right gene">
##INFO=<ID=Right_description,Number=1,Type=String,Description="Description of the right gene">
##INFO=<ID=Right_strand,Number=1,Type=String,Description="Strand of the right gene">
##INFO=<ID=Right_distance,Number=1,Type=Integer,Description="Distance to the right gene">
##INFO=<ID=Gene_interruption,Number=1,Type=String,Description="Gene interruption status">
##INFO=<ID=Insertion,Number=1,Type=String,Description="Details of the insertion">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GL,Number=G,Type=Float,Description="Log10-scaled genotype likelihoods for RR,RA,AA genotypes">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=FT,Number=1,Type=String,Description="Per-sample genotype filter">
##FORMAT=<ID=RC,Number=1,Type=Integer,Description="Raw high-quality read counts or base counts for the SV">
##FORMAT=<ID=RCL,Number=1,Type=Integer,Description="Raw high-quality read counts or base counts for the left control region">
##FORMAT=<ID=RCR,Number=1,Type=Integer,Description="Raw high-quality read counts or base counts for the right control region">
##FORMAT=<ID=RDCN,Number=1,Type=Integer,Description="Read-depth based copy-number estimate for autosomal sites">
##FORMAT=<ID=DR,Number=1,Type=Integer,Description="# high-quality reference pairs">
##FORMAT=<ID=DV,Number=1,Type=Integer,Description="# high-quality variant pairs">
##FORMAT=<ID=RR,Number=1,Type=Integer,Description="# high-quality reference junction reads">
##FORMAT=<ID=RV,Number=1,Type=Integer,Description="# high-quality variant junction reads">
##reference=NC-000962-3-H37Rv.fa
##contig=<ID=NC-000962-3-H37Rv,length=4411532>
"""


def get_base_name(file_path):
    base_name = os.path.basename(file_path)
    if '_' in base_name:
        split_base_name = base_name.split("_")
        return f"{split_base_name[0]}_{split_base_name[1]}"
    return base_name


# Function to extract sample name from the ISMapper file path
def extract_sample_name(file_path):
    extracted_sample_name = get_base_name(file_path)  # Extract the sample name (e.g., S011)
    sample_name = f"IS6110.{extracted_sample_name}"  # Prepend 'IS6110.'
    return sample_name


# Function to read the reference genome
def read_reference_genome(reference_file):
    reference_sequences = {}
    reference_sequence_length = 0
    with open(reference_file, 'r') as ref_file:
        # Read the reference genome using Biopython
        for record in SeqIO.parse(ref_file, "fasta"):
            reference_sequences[record.id] = record.seq
            reference_sequence_length = len(record.seq)
    return reference_sequences, reference_sequence_length

# Function to read the transposable element information
def read_transposable_elements(query_file):
    te_info = {}
    with open(query_file, 'r') as query_file:
        # Read the multifasta using Biopython
        for record in SeqIO.parse(query_file, "fasta"):
            te_name = record.id
            te_info[te_name] = len(record.seq)

    return te_info

# Function to convert ISMapper data to VCF format
def convert_is_mapper_to_vcf(is_mapper_dir, vcf_file, reference_sequences, te_info):

    is_mapper_txt_file = glob.glob(is_mapper_dir + "*.txt")[0]
    sample_name = extract_sample_name(is_mapper_txt_file)

    with open(is_mapper_txt_file, 'r') as infile_is_mapper, open(vcf_file, 'w') as outfile:
        # Write the VCF header with the correct sample name
        outfile.write(vcf_header + f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample_name}\n")

        # Read the ISMapper TSV file
        reader_is_mapper = csv.DictReader(infile_is_mapper, delimiter='\t')

        for row in reader_is_mapper:
            chrom = "NC-000962-3-H37Rv"  # Using the specified chromosome
            pos = int(row['x'])  # Convert position to integer
            end = int(row['y'])  # Convert end position to integer
            region_id = row['region']
            ref = reference_sequences[chrom][pos - 1]  # Extract the reference allele from the reference sequence
            orientation = row['orientation']
            te_name = 'IS6110'
            te_length = te_info.get(te_name, 'NA')
            qual = '1000'  # Set QUAL to 1000
            filter_status = 'PASS'
            info = (
                f"SVTYPE=INS;"
                f"Orientation={orientation};"
                f"Gap={row['gap']};"
                f"Call={row['call']};"
                f"Percent_ID={row['percent_ID']};"
                f"Percent_cov={row['percent_cov']};"
                f"Left_gene={row['left_gene']};"
                f"Left_description={row['left_description']};"
                f"Left_strand={row['left_strand']};"
                f"Left_distance={row['left_distance']};"
                f"Right_gene={row['right_gene']};"
                f"Right_description={row['right_description']};"
                f"Right_strand={row['right_strand']};"
                f"Right_distance={row['right_distance']};"
                f"Gene_interruption={row['gene_interruption']};"
                f"SVLEN={te_length};"
                f"END={end};"
                f"Insertion=<{te_name},{orientation}>:{te_length}"
            )
            format_field = "GT:GL:GQ:FT:RCL:RC:RCR:RDCN:DR:DV:RR:RV"
            gl_value = "-665.375,-56.2674,0" # Example values for GL
            gq_value = "10000"     # Genotype Quality
            ft_value = "PASS"      # Per-sample genotype filter
            rcl_value = "60468"    # Raw high-quality read counts or base counts for the left control region
            rc_value = "122933"    # Raw high-quality read counts or base counts for the SV
            rcr_value = "62465"    # Raw high-quality read counts or base counts for the right control region
            rdcn_value = "2"       # Read-depth based copy-number estimate for autosomal sites
            dr_value = "0"         # High-quality reference pairs
            dv_value = "0"         # High-quality variant pairs
            rr_value = "0"         # High-quality reference junction reads
            rv_value = "187"       # High-quality variant junction reads
            gt_value = "1/1"       # Genotype for homozygous alternate allele

            sample_field = f"{gt_value}:{gl_value}:{gq_value}:{ft_value}:{rcl_value}:{rc_value}:{rcr_value}:{rdcn_value}:{dr_value}:{dv_value}:{rr_value}:{rv_value}"

            # Write the VCF entry
            vcf_entry = f"{chrom}\t{pos}\t{region_id}\t{ref}\t<INS>\t{qual}\t{filter_status}\t{info}\t{format_field}\t{sample_field}\n"
            outfile.write(vcf_entry)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Convert ISMapper output to VCF format.')
    parser.add_argument('--is_mapper_dir', required=True, help='Path to the ISMapper output directory.')
    parser.add_argument('--reference_file', required=True, help='Path to the reference genome file.')
    parser.add_argument('--query_file', required=True, help='Path to the transposable elements multifasta file.')
    parser.add_argument('--output_vcf_file', required=True, help='Path to the output VCF file.')

    # Step 3: Parse the command-line arguments
    args = parser.parse_args()

    # Step 4: Replace hardcoded file paths with variables
    is_mapper_dir = args.is_mapper_dir
    vcf_file = args.output_vcf_file
    reference_file = args.reference_file
    query_file = args.query_file

    # Read the reference genome sequence
    reference_sequences, _ = read_reference_genome(reference_file)

    # Read the transposable element information
    query_info = read_transposable_elements(query_file)

    # Run the conversion function
    convert_is_mapper_to_vcf(is_mapper_dir, vcf_file, reference_sequences, query_info)
