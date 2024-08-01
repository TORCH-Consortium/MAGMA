#!/usr/bin/env python3

import argparse
import csv
import glob
from Bio import SeqIO

# Define the VCF header
vcf_header = """##fileformat=VCFv4.2
##source=ISMapper
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
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
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE
"""

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

        # NOTE: Renable in case we need to dump this on disk.
        # with open("temp_transposable_elements.csv", mode='w', newline='') as file:
        #     writer = csv.writer(file)
        #     writer.writerow(["name", "length"])
        #     for te_name, te_length in te_info.items():
        #         writer.writerow([te_name, te_length])
    return te_info

# Function to convert ISMapper data to VCF format
def convert_is_mapper_to_vcf(is_mapper_dir, vcf_file, reference_sequences, te_info):

    is_mapper_txt_file = glob.glob(is_mapper_dir + "*.txt")[0]

    with open(is_mapper_txt_file, 'r') as infile_is_mapper, open(vcf_file, 'w') as outfile:
        # Write the VCF header
        outfile.write(vcf_header)

        # Read the ISMapper TSV file
        reader_is_mapper = csv.DictReader(infile_is_mapper, delimiter='\t')

        for row in reader_is_mapper:
            chrom = "NC-000962-3-H37Rv"  # Using the specified chromosome
            pos = int(row['x'])  # Convert position to integer
            region_id = row['region']
            ref = reference_sequences[chrom][pos - 1]  # Extract the reference allele from the reference sequence
            orientation = row['orientation']
            # FIXME Hard-code the name of this specific element
            te_name = 'IS6110'
            te_length = te_info.get(te_name, 'NA')
            alt = f"{ref}[<{te_name},{orientation}>:{te_length}["  # Use transposable element, orientation, and its length
            qual = '.'
            filter_status = 'PASS'
            info = (
                f"SVTYPE=BND;"
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
                f"Gene_interruption={row['gene_interruption']}"
            )
            format_field = "GT:AD:DP:GQ:PL"
            ad_value = "10,10"  # Allelic depths for the ref and alt alleles
            dp_value = "10"     # Approximate read depth
            gq_value = "99"     # Genotype Quality
            gt_value = "1"      # Genotype
            pl_value = "1800"   # Normalized, Phred-scaled likelihoods

            sample_field = f"{gt_value}:{ad_value}:{dp_value}:{gq_value}:{pl_value}"

            # Write the VCF entry
            vcf_entry = f"{chrom}\t{pos}\t{region_id}\t{ref}\t{alt}\t{qual}\t{filter_status}\t{info}\t{format_field}\t{sample_field}\n"
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
