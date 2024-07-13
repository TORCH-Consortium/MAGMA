import csv
from Bio import SeqIO

# Define the input files
is_mapper_file = 'S011__NC_000962.3_table.txt'
vcf_file = 'S011_ismapper.vcf'
reference_file = 'NC-000962-3-H37Rv.fa'
te_file = 'transposable_elements.csv'

# Define the VCF header
vcf_header = """##fileformat=VCFv4.2
##source=ISMapper
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
"""

# Function to read the reference genome
def read_reference_genome(reference_file):
    reference_sequences = {}
    with open(reference_file, 'r') as ref_file:
        # Read the reference genome using Biopython
        for record in SeqIO.parse(ref_file, "fasta"):
            reference_sequences[record.id] = record.seq
    return reference_sequences

# Function to read the transposable element information
def read_transposable_elements(te_file):
    te_info = {}
    with open(te_file, 'r') as te_csv:
        reader = csv.DictReader(te_csv)
        for row in reader:
            te_name = row['name']
            te_length = int(row['length'])
            te_info[te_name] = te_length
    return te_info

# Function to convert ISMapper data to VCF format
def convert_is_mapper_to_vcf(is_mapper_file, vcf_file, reference_sequences, te_info):
    with open(is_mapper_file, 'r') as infile, open(vcf_file, 'w') as outfile:
        # Write the VCF header
        outfile.write(vcf_header)

        # Read the ISMapper TSV file
        reader = csv.DictReader(infile, delimiter='\t')

        for row in reader:
            chrom = "NC-000962-3-H37Rv"  # Using the specified chromosome
            pos = int(row['x'])  # Convert position to integer
            region_id = row['region']
            ref = reference_sequences[chrom][pos - 1]  # Extract the reference allele from the reference sequence
            orientation = row['orientation']
            te_name = row['transposable_element']
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

            # Write the VCF entry
            vcf_entry = f"{chrom}\t{pos}\t{region_id}\t{ref}\t{alt}\t{qual}\t{filter_status}\t{info}\n"
            outfile.write(vcf_entry)

# Read the reference genome sequence
reference_sequences = read_reference_genome(reference_file)

# Read the transposable element information
te_info = read_transposable_elements(te_file)

# Run the conversion function
convert_is_mapper_to_vcf(is_mapper_file, vcf_file, reference_sequences, te_info)
