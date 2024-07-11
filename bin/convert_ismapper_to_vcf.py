import csv

# Define the input ISMapper file and output VCF file
is_mapper_file = 'is_mapper_output.tsv'
vcf_file = 'output.vcf'

# Define the VCF header
vcf_header = """##fileformat=VCFv4.2
##source=ISMapper
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO
"""

# Function to convert ISMapper data to VCF format
def convert_is_mapper_to_vcf(is_mapper_file, vcf_file):
    with open(is_mapper_file, 'r') as infile, open(vcf_file, 'w') as outfile:
        # Write the VCF header
        outfile.write(vcf_header)

        # Read the ISMapper TSV file
        reader = csv.DictReader(infile, delimiter='\t')

        for row in reader:
            chrom = "NC-000962-3-H37Rv"  # Using the specified chromosome
            pos = row['x']
            region_id = row['region']
            ref = 'N'
            alt = f"N[<ctg1>:1473["
            qual = '.'
            filter_status = 'PASS'
            info = (
                f"SVTYPE=BND;"
                f"Orientation={row['orientation']};"
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

# Run the conversion function
convert_is_mapper_to_vcf(is_mapper_file, vcf_file)
