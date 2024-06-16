import vcf

from scipy.stats import binom_test

# Path to your VCF file

vcf_file = "your_file.vcf"  # Replace "your_file.vcf" with the actual file name

# Open the VCF file for reading

vcf_reader = vcf.Reader(open(vcf_file, 'r'))

# Create a list to store the filtered records

filtered_records = []

# Iterate over each record (variant) in the VCF file

for record in vcf_reader:

    # Extract the DP4 field from the INFO dictionary

    dp4_values = record.INFO.get('DP4')

    # Perform binomial test on the last two values

    if dp4_values:

        p_value = binom_test([dp4_values[-2], dp4_values[-1]], n=sum(dp4_values[-2:]), p=0.5)

        # Check if p-value is above or equal to 0.05

        if p_value >= 0.05:

            filtered_records.append(record)

# Write the filtered records to a new VCF file

output_vcf = "filtered_output.vcf"

vcf_writer = vcf.Writer(open(output_vcf, 'w'), vcf_reader)


for record in filtered_records:

    vcf_writer.write_record(record)

# Close the VCF writer

vcf_writer.close()


print("Filtered VCF file has been created successfully.")


