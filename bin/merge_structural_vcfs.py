import pysam
import pandas as pd
import json
import subprocess  # Import subprocess module

# Define file paths
delly_vcf = '/mnt/c/Users/peppe/OneDrive/Documentos/Escritorio/downloads/int/pipe2/IS6110.S011.filtered.delly.bcf'
ismapper_vcf = '/mnt/c/Users/peppe/OneDrive/Documentos/Escritorio/downloads/int/IS6110.S011.L1.A1.1.1.1.ismapper.vcf'
bed_file = '/mnt/c/Users/peppe/OneDrive/Documentos/Escritorio/downloads/int/pipe2/who_plus.bed'
existing_json_file = '/mnt/c/Users/peppe/OneDrive/Documentos/Escritorio/downloads/int/pipe2/IS6110.S011.results-4_TB_profiler.json'
updated_json_file = '/mnt/c/Users/peppe/OneDrive/Documentos/Escritorio/downloads/int/pipe2/IS6110.S011.results-4_TBprofiler_updated_fromdelly.json'
cleaned_json_file = '/mnt/c/Users/peppe/OneDrive/Documentos/Escritorio/downloads/int/pipe2/IS6110.S011.results-4_TBprofiler_cleaned_fromdelly.json'

# Output concatenated VCF file path
output_vcf = '/mnt/c/Users/peppe/OneDrive/Documentos/Escritorio/downloads/int/concatenated_output.vcf'

# Command to concatenate BCF/VCF files using bcftools
bcftools_command = [
    'bcftools', 'concat', '-o', output_vcf, '-O', 'v', delly_vcf, ismapper_vcf  # Corrected variable name
]

# Run the command
subprocess.run(bcftools_command, check=True)

print(f'Concatenated VCF file created at: {output_vcf}')


# Step 1: Load BED file into a DataFrame
bed_df = pd.read_csv(bed_file, sep='\t', header=None, names=['Chromosome', 'start', 'end', 'gene_code', 'gene_name', 'drug'])

# Create a dictionary to store BED intervals by chromosome
bed_intervals = {}
for _, row in bed_df.iterrows():
    chrom = row['Chromosome']
    if chrom not in bed_intervals:
        bed_intervals[chrom] = []
    bed_intervals[chrom].append(row)

# Step 2: Process the VCF file and create JSON output
vcf_in = pysam.VariantFile(output_vcf)
json_output = []

for record in vcf_in:
    vcf_chrom = record.chrom
    vcf_start = record.start
    vcf_end = record.stop  # Use record.stop to get the end position

    if vcf_chrom in bed_intervals:
        for bed_entry in bed_intervals[vcf_chrom]:
            bed_start = bed_entry['start']
            bed_end = bed_entry['end']

            # Check for overlap between VCF and BED intervals
            overlap = (
                (vcf_start <= bed_end and vcf_end >= bed_start) or
                (bed_start <= vcf_end and bed_end >= vcf_start)
            )
            
            if overlap:
                sv_type = record.info.get('SVTYPE', None)
                
                # Compute sv_len
                if sv_type in ['DEL', 'DUP', 'INV']:
                    vcf_svlen = vcf_end - vcf_start + 1
                elif sv_type == 'INS':
                    vcf_svlen = abs(record.info.get('SVLEN', record.info.get('INSLEN', 0)))
                else:
                    vcf_svlen = None

                # Determine change format
                change_format = sv_type.lower() if sv_type in ['INS', 'DEL', 'DUP', 'INV'] else 'unknown'

                # Create JSON record
                json_record = {
                    "chrom": vcf_chrom,
                    "pos": vcf_start + 1,  # VCF is 0-based, JSON is 1-based
                    "ref": record.ref,
                    "alt": record.alts[0],
                    "depth": None,
                    "freq": None,
                    "sv": True,
                    "filter": list(record.filter.keys())[0] if record.filter else "PASS",
                    "forward_reads": None,
                    "reverse_reads": None,
                    "sv_len": vcf_svlen,
                    "gene_id": bed_entry['gene_code'],
                    "gene_name": bed_entry['gene_name'],
                    "feature_id": None,
                    "type": None,
                    "change": f"c.{vcf_start + 1}_{vcf_end}{change_format}",
                    "nucleotide_change": f"c.{vcf_start + 1}_{vcf_end}{change_format}",
                    "protein_change": None,
                    "annotation": [
                        {
                            "type": "drug_resistance",
                            "drug": bed_entry['drug'],
                            "original_mutation": "LoF",
                            "confidence": " ",
                            "source": "",
                            "comment": ""
                        }
                    ],
                    "consequences": [
                        {
                            "gene_id": bed_entry['gene_code'],
                            "gene_name": bed_entry['gene_name'],
                            "feature_id": None,
                            "type": None,
                            "nucleotide_change": f"c.{vcf_start + 1}_{vcf_end}{change_format}",
                            "protein_change": None,
                            "annotation": [
                                {
                                    "type": "drug_resistance",
                                    "drug": bed_entry['drug'],
                                    "original_mutation": "LoF",
                                    "confidence": "",
                                    "source": "",
                                    "comment": ""
                                }
                            ]
                        }
                    ],
                    "drugs": [
                        {
                            "type": "drug_resistance",
                            "drug": bed_entry['drug'],
                            "original_mutation": "LoF",
                            "confidence": "",
                            "source": "",
                            "comment": ""
                        }
                    ],
                    "locus_tag": bed_entry['gene_code'],
                    "gene_associated_drugs": [bed_entry['drug']] if bed_entry['drug'] else []
                }

                # Append the record to the output list
                json_output.append(json_record)
                break

# Step 3: Update the existing JSON with new data
try:
    with open(existing_json_file, 'r') as file:
        existing_data = json.load(file)
except IOError as e:
    print(f"Error reading file {existing_json_file}: {e}")
    existing_data = {}

# Mark records from TBprofiler as "source: TBprofiler"
for record in existing_data.get('dr_variants', []):
    record['source'] = 'TBprofiler'

# Mark new records from Delly/Ismapper as "source: Delly/Ismapper"
for record in json_output:
    record['source'] = 'Delly/Ismapper'

# Combine both sets of records
combined_data = existing_data.get('dr_variants', []) + json_output

# Step 4: Remove duplicates, prioritizing TBprofiler records
unique_records = {}
for record in combined_data:
    key = (record["chrom"], record["pos"], record["alt"], record["sv_len"], record["gene_id"])
    
    # If the key exists and the existing record is from TBprofiler, skip adding the new one
    if key in unique_records:
        if unique_records[key]['source'] == 'TBprofiler':
            continue
    
    unique_records[key] = record

# Update with unique records, removing the 'source' key
existing_data["dr_variants"] = [v for k, v in unique_records.items()]

# Step 5: Save the cleaned data to a new JSON file
try:
    with open(cleaned_json_file, 'w') as file:
        json.dump(existing_data, file, indent=4)
    print(f"Cleaned JSON file saved: {cleaned_json_file}")
except IOError as e:
    print(f"Error writing file {cleaned_json_file}: {e}")
