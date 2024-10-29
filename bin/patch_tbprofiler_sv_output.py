#!/usr/bin/env python3


import argparse
import subprocess
import pandas as pd
import pysam
import os
import glob
import json

# def concatenate_vcf_files(delly_vcf, ismapper_vcf, output_vcf):
#     bcftools_command = [

#         'bcftools', 'concat', '-o', output_vcf, '-O', 'v', delly_vcf, ismapper_vcf
#     ]
#     subprocess.run(bcftools_command, check=True)
#     print(f'Concatenated VCF file created at: {output_vcf}')



def extract_unique_prefixes(file_list):
  """Extracts unique prefixes from a list of file paths.

  Args:
    file_list: A list of file paths.

  Returns:
    A list of unique prefixes.
  """

  prefixes = set()
  for file_path in file_list:
    # Split the path based on "." (dot)
    parts = file_path.split(".")

    # Check if there are at least 3 parts (NICD, prefix, extension)
    if len(parts) >= 3:
      # The prefix is the second part (excluding the leading "NICD")
      prefix = f"{parts[0]}.{parts[1]}"
      prefixes.add(prefix)
  return list(prefixes)  # Convert set back to a list

# json_files = glob.glob(args['indir'] + "/" + "*.json")

def load_bed_file(bed_file):
    bed_df = pd.read_csv(bed_file, sep='\t', header=None, names=['Chromosome', 'start', 'end', 'gene_code', 'gene_name', 'drug'])
    bed_intervals = {}
    for _, row in bed_df.iterrows():
        chrom = row['Chromosome']
        if chrom not in bed_intervals:
            bed_intervals[chrom] = []
        bed_intervals[chrom].append(row)
    return bed_intervals

def process_vcf_file(output_vcf, bed_intervals):
    vcf_in = pysam.VariantFile(output_vcf)
    json_output = []

    for record in vcf_in:
        vcf_chrom = 'Chromosome'
        vcf_start = record.start
        vcf_end = record.stop

        if vcf_chrom in bed_intervals:
            for bed_entry in bed_intervals[vcf_chrom]:
                bed_start = bed_entry['start']
                bed_end = bed_entry['end']

                overlap = (
                    (vcf_start <= bed_end and vcf_end >= bed_start) or
                    (bed_start <= vcf_end and bed_end >= vcf_start)
                )

                if overlap:
                    sv_type = record.info.get('SVTYPE', None)
                    if sv_type in ['DEL', 'DUP', 'INV']:
                        vcf_svlen = vcf_end - vcf_start + 1
                    elif sv_type == 'INS':
                        vcf_svlen = abs(record.info.get('SVLEN', record.info.get('INSLEN', 0)))
                    else:
                        vcf_svlen = None

                    change_format = sv_type.lower() if sv_type in ['INS', 'DEL', 'DUP', 'INV'] else 'unknown'

                    json_record = {
                        "chrom": vcf_chrom,
                        "pos": vcf_start + 1,
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
                                        "original_mutation": "TE_insertion",
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

                    json_output.append(json_record)
                  
    return json_output

def update_json(processed_bed_file, existing_json_file, concat_vcf, output_json_file):
    json_output= process_vcf_file(concat_vcf, processed_bed_file)

    try:
        with open(existing_json_file, 'r') as file:
            existing_data = json.load(file)
    except IOError as e:
        print(f"Error reading file {existing_json_file}: {e}")
        existing_data = {}

    for record in existing_data.get('dr_variants', []):
        record['source'] = 'TBprofiler'

    for record in json_output:
        record['source'] = 'Delly/Ismapper'

    combined_data = existing_data.get('dr_variants', []) + json_output

    unique_records = {}
    for record in combined_data:
        key = (record["chrom"], record["pos"], record["alt"], record["sv_len"], record["gene_id"])
        if key in unique_records:
            if unique_records[key]['source'] == 'TBprofiler':
                continue
        unique_records[key] = record

    existing_data["dr_variants"] = [v for k, v in unique_records.items()]

    results_directory = 'results'
    os.makedirs(results_directory, exist_ok=True)
    output_json_file = os.path.join(results_directory, output_json_file)

    try:
        with open(output_json_file, 'w') as file:
            json.dump(existing_data, file, indent=4)
        print(f"Cleaned JSON file saved: {output_json_file}")
    except IOError as e:
        print(f"Error writing file {output_json_file}: {e}")


def main():
    parser = argparse.ArgumentParser(description='Process VCF and JSON files.')
    # parser.add_argument('--delly_vcf', required=True, help='Path to the Delly VCF file')
    # parser.add_argument('--ismapper_vcf', required=True, help='Path to the ISMapper VCF file')
    parser.add_argument('--bed_file', required=True, help='Path to the BED file')
    parser.add_argument('--existing_json_dir', required=True, help='Path to the existing TBprofiler JSON output files')
    parser.add_argument('--concat_vcf_dir', required=True, help='Path to the concatenated VCF files')

    args = parser.parse_args()

    unique_prefixes = extract_unique_prefixes(os.listdir(args.existing_json_dir))

    processed_bed_file = load_bed_file(args.bed_file)

    for p in unique_prefixes:
        concat_vcf_file = glob.glob(args.concat_vcf_dir + "/" + p  + "*")[0]
        existing_json_file = glob.glob(args.existing_json_dir + "/" + p  + "*")[0]
        output_json_file = f"{p}.patched.json"

        print(existing_json_file, concat_vcf_file, output_json_file)

        update_json( processed_bed_file , existing_json_file, concat_vcf_file, output_json_file)

if __name__ == "__main__":
    main()
