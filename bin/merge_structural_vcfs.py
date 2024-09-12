import argparse
import subprocess
import pandas as pd
import pysam
import json

def concatenate_vcf_files(delly_vcf, ismapper_vcf, output_vcf):
    bcftools_command = [
        'bcftools', 'concat', '-o', output_vcf, '-O', 'v', delly_vcf, ismapper_vcf
    ]
    subprocess.run(bcftools_command, check=True)
    print(f'Concatenated VCF file created at: {output_vcf}')

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
        vcf_chrom = record.chrom
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

                    json_output.append(json_record)
                    break
    return json_output

def update_json(existing_json_file, json_output, cleaned_json_file):
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

    try:
        with open(cleaned_json_file, 'w') as file:
            json.dump(existing_data, file, indent=4)
        print(f"Cleaned JSON file saved: {cleaned_json_file}")
    except IOError as e:
        print(f"Error writing file {cleaned_json_file}: {e}")

def main():
    parser = argparse.ArgumentParser(description='Process VCF and JSON files.')
    parser.add_argument('--delly_vcf', required=True, help='Path to the Delly VCF file')
    parser.add_argument('--ismapper_vcf', required=True, help='Path to the ISMapper VCF file')
    parser.add_argument('--bed_file', required=True, help='Path to the BED file')
    parser.add_argument('--existing_json_file', required=True, help='Path to the existing JSON file')
    parser.add_argument('--cleaned_json_file', required=True, help='Path to the cleaned JSON file')
    parser.add_argument('--output_vcf', required=True, help='Path to the output concatenated VCF file')

    args = parser.parse_args()

    concatenate_vcf_files(args.delly_vcf, args.ismapper_vcf, args.output_vcf)
    bed_intervals = load_bed_file(args.bed_file)
    json_output = process_vcf_file(args.output_vcf, bed_intervals)
    update_json(args.existing_json_file, json_output, args.cleaned_json_file)

if __name__ == "__main__":
    main()