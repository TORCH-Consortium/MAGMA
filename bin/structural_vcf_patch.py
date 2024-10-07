#!/usr/bin/env python3

import os
import tempfile
import subprocess
from glob import glob

# Parameters
vcf_files =  glob.glob("*vcf.gz")  # Example list of VCF files
bcftools_path = "bcftools"  # Path to bcftools
vcf_name = "joint"  # Output VCF name
file_format = "gz"  # Output file format

# Create temporary directory
tmp_dir = tempfile.mkdtemp(prefix="bcftools-")

try:
    # Extract sample prefixes
    prefixes = sorted(set(os.path.basename(file).split('.')[1] for file in vcf_files))

    # Concatenate files for each sample prefix
    concat_files = []
    for prefix in prefixes:
        files = glob(f"{prefix}.*.bcf.gz") + glob(f"{prefix}.*.vcf.gz")
        if files:
            concat_file = os.path.join(tmp_dir, f"{prefix}.concat.vcf")
            subprocess.run([bcftools_path, "concat"] + files + ["-o", concat_file], check=True)
            subprocess.run(["bgzip", concat_file], check=True)
            sorted_file = os.path.join(tmp_dir, f"{prefix}.sorted.vcf.gz")
            subprocess.run([bcftools_path, "sort", f"{concat_file}.gz", "-o", sorted_file], check=True)
            subprocess.run([bcftools_path, "index", sorted_file], check=True)
            concat_files.append(sorted_file)

    # Merge the concatenated files
    if concat_files:
        output_vcf = f"{vcf_name}.{file_format}.vcf"
        subprocess.run([bcftools_path, "merge"] + concat_files + ["-o", output_vcf], check=True)
        subprocess.run(["bgzip", output_vcf], check=True)
        subprocess.run([bcftools_path, "index", f"{output_vcf}.gz"], check=True)
finally:
    # Clean up the temporary directory
    for root, dirs, files in os.walk(tmp_dir, topdown=False):
        for name in files:
            os.remove(os.path.join(root, name))
        for name in dirs:
            os.rmdir(os.path.join(root, name))
    os.rmdir(tmp_dir)
