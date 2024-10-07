#!/usr/bin/env python3

import os
import subprocess
import shlex

# Create temporary directory
tmp_dir = "interim"
os.makedirs(tmp_dir, exist_ok=True)

# Extract unique sample prefixes from filenames
prefixes = [f"{os.path.basename(file).split('.')[0]}.{os.path.basename(file).split('.')[1]}" for file in os.listdir('.') if file.endswith('.bcf.gz') or file.endswith('.vcf.gz')]


prefixes = list(set(prefixes))

# Process each sample prefix
for prefix in prefixes:
    # Find related files for the prefix
    files = [f for f in os.listdir('.') if f.startswith(prefix) and (f.endswith('.bcf.gz') or f.endswith('.vcf.gz'))]

    if files:
        # Concatenate files
        concat_file = os.path.join(tmp_dir, f"{prefix}.concat.vcf.gz")
        subprocess.run(shlex.split(f"bcftools concat {' '.join(files)} -o {concat_file}"))
        subprocess.run(shlex.split(f"bcftools view {concat_file} -Oz -o {concat_file}"))

        # Sort and index the concatenated file
        sorted_file = os.path.join(tmp_dir, f"{prefix}.concat.sorted.vcf.gz")
        subprocess.run(shlex.split(f"bcftools sort {concat_file} -Oz -o {sorted_file}"))
        subprocess.run(shlex.split(f"bcftools index {sorted_file}"))

        # Add sorted file to merged list
        concat_files.append(sorted_file)

# Merge concatenated files (if any)
if concat_files:
    subprocess.run(shlex.split(f"bcftools merge {' '.join(concat_files)} -o joint.delly.vcf"))
    subprocess.run(shlex.split("bgzip joint.delly.vcf"))
    subprocess.run(shlex.split("bcftools index joint.delly.vcf.gz"))

# # NOTE For now commenting this for debugging.
# # Clean up temporary directory
# os.rmdir(tmp_dir)
