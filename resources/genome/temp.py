# from Bio import SeqIO
# record = SeqIO.read("NC_002945v4.gb", "genbank")
# for feature in record.features:
#     if feature.type == "rRNA":
#         print(f"{feature.qualifiers['product'][0]}: {feature.location}")



#!/usr/bin/env python3
"""
Extract rRNA gene coordinates from Mycobacterium bovis AF2122/97 genome
Reads from local GenBank file: NC_002945v4.gb
Outputs coordinates in BED format suitable for MAGMA pipeline
"""

from Bio import SeqIO
import sys
import os

def extract_rrna_coordinates(record):
    """Extract rRNA gene coordinates from GenBank record"""
    rrna_features = []

    for feature in record.features:
        if feature.type == "rRNA":
            # Extract gene information
            product = feature.qualifiers.get('product', ['Unknown'])[0]
            gene = feature.qualifiers.get('gene', [''])[0]
            locus_tag = feature.qualifiers.get('locus_tag', [''])[0]

            # Get coordinates (convert to 0-based for BED format)
            start = int(feature.location.start)  # Already 0-based in BioPython
            end = int(feature.location.end)
            strand = '+' if feature.location.strand == 1 else '-'

            rrna_features.append({
                'gene': gene,
                'product': product,
                'locus_tag': locus_tag,
                'start': start,
                'end': end,
                'strand': strand,
                'length': end - start
            })

    return rrna_features

def write_bed_file(rrna_features, output_file, genome_name):
    """Write coordinates to BED format file"""
    with open(output_file, 'w') as f:
        # Write header (optional, comment out if not needed)
        f.write(f"# rRNA regions for {genome_name}\n")
        f.write(f"# Format: chromosome start end name score strand\n")

        for feature in sorted(rrna_features, key=lambda x: x['start']):
            # BED format: chrom start end name score strand
            f.write(f"{genome_name}\t{feature['start']}\t{feature['end']}\t"
                   f"{feature['gene']}\t0\t{feature['strand']}\n")

    print(f"\nWrote BED file: {output_file}")

def write_simple_list(rrna_features, output_file):
    """Write coordinates to simple list format (similar to MAGMA)"""
    with open(output_file, 'w') as f:
        for feature in sorted(rrna_features, key=lambda x: x['start']):
            # Format: start-end (or gene:start-end)
            f.write(f"{feature['gene']}:{feature['start']}-{feature['end']}\n")

    print(f"Wrote list file: {output_file}")

def write_magma_format(rrna_features, output_file):
    """Write coordinates in MAGMA rRNA.list format"""
    with open(output_file, 'w') as f:
        for feature in sorted(rrna_features, key=lambda x: x['start']):
            # MAGMA format appears to be: start-end
            f.write(f"{feature['start']}-{feature['end']}\n")

    print(f"Wrote MAGMA format file: {output_file}")

def main():
    genbank_file = "NC_002945v4.gb"
    genome_name = "NC_002945"  # Chromosome name for BED file

    print("="*60)
    print("Extracting rRNA coordinates from M. bovis AF2122/97")
    print("="*60)

    # Check if file exists
    if not os.path.exists(genbank_file):
        print(f"\nError: File '{genbank_file}' not found!")
        print("Please ensure the GenBank file is in the current directory.")
        sys.exit(1)

    # Read GenBank file
    print(f"\nReading {genbank_file}...")
    try:
        record = SeqIO.read(genbank_file, "genbank")
    except Exception as e:
        print(f"Error reading GenBank file: {e}")
        sys.exit(1)

    print(f"\nGenome: {record.description}")
    print(f"Length: {len(record.seq):,} bp")
    print(f"Accession: {record.id}")

    # Extract rRNA coordinates
    rrna_features = extract_rrna_coordinates(record)

    if not rrna_features:
        print("\nWarning: No rRNA features found!")
        sys.exit(1)

    # Display results
    print(f"\nFound {len(rrna_features)} rRNA genes:")
    print("-"*80)
    print(f"{'Gene':<10} {'Product':<20} {'Start':<10} {'End':<10} {'Length':<8} {'Strand'}")
    print("-"*80)

    for feature in sorted(rrna_features, key=lambda x: x['start']):
        print(f"{feature['gene']:<10} {feature['product']:<20} "
              f"{feature['start']:<10} {feature['end']:<10} "
              f"{feature['length']:<8} {feature['strand']}")

    # Write output files
    write_bed_file(rrna_features, "Mbovis_rRNA.bed", genome_name)
    write_simple_list(rrna_features, "Mbovis_rRNA.list")
    write_magma_format(rrna_features, "rRNA.list")  # MAGMA-compatible format

    print("\n" + "="*60)
    print("Complete! Files created:")
    print("  - Mbovis_rRNA.bed     (BED format)")
    print("  - Mbovis_rRNA.list    (Simple list with gene names)")
    print("  - rRNA.list           (MAGMA format - ready to use!)")
    print("="*60)
    print("\nThe 'rRNA.list' file is ready to use with MAGMA pipeline")
    print("for M. bovis AF2122/97")

if __name__ == "__main__":
    main()
