process SAMTOOLS_COVERAGE_STATS_DR_REGIONS {
    tag "${sampleName}"

    input:
        tuple val(sampleName), path(bam), path(bai)
        path regions

    output:
        tuple val(sampleName), path("${sampleName}.dr_region_coverage.tsv"), emit: coverage

    script:
        """
        awk -F '[:-]' '
        BEGIN {
            OFS="\\t"
        }
        {
            # Input format:
            # chrom:start-end
            #
            # Convert to BED:
            # chrom  start-1  end
            print \$1, \$2 - 1, \$3
        }
        ' ${regions} > dr_regions.bed

        samtools bedcov dr_regions.bed ${bam} > ${sampleName}.dr_region_coverage.long.tsv

        awk -v sample="${sampleName}" '
        BEGIN {
            OFS="\\t"
            header="sample"
            values=sample
        }
        {
            chrom=\$1
            start=\$2
            end=\$3
            summed_depth=\$NF
            length=end-start

            # Convert BED start back to 1-based display coordinate
            display_start=start + 1

            region=chrom "_" display_start "_" end
            mean_depth=(length > 0 ? summed_depth / length : "NA")

            header=header OFS "dr_region_" region "_mean_depth"
            values=values OFS mean_depth
        }
        END {
            print header
            print values
        }
        ' ${sampleName}.dr_region_coverage.long.tsv > ${sampleName}.dr_region_coverage.tsv
        """
}
