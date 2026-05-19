process SAMTOOLS_COVERAGE_STATS_DR_REGIONS {
    tag "${sampleName}"

    input:
        tuple val(sampleName), path(bam), path(bai)
        path regions

    output:
        tuple val(sampleName), path("${sampleName}.dr_region_coverage.tsv"), emit: coverage

    script:
        """
        samtools bedcov ${regions} ${bam} > ${sampleName}.dr_region_coverage.long.tsv

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

            region=chrom "_" start "_" end
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
