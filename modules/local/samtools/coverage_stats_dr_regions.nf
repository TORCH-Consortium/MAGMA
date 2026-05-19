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
            chrom=\\$1
            bed_start=\\$2
            bed_stop=\\$3
            summed_depth=\\$NF
            region_size=bed_stop-bed_start

            display_start=bed_start + 1
            region_name=chrom "_" display_start "_" bed_stop
            mean_depth=(region_size > 0 ? summed_depth / region_size : "NA")

            header=header OFS "dr_region_" region_name "_mean_depth"
            values=values OFS mean_depth
        }
        END {
            print header
            print values
        }
        ' ${sampleName}.dr_region_coverage.long.tsv > ${sampleName}.dr_region_coverage.tsv
        """
}
