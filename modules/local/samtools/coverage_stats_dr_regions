process SAMTOOLS_COVERAGE_STATS_DR_REGIONS {

    tag "$sampleName"

    label 'process_low'

    input:
    tuple val(sampleName), path(bam), path(bai)
    path regions

    output:
    tuple val(sampleName), path("${sampleName}.tier1_tier2_depth.tsv"), emit: depth_tsv
    path "versions.yml", emit: versions

    script:
    """
    samtools bedcov ${regions} ${bam} > ${sampleName}.tier1_tier2_depth.long.tsv

    awk -v sample="${sampleName}" '
    BEGIN {
        OFS="\\t"
        header="sample"
        values=sample
    }
    {
        region=\$1 ":" \$2 "-" \$3
        len=\$3-\$2
        mean_depth=(len > 0 ? \$NF / len : "NA")

        header=header OFS region
        values=values OFS mean_depth
    }
    END {
        print header
        print values
    }
    ' ${sampleName}.tier1_tier2_depth.long.tsv > ${sampleName}.tier1_tier2_depth.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | sed -n '1s/samtools //p')
    END_VERSIONS
    """
}
