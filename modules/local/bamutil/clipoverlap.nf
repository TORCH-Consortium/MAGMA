process BAMUTIL_CLIPOVERLAP {
    tag "${sampleName}"
    label 'cpu_2_memory_2'

    container "community.wave.seqera.io/library/bamutil:1.0.15--1582863f9eda77c6"

    input:
        tuple val(sampleName), path(bam)

    output:
        tuple val(sampleName), path("*.clipped.bam")
        path("*.clipOverlap.log"), emit: stats

    script:
        """
        bam clipOverlap \\
            --in ${bam} \\
            --out ${sampleName}.clipped.bam \\
            --stats \\
            2> ${sampleName}.clipOverlap.log
        """
}
