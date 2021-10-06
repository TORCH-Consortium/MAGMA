
RG="@RG\tID:$FLOWCELL.$LANE\tSM:$STUDY.$SAMPLE\tPL:illumina\tLB:lib$LIBRARY\tPU:$FLOWCELL.$LANE.$INDEX"

$BWA mem -M -t $BWA_THREADS -R $RG $REFERENCE $SEQ_R1 $SEQ_R2 | $SAMTOOLS sort -@ $SAMTOOLS_THREADS -O BAM -o $OUT_DIR/mapped_singles/$UNIQUE_ID.sorted_reads.bam -



nextflow.enable.dsl = 2

params.results_dir = "${params.outdir}/tbbwa"
params.save_mode = 'copy'
params.should_publish = true



process process_name{
    tag "something"
    publishdir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
    path(somefile)

    output:
    path("pattern"),  emit: "ch_output"

    script:

    """
    echo "nothing"
    """

    stub:

    """
    echo "nothing on stub"
    """

}
