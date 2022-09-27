nextflow.enable.dsl = 2

include { MULTIQC } from '../modules/multiqc/multiqc.nf' addParams (params.MULTIQC)

workflow REPORTS_WF {
    take:
         reports_fastqc_ch

    main:
        MULTIQC(reports_fastqc_ch)

}
