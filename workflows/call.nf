nextflow.enable.dsl = 2

include { GATK_HAPLOTYPE_CALLER } from "../../modules/gatk/haplotype_caller.nf" addParams(params.GATK_HAPLOTYPE_CALLER)

workflow CALL {

}
