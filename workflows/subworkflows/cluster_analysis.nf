include { CLUSTERPICKER } from "../../modules/clusterpicker/clusterpicker.nf" addParams ( params.CLUSTERPICKER )

workflow CLUSTER_ANALYSIS {

    take:
        incComplex_ch //tuple path(incComplexFasta), path(incComplexNewickTree)
        exComplex_ch //tuple path(exComplexFasta), path(exComplexNewickTree)


    main:
        CLUSTERPICKER(incComplex_ch, 5)
        // CLUSTERPICKER(incComplex_ch, 12)
        // CLUSTERPICKER(exComplex_ch, 5)
        // CLUSTERPICKER(exComplex_ch, 12)
}
