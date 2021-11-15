include { CLUSTERPICKER as CLUSTERPICKER__5SNP  } from "../../modules/clusterpicker/clusterpicker.nf" addParams ( params.CLUSTERPICKER )
include { CLUSTERPICKER as CLUSTERPICKER__12SNP  } from "../../modules/clusterpicker/clusterpicker.nf" addParams ( params.CLUSTERPICKER )

workflow CLUSTER_ANALYSIS {

    take:
        cluster_files_ch
        prefix


    main:
        CLUSTERPICKER__5SNP(cluster_files_ch, 5, prefix)
        CLUSTERPICKER__12SNP(cluster_files_ch, 12, prefix)

}
