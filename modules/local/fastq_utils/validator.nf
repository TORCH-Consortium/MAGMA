/*
 * Copyright (c) 2021-2024 MAGMA pipeline authors, see https://doi.org/10.1371/journal.pcbi.1011648
 *
 * This file is part of MAGMA pipeline, see https://github.com/TORCH-Consortium/MAGMA
 *
 * For quick overview of GPL-3 license, please refer
 * https://www.tldrlegal.com/license/gnu-general-public-license-v3-gpl-3
 *
 * - You MUST keep this license with original authors in your copy
 * - You MUST acknowledge the original source of this software
 * - You MUST state significant changes made to the original software
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program . If not, see <http://www.gnu.org/licenses/>.
 */
process FASTQ_VALIDATOR {
    tag "${sampleName}"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish


    input:
        tuple val(sampleName), path(sampleRead)
        val ready

    output:
        path("*.fastq_report.csv")                                       , emit: fastq_report
        tuple val(sampleName), path(sampleRead)                          , emit: reads

    shell:

        '''
        md5sum !{sampleRead} > !{sampleRead.simpleName}.md5sum_out.txt
        cat *md5sum_out.txt | csvtk space2tab | csvtk tab2csv | csvtk add-header -n md5sum,file > !{sampleRead.simpleName}.md5sum_stats.csv

        du -shL !{sampleRead} > !{sampleRead.simpleName}.du_out.txt
        cat *du_out.txt | csvtk tab2csv | csvtk add-header -n size,file > !{sampleRead.simpleName}.du_stats.csv



        csvtk join -f file \\
        !{sampleRead.simpleName}.md5sum_stats.csv \\
        !{sampleRead.simpleName}.du_stats.csv \\
        > !{sampleRead.simpleName}.du_md5sum_stats.csv

        rm *_out.txt !{sampleRead.simpleName}.du_stats.csv !{sampleRead.simpleName}.md5sum_stats.csv


        seqkit stats -a -T  !{sampleRead}  > !{sampleRead.simpleName}.seqkit_out.txt
        cat *seqkit_out.txt | csvtk space2tab | csvtk tab2csv > !{sampleRead.simpleName}.seqkit_stats.csv


        csvtk join -f file \\
        !{sampleRead.simpleName}.seqkit_stats.csv \\
        !{sampleRead.simpleName}.du_md5sum_stats.csv \\
        > !{sampleRead.simpleName}.fastq_statistics.csv

        rm *_out.txt


        !{params.fastq_validator_path} !{sampleRead} \\
        2>!{sampleRead.simpleName}.command.log || true

        cp !{sampleRead.simpleName}.command.log .command.log



        TEMP=$(tail -n 1 !{sampleRead.simpleName}.command.log)


        if [ "$(echo "$TEMP")" == "OK" ]; then
            VALIDATED=1
            STATUS="passed"
            echo -e "file,magma_name,fastq_utils_check" > !{sampleRead.simpleName}.check.${STATUS}.csv
            echo -e "!{sampleRead},!{sampleName},${STATUS}" >> !{sampleRead.simpleName}.check.${STATUS}.csv

            csvtk join -f file  !{sampleRead.simpleName}.fastq_statistics.csv !{sampleRead.simpleName}.check.${STATUS}.csv >  !{sampleRead.name}.fastq_report.csv



            rm !{sampleRead.simpleName}.fastq_statistics.csv
            exit 0

        else
            VALIDATED=0
            STATUS="failed"
            echo -e "file,magma_name,fastq_utils_check" > !{sampleRead.simpleName}.check.${STATUS}.csv
            echo -e "!{sampleRead},!{sampleName},${STATUS}" >> !{sampleRead.simpleName}.check.${STATUS}.csv

            csvtk join -f file  !{sampleRead.simpleName}.fastq_statistics.csv !{sampleRead.simpleName}.check.${STATUS}.csv  > !{sampleRead.name}.fastq_report.csv


            rm !{sampleRead.simpleName}.fastq_statistics.csv
            exit 1
        fi


        '''

    stub:

        """
        touch ${sampleRead.simpleName}.check.csv
        """

}
