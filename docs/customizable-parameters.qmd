# MAGMA Customizable Parameters


This document provides an overview of the customizable parameters for the MAGMA pipeline. Each parameter is listed with its default value, description.

---

## Common Parameters

### Input Samplesheet
| Parameter             | Default Value              | Description                                                                                     |
|-----------------------|----------------------------|-------------------------------------------------------------------------------------------------|
| `input_samplesheet`   | `"samplesheet.magma.csv"`  | The input CSV file containing sample information. The study ID cannot start with `XBS_REF_`.   |

> 💡 **Hint**: The samplesheet should include the fields `[Sample, R1, R2]`. Optionally, you can add `[study, library, attempt, flowcell, lane, index_sequence]`.

---

### Output Directory
| Parameter   | Default Value         | Description                                                                 |
|-------------|-----------------------|-----------------------------------------------------------------------------|
| `outdir`    | `"magma-results"`     | The directory where all output files will be written.                      |
| `vcf_name`  | `"joint"`             | The name of the output folder for results. Used to derive `JOINT_NAME`.    |

> 💡 **Note**: The `vcf_name` parameter is critical for naming conventions in downstream processes.

---

## Additional Samples Addition

| Parameter         | Default Value                                                                 | Description                                                                                     |
|-------------------|-------------------------------------------------------------------------------|-------------------------------------------------------------------------------------------------|
| `use_ref_gvcf`    | `true`                                                                      | Whether to use a reference GVCF file to include additional samples.                            |
| `ref_gvcf`        | `"${projectDir}/resources/ref_gvcfs/LineagesAndOutgroupV2.g.vcf.gz"`         | Path to the reference GVCF file.                                                              |
| `ref_gvcf_tbi`    | `"${projectDir}/resources/ref_gvcfs/LineagesAndOutgroupV2.g.vcf.gz.tbi"`     | Path to the index file for the reference GVCF.                                                |

> 💡 **Hint**: Use this feature if your dataset has low genetic diversity (e.g., clonal or fewer than 20 samples).

---

## Quality Control Parameters

| Parameter                     | Default Value | Description                                                                                     |
|-------------------------------|---------------|-------------------------------------------------------------------------------------------------|
| `cutoff_median_coverage`      | `10`          | The minimal median coverage required to process the sample.                                    |
| `cutoff_breadth_of_coverage`  | `0.90`        | The minimal breadth of coverage required to process the sample.                                |
| `cutoff_rel_abundance`        | `0.70`        | The minimal relative abundance of the majority strain required to process the sample.          |
| `cutoff_ntm_fraction`         | `0.20`        | The maximum fraction of NTM DNA allowed to process the sample.                                 |

> ⚠️ **Attention**: Ensure these values are adjusted based on the quality of your input data to avoid processing errors.

---

## Skipping Pipeline Steps

| Parameter                        | Default Value | Description                                                                                     |
|----------------------------------|---------------|-------------------------------------------------------------------------------------------------|
| `only_validate_fastqs`           | `false`       | Set to `true` to only validate input FASTQs and check their FASTQC reports.                    |
| `skip_merge_analysis`            | `false`       | Skip the final merge analysis step.                                                            |
| `skip_variant_recalibration`     | `false`       | Skip variant quality score recalibration (VQSR).                                               |
| `skip_base_recalibration`        | `true`        | Skip base quality score recalibration (BQSR). Not suitable for low-coverage Mtb genomes.       |
| `skip_minor_variants_gatk`       | `true`        | Skip minor variants detection with GATK. LoFreq is recommended for most purposes.              |
| `skip_phylogeny_and_clustering`  | `false`       | Disable downstream phylogenetic analysis of merged GVCF.                                       |
| `skip_complex_regions`           | `false`       | Disable downstream complex region analysis of merged GVCF.                                     |
| `skip_ntmprofiler`               | `false`       | Disable execution of `ntmprofiler` on FASTQ files.                                             |
| `skip_tbprofiler_fastq`          | `true`        | Disable `tbprofiler` analysis on FASTQ files.                                                  |
| `skip_spotyping`                 | `false`       | Disable spoligotyping analysis.                                                                |

> 💡 **Hint**: Use these flags to customize the pipeline execution based on your specific requirements.

---

## Reference Files

| Parameter               | Default Value                                              | Description                                                                                     |
|-------------------------|------------------------------------------------------------|-------------------------------------------------------------------------------------------------|
| `ref_fasta_basename`    | `"NC-000962-3-H37Rv"`                                      | Basename of the reference FASTA file.                                                          |
| `ref_fasta_dir`         | `"${projectDir}/resources/genome"`                         | Directory containing the reference FASTA file.                                                 |
| `ref_fasta`             | `"${params.ref_fasta_dir}/${params.ref_fasta_basename}.fa"`| Full path to the reference FASTA file.                                                         |
| `ref_fasta_dict`        | `"${params.ref_fasta_dir}/${params.ref_fasta_basename}.dict"`| Path to the reference FASTA dictionary file.                                                  |
| `ref_fasta_gb`          | `"${params.ref_fasta_dir}/${params.ref_fasta_basename}.gb"`| Path to the reference GenBank file.                                                            |

> ⚠️ **Warning**: It is recommended to use the provided reference files to ensure compatibility with the pipeline.

---

