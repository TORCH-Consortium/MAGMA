# Input and Output
The input CSV sample file (the study id cannot start with 'XBS_REF_')
> :bulb: **Hint**: <br>
> The samplesheet should have the following fields [Sample, R1, R2], optionally you could add other fields [study, library, attempt, flowcell, lane, index_sequence]

> :bulb: **NOTE**: <br>
> Most of these parameters are used to create unique_id in XBS_main.py

`input_samplesheet`: "samplesheet.magma.csv"

The directory to which all output files should be written

`outdir`: "magma-results"

The name of the output folder where all results will be allocated

`vcf_name`: "joint"
> :bulb: **NOTE**: <br> This parameter is used to derive the JOINT_NAME in XBS_main.py

# Additional samples addition

> :bulb: **Hint**: <br> Got little genetic diversity in your dataset? (.e.g clonal or <20 samples) - use the EXIT-RIF GVCF file to include additional samples.

`use_ref_gvcf`: true  
`ref_gvcf`:  "${projectDir}/resources/ref_gvcfs/LineagesAndOutgroupV2.g.vcf.gz"  
`ref_gvcf_tbi`:  "${projectDir}/resources/ref_gvcfs/LineagesAndOutgroupV2.g.vcf.gz.tbi"

# Optional parameters

## Adjusting Quality control parameters
`cutoff_median_coverage`: 10  

The minimal median coverage required to process the sample.

`cutoff_breadth_of_coverage`: 0.90  

The minimal breadth of coverage required to process the sample.

`cutoff_rel_abundance`: 0.70  
The minimal relative abundance of the majority strain required to process the sample.

`cutoff_ntm_fraction`: 0.20  

The maximum fraction of NTM DNA allowed to process the sample.

`cutoff_site_representation`: 0.95  

The minimum fraction of samples that need to have a call at a site before the site is considered in phylogeny.

`cutoff_strand_bias`: 0  

The p-value below which binomial strand bias is considered significant (switched to zero due to FP mixed infection issues).

## Skipping pipeline steps

`only_validate_fastqs`: false OR true  

Set this to true if you'd like to only validate input FASTQs and check their FASTQC reports.

`skip_merge_analysis`: false OR true  

Use this flag to skip the final merge analysis.

`skip_variant_recalibration`: false OR true  

MAGMA optimizes VQSR. If this messes up, use the default settings for VQSR.

`skip_base_recalibration`: true  
Do NOT enable for BQSR Mtb. BQSR is not suitable for low-coverage Mtb genomes or when contamination is present.

`skip_minor_variants_gatk`: true  
Output file for minor variants detection from BAM with GATK. LoFreq does a better job for most purposes.

`skip_phylogeny_and_clustering`: false OR true  
Use this flag to disable downstream phylogenetic analysis of merged GVCF.

`skip_complex_regions`: false OR true  
Use this flag to enable downstream complex region analysis of merged GVCF.

`skip_ntmprofiler`: false OR true  
Enable execution of ntmprofiller on FASTQ files.

`skip_tbprofiler_fastq`: true OR false  
Enable or disable tbprofiler analysis on FASTQ files (with WHO+ database) on FASTQ files.

`skip_spotyping`: false OR true  
Enable or disable spoligotyping analysis.

## Flags for experimental features

### RDAnalyzer
`skip_rdanalyzer`: true  
Enable or disable RDAnalyzer analysis for Region of Difference identification (*NOT working yet*).  

`ref_fasta_rdanalyzer`: "${projectDir}/resources/rdanalyzer/RDs30.fasta"  
Reference FASTA file for RDAnalyzer.

## IQTREE Parameters

> :bulb: **NOTE**: PICK ONE of the following parameters related to IQTREE.

`iqtree_standard_bootstrap`: true  
Enable standard bootstrap analysis.  

`iqtree_fast_ml_only`: false  
Enable fast maximum likelihood analysis only.  

`iqtree_fast_bootstrapped_phylogeny`: false  
Enable fast bootstrapped phylogeny analysis.  

`iqtree_accurate_ml_only`: false  
Enable accurate maximum likelihood analysis only.  

`iqtree_custom_argument`:  
Custom arguments for IQTREE can be specified here.  

# Specific paths to reference files

> :bulb: **NOTE**: You can customize reference files used to the reads mapping
> :warn: **WARN**: It is best not to change this parameters and to rely upon the provided reference files that are packed within the pipeline


`ref_fasta_basename`: "NC-000962-3-H37Rv"
`ref_fasta_dir`: "${projectDir}/resources/genome"
`ref_fasta_gb`: "${params.ref_fasta_dir}/${params.ref_fasta_basename}.gb"
`ref_fasta_dict`: "${params.ref_fasta_dir}/${params.ref_fasta_basename}.dict"
`ref_fasta`: "${params.ref_fasta_dir}/${params.ref_fasta_basename}.fa"
`ref_fasta_amb`: "${params.ref_fasta}.amb"
`ref_fasta_ann`: "${params.ref_fasta}.ann"
`ref_fasta_bwt`: "${params.ref_fasta}.bwt"
`ref_fasta_fai`: "${params.ref_fasta}.fai"
`ref_fasta_pac`: "${params.ref_fasta}.pac"
`ref_fasta_sa`: "${params.ref_fasta}.sa"

`queries_multifasta` = "${projectDir}/resources/regions/IS_Queries_Mtb.multi.fasta"


`drgenes_list`: "${projectDir}/resources/regions/WHO_Tier1_Tier2_DR.list"

`rrna_list`: "${projectDir}/resources/regions/rRNA.list"

`dbsnp_vcf`: "${projectDir}/resources/known/Benavente2015.UVPapproved.rRNAexcluded.vcf.gz"

`dbsnp_vcf_tbi`:  "${params.dbsnp_vcf}.tbi"

`excluded_loci_list`: "${projectDir}/resources/regions/UVP_List_of_Excluded_loci.list"



# Process level configuration

> :bulb: **NOTE**: You can customize parameters specific to process, modifying the standard behavior of a tool
> :warn: **WARN**: Advanced parameters


### BWA_MEM
`BWA_MEM.arguments`: "-k 100"  

### BWA_MEM__DELLY
  
`BWA_MEM__DELLY.arguments`: ""  

### GATK_HAPLOTYPE_CALLER

`GATK_HAPLOTYPE_CALLER.arguments`: "-ploidy 1 --dont-use-soft-clipped-bases --read-filter MappingQualityNotZeroReadFilter -G StandardAnnotation -G AS_StandardAnnotation"  

### GATK_HAPLOTYPE_CALLER__MINOR_VARIANTS

`GATK_HAPLOTYPE_CALLER__MINOR_VARIANTS.arguments`: "-ploidy 1 \
  --minimum-mapping-quality 60 \
  --min-base-quality-score 20 \
  --read-filter MappingQualityNotZeroReadFilter \
  -G StandardAnnotation \
  --dont-use-soft-clipped-bases \
  --output-mode EMIT_ALL_ACTIVE_SITES"  

### LOFREQ_CALL__NTM
`LOFREQ_CALL__NTM.region`: "1472307-1472307"  
`LOFREQ_CALL__NTM.arguments`: "-m 60 -Q 20 -a 1"  

### LOFREQ_INDELQUAL
`LOFREQ_INDELQUAL.arguments`: "-m 60"  

### LOFREQ_CALL 
`LOFREQ_CALL.arguments`: "-m 60 --call-indels -Q 30"  

### LOFREQ_FILTER
`LOFREQ_FILTER.arguments`: "-a 0.20"  

### DELLY_CALL
`DELLY_CALL.arguments`: "-u 30"  

### BCFTOOLS_VIEW__DELLY 
`BCFTOOLS_VIEW__DELLY.arguments`: """-i 'GT="1/1"'"""  

### SAMTOOLS_STATS
`SAMTOOLS_STATS.arguments`: "-F DUP,SUPPLEMENTARY,SECONDARY,UNMAP,QCFAIL"  

### GATK_COLLECT_WGS_METRICS
`GATK_COLLECT_WGS_METRICS.arguments`: "--READ_LENGTH 0 --COVERAGE_CAP 10000 --COUNT_UNPAIRED"  

### GATK_SELECT_VARIANTS__SNP
`GATK_SELECT_VARIANTS__SNP.arguments`: "--remove-unused-alternates --exclude-non-variants"  

### GATK_SELECT_VARIANTS__INDEL
`GATK_SELECT_VARIANTS__INDEL.arguments`: "--remove-unused-alternates --exclude-non-variants --select-type-to-include MNP --select-type-to-include MIXED"  

### GATK_VARIANT_RECALIBRATOR__SNP
`GATK_VARIANT_RECALIBRATOR__SNP.arguments`: "--use-allele-specific-annotations \
  -AS \
  --target-titv 1.7 \
  --truth-sensitivity-tranche 100.0 \
  --truth-sensitivity-tranche 99.9 \
  --truth-sensitivity-tranche 99.8 \
  --truth-sensitivity-tranche 99.7 \
  --truth-sensitivity-tranche 99.6 \
  --truth-sensitivity-tranche 99.5 \
  --truth-sensitivity-tranche 99.4 \
  --truth-sensitivity-tranche 99.3 \
  --truth-sensitivity-tranche 99.2 \
  --truth-sensitivity-tranche 99.1 \
  --truth-sensitivity-tranche 99.0 \
  --max-gaussians 4 \
  -mq-cap 60"  

### GATK_VARIANT_RECALIBRATOR__INDEL
`GATK_VARIANT_RECALIBRATOR__INDEL.arguments`: "-AS \
  --target-titv 1.8 \
  --truth-sensitivity-tranche 100.0 \
  --truth-sensitivity-tranche 99.9 \
  --truth-sensitivity-tranche 99.8 \
  --truth-sensitivity-tranche 99.7 \
  --truth-sensitivity-tranche 99.6 \
  --truth-sensitivity-tranche 99.5 \
  --truth-sensitivity-tranche 99.4 \
  --truth-sensitivity-tranche 99.3 \
  --truth-sensitivity-tranche 99.2 \
  --truth-sensitivity-tranche 99.1 \
  --truth-sensitivity-tranche 99.0 \
  --max-gaussians 4 \
  -mq-cap 60"  

### GATK_APPLY_VQSR__SNP
`GATK_APPLY_VQSR__SNP.arguments`: "--ts-filter-level 99.90 -AS --exclude-filtered"  

### GATK_APPLY_VQSR__INDEL
`GATK_APPLY_VQSR__INDEL.arguments`: ""  

### GATK_SELECT_VARIANTS__EXCLUSION__SNP
`GATK_SELECT_VARIANTS__EXCLUSION__SNP.arguments`: ""  

### GATK_SELECT_VARIANTS__EXCLUSION__INDEL
`GATK_SELECT_VARIANTS__EXCLUSION__INDEL.arguments`: "--select-type-to-include MNP --select-type-to-include MIXED"  

### BCFTOOLS_MERGE__LOFREQ
`BCFTOOLS_MERGE__LOFREQ.file_format`: "lofreq"  

### BCFTOOLS_MERGE__DELLY
`BCFTOOLS_MERGE__DELLY.file_format`: "delly"  

### TBPROFILER_VCF_PROFILE__COHORT 
`TBPROFILER_VCF_PROFILE__COHORT.arguments`: "--depth 0,0 --af 0,0 --strand 0 --sv_depth 0,0 --sv_af 0,0 --sv_len 100000,50000"  

### TBPROFILER_COLLATE__COHORT
`TBPROFILER_COLLATE__COHORT.prefix`: "major_variants"  

### TBPROFILER_FASTQ_PROFILE
`TBPROFILER_FASTQ_PROFILE.arguments`: "--csv"  

### TBPROFILER_FASTQ_COLLATE
`TBPROFILER_FASTQ_COLLATE.prefix`: "fastq"  

### SPOTYPING
`SPOTYPING.arguments`: "" OR "--noQuery"  

### UTILS_CAT_SPOTYPING 
`UTILS_CAT_SPOTYPING.arguments`: ""  

### RDANALYZER
`RDANALYZER.arguments`: ""  

### TBPROFILER_VCF_PROFILE__LOFREQ 
`TBPROFILER_VCF_PROFILE__LOFREQ.arguments`: "--depth 0,0 --af 0,0 --strand 0 --sv_depth 0,0 --sv_af 0,0 --sv_len 100000,50000"  

### TBPROFILER_COLLATE__LOFREQ
`TBPROFILER_COLLATE__LOFREQ.prefix`: "minor_variants"  

### TBPROFILER_VCF_PROFILE__DELLY 
`TBPROFILER_VCF_PROFILE__DELLY.arguments`: "--depth 0,0 --af 0,0 --strand 0 --sv_depth 0,0 --sv_af 0,0 --sv_len 100000,50000"  

### TBPROFILER_COLLATE__DELLY 
`TBPROFILER_COLLATE__DELLY.prefix`: "structural_variants"  

### GATK_SELECT_VARIANTS__PHYLOGENY
`GATK_SELECT_VARIANTS__PHYLOGENY.arguments`: "--remove-unused-alternates --exclude-non-variants"   

### GATK_VARIANTS_TO_TABLE
`GATK_VARIANTS_TO_TABLE.arguments`: "-GF GT"   

### CLUSTERPICKER
`CLUSTERPICKER.bootstrap_1`: 0  
`CLUSTERPICKER.bootstrap_2`: 0  
`CLUSTERPICKER.max_cluster_size`: 0  
`CLUSTERPICKER.algorithm`: "gap"
