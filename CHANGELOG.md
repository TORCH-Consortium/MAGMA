# CHANGELOG FOR THE MAGMA PIPELINE VERSIONS
<!-- https://keepachangelog.com/en/1.1.0/ -->


## unreleased

1. **Optional BAM → CRAM checkpoint after MarkDup / BQSR** (Pattern 1 of the BAM-checkpoints spec; see `abc-universe/specs/active/magma-bam-checkpoints.md`). Opt-in via `--bam_checkpoint_compression cram` (default `'none'` keeps current behaviour exactly). When enabled, a new local module `SAMTOOLS_VIEW_CRAM` replaces `SAMTOOLS_INDEX` at the post-MarkDup-or-BQSR junction in `workflows/call_wf.nf` and emits the identical channel shape `(sampleName, *.cram.crai, *.cram)` instead of `(sampleName, *.bai, *.bam)`. GATK HaplotypeCaller, the minor-variants HC, and LoFreq (htslib-backed) read the `.cram` directly via `-R / -f <ref>` — zero round-trip cost at consume time. Lossless CRAM 3.1 by default; quality-score binning is opt-in via `--cram_lossy_qualities true` (use only for archival of already-called cohorts).

2. New params: `bam_checkpoint_compression` (`'none' | 'cram'`), `cram_lossy_qualities` (boolean). Both default to safe values; no behaviour change for existing invocations.


## v2.0.0

1. Update `tb-profiler` to `v6.2.1` along with the `WHO mutations catalog v2` with the database version `30f8bc37df15affa378ebbfbd3e1eb4c5903056e`
 
2. Addition of `ntm-profiler` with the database FIXME

3. Two new `conda-envs` and containers have been added to the pipeline, relying directly on upstream biocontainers for both `ntm-profiler` and `tb-profiler`

4. Changes to the directory structure

5. Addition of mixed-infections resistance summary script

6. Addition of `low_memory` profile to accommodate lower end infrastructure

7. Updates to the names of parameters for triggering partial/optional workflows


## v1.1.1

Created a parallel workflow for mapping without using the strict seed lenght for use in the structural variant workflow.

Updated TBProfiler to version 5.0.0 and recreated the resistance database to work with the the new version

Updated the summarize resistance script to include the structural variants in the excel output


## v1.0.0

Initial release of the pipeline
