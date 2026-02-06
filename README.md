# RNAQuant_NF

RNA-seq analysis pipeline based on Nextflow (DSL2). It supports QC, read filtering, alignment, strand inference, gene counting, merge/TPM/FPKM, and MultiQC reporting. It can also skip alignment and start from BAM inputs.

## Requirements
- Nextflow (DSL2)
- Conda/Singularity/Docker with tools: `fastqc`, `fastp`, `hisat2`, `STAR`, `samtools`, `infer_experiment.py`, `featureCounts`, `multiqc`

## Quick Start
```bash
nextflow run main.nf
```

## Parameters
- `--input` : samplesheet TSV with columns `ID`, `R1`, `R2`
- `--project` : project name prefix
- `--outdir` : output directory
- `--ref` : reference key in `nextflow.config` (`hsa`/`mmu`)
- `--aligner` : `star` (default) or `hisat2`
- `--skip_align` : skip alignment and start from BAMs
- `--bam_list` : TSV with columns `ID`, `BAM` (required if `--skip_align` is set)

## Samplesheet Format
```tsv
ID\tR1\tR2
sample1\t/path/sample1_R1.fastq.gz\t/path/sample1_R2.fastq.gz
```

## BAM List Format (skip alignment)
```tsv
ID\tBAM
sample1\t/path/sample1.bam
sample2\t/path/sample2.bam
```

## Run Examples
```bash
# Standard run
nextflow run main.nf --input /path/to/samplesheet.tsv --ref hsa --aligner star

# Start from BAM
nextflow run main.nf --skip_align --bam_list /path/to/bam_list.tsv
```

## Pipeline Structure
- `modules/` : process modules
- `workflows/` : sub-workflows
- `main.nf` : entry workflow
- `nextflow.config` : reference and process configuration

## Outputs
- `results/fastqc` : FastQC reports
- `results/clean` : fastp reports
- `results/align` : alignment outputs
- `results/count` : featureCounts outputs
- `results/merge` : merged count/TPM/FPKM
- `results/reports` : MultiQC reports

## Notes
- Reference paths are configured in `nextflow.config` under `params.references`.
- When `--skip_align` is used, BAM files will be indexed if `.bai` is missing.
