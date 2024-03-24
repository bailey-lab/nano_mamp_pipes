# NanoMAMP Pipes

## Overview

NanoMAMP Pipes is a fully automated bioinformatics pipeline built in snakemake designed for comprehensive analysis of short read targeted sequencing data from Nanopore platforms. This streamlined tool takes you from raw Nanopore output in FAST5 or POD5 formats through basecalling, custom sample barcode demultiplexing, alignment, variant calling, and post-variant analysis with integrated quality control checkpoints. NanoMAMP Pipes is designed to be run on a shared computing cluster running the SLURM Workload Manager.
## Features

**Single-Step Execution:** Run your entire analysis with one command.

**Basecalling with Dorado:** Performs simplex basecalling with Dorado using a model tailored to your ONT flowcell and chemistry versions.

**Custom Demultiplexing:** Supports custom sample barcodes for accurately demultiplexing pooled samples.

**Minimap2 Alignment:** Utilizes the speed and efficiency of latest Minimap2 using models tailed to  for read alignment in short read mode.

**Clair3 Variant Calling:** Employs Clair3 with pretrained models for high-precision variant calling.

**Variant Post Analysis:** Includes prevalence calculation for actionable insights.

**Quality Control:** Provides sequencing quality summaries for confidence in data integrity.
## Quick Start

### Prerequisites

* Snakemake (tested on v7.28.3) 
```bash
mamba create -c conda-forge -c bioconda -n nanomamp snakemake
```
* Minimap2 (tested on v2.26 (r1175)) 
```bash
mamba install -n nanomamp minimap2
```
* Clair3 (tested on v1.0.5) 

Install the Clair3 Singularity container following [these instructions](https://github.com/HKU-BAL/Clair3?tab=readme-ov-file#option-2-singularity) and store the sif file in a convenient place. You'll specify this sif path in the yaml file
### Install NanoMAMP Pipes

```bash
git clone https://github.com/bailey-lab/nano_mamp_pipes.git
```
### Run NanoMAMP Pipeline

Run the full automated pipeline from POD5 to VCF post-analysis in a single command following these steps:

1. Set input paths and parameters in the YAML configuration file stored in  /nano_mamp_pipes/analysis/ 

The `mips_to_nanopore.yaml` file is the central place for setting up your NanoMAMP Pipes analysis. Carefully check that you've input the correct paths and resources before running the pipeline.

In the YAML file, you'll specify:
* `fastq_dir: /path/to/demuxed_fastq/fastq`
* `output_dir: /path/to/pipeline_output`
* `project_resources: /path/to/project_resources`
* `genome_directory: /path/to/genome_resources`
* `clair3_sif_path: /path/to/clair3_latest.sif`
* `src_path: /path/to/pipeline_scripts`

2. Run the automated one-step bash script which launches sequential snakemake scripts with automated parallelization and sample tracking from the /nano_mamp_pipes/analysis/ directory
```bash
bash run_pipeline.sh
```

### Run Modular Steps of the Pipeline

1. Activate the nanomamp conda environment

```
mamba activate nanomamp
```

2. Run the snakemake script for your step of choice from the /nano_mamp_pipes/analysis/ directory 
```
snakemake -s demux_nanopore_dorado_version.smk --profile slurm
```
```
snakemake -s mips_to_nanopore_alignment.smk --profile slurm
```
```
snakemake -s mips_to_nanopore_variant_calling.smk --profile slurm
```

# NanoMAMP Pipeline Diagram
![pipeline_wireframe_diagram](https://github.com/bailey-lab/nano_mamp_pipes/blob/main/src/resources/mips_to_nanopore.drawio.svg?raw=true "NanoMAMP Pipeline Wireframe")





