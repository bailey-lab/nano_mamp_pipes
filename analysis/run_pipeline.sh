#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=4G
#SBATCH -t 96:00:00
#SBATCH -N 1
#SBATCH -J nano_mamp_pipes
#SBATCH --mail-type=END

################################################################################
#
# UNDER DEVELOPMENT: THIS SCRIPT IS NOT FULLY FUNCTIONING AND TESTED YET
#
################################################################################

# Batch script to run all steps of start to finish mips to nanopore pipeline including:
# This does not need lots of resources, but does need to run for as long as possible to allow submitting jobs for each smk step (96 hour for priority users)
# 1. Run dorado basecalling
# 2. Run demultiplexing
# 3. Run minimap2 alignment
# 4. Run clair3 variant calling
# 5. Run vcf_to_tables and miplithon combined .py script

# define directory of batch script so dorado_bashcalling.sh can be found and called
# Inside the script
if [ -n "${SLURM_JOB_ID:-}" ]; then
  COMMAND=$(scontrol show job "$SLURM_JOB_ID" | grep Command= | sed 's/Command=//')
  # Assuming the script's path is two levels deep from the desired "root"
  THEPATH=$(dirname "$(dirname "$(dirname "$COMMAND")")")
else
  # Similar assumption for direct execution
  THEPATH=$(dirname "$(dirname "$(dirname "$(realpath "$0")")")")
fi

# this find dirname one back from current dir holding this batch script. 
# assuming here that this run_pipeline.sh is in /nano_mamp_pipes/analysis and smk files are in /nano_mamp_pipes/src
echo "$THEPATH"/analysis

mamba activate snakemake

# run dorado basecalling
# Note: need to run dorado with gpu but not the rest, so dorado/guppy call is to be made as a batch script stored in src that requests GPU resources
sbatch "$THEPATH"/src/dorado_basecalling_dev.sh 

# run demultiplexing
snakemake -s "$THEPATH"/src/demux_nanopore_dorado_version.smk --profile slurm

# run QC script - find_failed_reads.py
sbatch "$THEPATH"/src/find_failed_reads.sh

# run minimap2 alignment
snakemake -s  "$THEPATH"/src/mips_to_nanopore_alignment.smk --profile slurm

# run check_bams_for_target_reads.sh (or could add this to the end of the alignment smk)
sbatch  "$THEPATH"/src/check_bams_for_target_reads.sh

# run clair3 variant calling
snakemake -s  "$THEPATH"/src/mips_to_nanopore_variant_calling.smk --profile slurm

# run vcf_to_tables and miplithon combined .py script
snakemake -s  "$THEPATH"/src/variant_post_analysis.smk --profile slurm # this smk still to be written
