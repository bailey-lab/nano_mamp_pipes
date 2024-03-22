#!/bin/bash
#SBATCH --job-name=dorado_basecalling
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=dorado_basecalling_%j.log

module load cuda/12.2.0 gcc/10.1.0

# convert fast5 to pod5 - do this once, alec can write these pod5 in future - fast5 is legacy in favor of pod5 according to ONT P2 docs
# pod5 convert fast5 is the command. Built into convert_fast5_to_pod5.sh
# pod5 is needed for full performance of dorado basecaller

POD5_DIR="/nfs/jbailey5/baileyweb/gtollefs/mips_to_nanopore_gt/240215_MtoN_P2/output/pod5/pod5_output.pod5"
OUTPUT_DIR="/nfs/jbailey5/baileyweb/gtollefs/mips_to_nanopore_gt/240215_MtoN_P2/output/dorado_output"
# want res_dna_r10.4.1_e8.2_400bps_sup@v4.0.1 - looks like  downloading dna_r10.4.1_e8.2_400bps_sup@v4.1.0 with httplib is downloaded automaticlly. 4.0.1 says not avaiable
# run dorado basecaller
/nfs/jbailey5/baileyweb/gtollefs/toolshed/dorado_basecaller/dorado-0.5.3-linux-x64/bin/dorado basecaller --device "cuda:all" sup $POD5_DIR --emit-fastq > ${OUTPUT_DIR}/basecalls.fastq
