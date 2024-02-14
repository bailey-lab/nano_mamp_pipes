#!/bin/bash
#SBATCH --mem=48G
#SBATCH -p gpu --gres=gpu:1 --gres-flags=enforce-binding
#SBATCH -N 1
#SBATCH -n 8
#SBATCH -t 48:00:00
input_dir=$1
output_dir=$2

export PATH=/nfs/jbailey5/baileyweb/bailey_share/bin/ont-guppy_6.4.8/bin:$PATH
module load cuda/12.2.0 gcc/10.1.0
guppy_basecaller -x "cuda:0" -c dna_r10.4.1_e8.2_400bps_sup.cfg --input_path $input_dir --save_path $output_dir --compress_fastq

