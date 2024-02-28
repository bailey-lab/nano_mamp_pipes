config=$1
input_dir=$2
output_dir=$3

export PATH=/nfs/jbailey5/baileyweb/bailey_share/bin/ont-guppy_6.4.8/bin:$PATH
module load cuda/12.2.0 gcc/10.1.0
guppy_basecaller -x "cuda:all" -c $config --input_path $input_dir --save_path $output_dir --compress_fastq

