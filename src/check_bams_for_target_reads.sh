#!/bin/bash
#SBATCH --job-name=check_bam_coverage
#SBATCH --output=check_bam_coverage_%j.out
#SBATCH --error=check_bam_coverage_%j.err
#SBATCH --time=12:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=4 

module load samtools

# testing: set paths
# out_dir=/nfs/jbailey5/baileyweb/gtollefs/mips_to_nanopore_gt/240215_MtoN_P2/second_output_rerun_for_completeness_check
# resource_dir=/nfs/jbailey5/baileyweb/gtollefs/mips_to_nanopore_gt/bams_to_vcf/resources

# get paths as arguments
out_dir=$1
resource_dir=$2

VCF_FILE=${resource_dir}/targets.vcf.gz 
BAM_DIR=${out_dir}/mapping/bam_files

# make directory for this check step output
mkdir -p ${out_dir}/samples_for_targeted_calling/

# get positions from compressed vcf
zcat ${VCF_FILE} | awk '{if (!/^#/) print $1":"$2"-"$2}' > ${out_dir}/samples_for_targeted_calling/target_positions.txt

# initialize txt log for bam files with coverage on target snps
>  ${out_dir}/samples_for_targeted_calling/bam_with_coverage.txt

# function to check coverage
check_coverage() {
    local bam_file=$1
    local positions_file=$2
    for position in $(cat ${positions_file}); do
        # Check if the depth at this position is greater than 0
        if samtools depth -r ${position} ${bam_file} | awk '$3 > 0 {print; exit}'; then
            echo found reads in ${bam_file}
            echo ${bam_file} >> ${out_dir}/samples_for_targeted_calling/bam_with_coverage.txt
            break
        fi
    done
}


# Export the function so it's available to parallel
export -f check_coverage

find ${BAM_DIR} -name "*.bam" | while read bam_file; do
  check_coverage "${bam_file}" ${out_dir}/samples_for_targeted_calling/target_positions.txt
done

# Check each BAM file in parallel
# find ${BAM_DIR} -name "*.bam" | parallel check_coverage {} positions.txt

echo "bam files with coverage on at least one snp in specified targets file have been logged in ${out_dir}/samples_for_targeted_calling/bam_with_coverage.txt"
