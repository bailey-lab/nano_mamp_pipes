configfile: 'mips_to_nanopore.yaml'

import os
import glob

# get sample prefixes
def get_samples():
    fastq_files = glob.glob(os.path.join(config['fastq_dir'], '*.fastq.gz'))
    sample_names = [
        os.path.basename(f).replace('.fastq.gz', '') for f in fastq_files
        if 'Undetermined' not in f and 'NTC' not in os.path.basename(f)
    ]
    return sample_names

# Rule all to specify final output files
rule all:
    input:
        sorted_bams=expand(config["output_dir"]+'/mapping/bam_files/{sample}.sorted.bam', sample=get_samples()),
        bam_index=expand(config["output_dir"]+'/mapping/bam_files/{sample}.sorted.bam.bai', sample=get_samples()),
        bams_with_target_coverage_list=config['output_dir'] + "/samples_for_targeted_calling/bam_with_coverage.txt"

rule make_sam:
	"""
	Rule for aligning reads with minimap2 to sam files
	"""
	input:
		indexed_genome=config['genome_directory']+'/genome.fa',
		sample_fastq=config['fastq_dir']+'/{sample}.fastq.gz'
	output:
		sample_sam=temp(config["output_dir"]+'/mapping/sam_files/{sample}.sam')
	shell:
		'minimap2 -ax sr {input.indexed_genome} {input.sample_fastq} -o {output.sample_sam}'

rule make_bam:
	"""
	Convert sam to bam with samtools view -b
	"""
	input:
		sample_sam=config["output_dir"]+'/mapping/sam_files/{sample}.sam'
	output:
		sample_bam=temp(config["output_dir"]+'/mapping/bam_files/{sample}.bam')
	shell:
		'samtools view -b -o {output.sample_bam} {input.sample_sam}'

rule sort_bam:
	"""
	Sort bams by coordinate with samtools sort
	"""
	input:
		sample_bam_to_sort=config["output_dir"]+'/mapping/bam_files/{sample}.bam'
	output:
		sorted_bam=config["output_dir"]+'/mapping/bam_files/{sample}.sorted.bam'
	shell:
		'samtools sort -o {output.sorted_bam} {input.sample_bam_to_sort}'

rule index_bam:
	"""
	index bams with samtools index
	"""
	input:
		sorted_bam=config["output_dir"]+'/mapping/bam_files/{sample}.sorted.bam'
	output:
		bam_index=config["output_dir"]+'/mapping/bam_files/{sample}.sorted.bam.bai'
	shell:
		'samtools index {input.sorted_bam}'

rule check_bam_coverage:
    input:
        bam_files=expand(config["output_dir"] + '/mapping/bam_files/{sample}.sorted.bam', sample=get_samples())
    output:
        bams_with_target_coverage_list=config['output_dir'] + "/samples_for_targeted_calling/bam_with_coverage.txt"
    params:
        out_dir=config["output_dir"],
        resource_dir=config["project_resources"],
        src_directory=config["src_path"]
    shell:
        """
        sbatch --export=ALL {params.src_directory}/check_bams_for_target_reads.sh {params.out_dir} {params.resource_dir}
        """
