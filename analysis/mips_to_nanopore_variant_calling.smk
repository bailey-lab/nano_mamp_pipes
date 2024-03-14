configfile: 'mips_to_nanopore.yaml'

"""
240207: Instead of splitting this smk and yaml into sections,  we may want to look into using mutiple smk files and calling them in 
this centralized smk using `include` or `subflow` to make the rule execution less complicated and code more 
modular(readable/debuggable/etc) (GT)
_______________________________________________________________________________________________________________________________________
"""

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
		finished_gvcf=config['output_dir']+'/variant_calling/nanomip_variants.gvcf.gz',
		finished_vcf=config['output_dir']+'/variant_calling/nanomip_variants.vcf.gz'
"""

"""
rule run_clair3_raw:
	"""
	call variants in all regions in sorted bams with clair3 and output gvcf and vcf files. Using 0.0001 as the minimum allele frequency for both snps and indels.
	"""
	input:
		samples_to_run=config["output_dir"]+'/mapping/bam_files/{sample}.sorted.bam',
		REF = config['genome_directory']+'/genome.fa', # TODO:make demo resources with Pf3D7 genome for built in tutorial demo data
		sif_path = config["clair3_sif_path"] # TODO: this can be handled by the singularity container eventually so users dont install clair3 themselves.
	output:
		vcf=config["output_dir"]+'/variant_calling/{sample}_raw'+'/merge_output.vcf.gz',
		gvcf=config["output_dir"]+'/variant_calling/{sample}_raw'+'/merge_output.gvcf.gz'
	params:
		INPUT_DIR = config["output_dir"]+'/mapping/bam_files',
		OUTPUT_DIR = config["output_dir"]+'/variant_calling/{sample}_raw',
		REF_DIR = config['genome_directory'],
		RESOURCE_DIR = config["project_resources"],
		THREADS = 4,
		MODEL_NAME = 'r941_prom_sup_g5014',
		targets_vcf = config["project_resources"]+'/targets.vcf.gz'
		
	shell:
		'''
		echo 'Processing sample: {wildcards.sample}'

		singularity exec \
			-B {params.INPUT_DIR} \
			-B {params.OUTPUT_DIR} \
			-B {params.REF_DIR} \
			-B {params.RESOURCE_DIR} \
			{input.sif_path} \
			/opt/bin/run_clair3.sh \
			--bam_fn={input.samples_to_run} \
			--ref_fn={input.REF} \
			--model_path=/opt/models/{params.MODEL_NAME} \
			--output={params.OUTPUT_DIR} \
			--threads={params.THREADS} \
			--platform=ont \
			--include_all_ctgs \
			--no_phasing_for_fa \
			--sample_name={wildcards.sample} \
			--gvcf \
			--snp_min_af=0.001 \
			--indel_min_af=0.001 \
			--include_all_ctgs \
			--threads=1 \
			--print_ref_calls
		'''

rule run_clair3_targets:
	"""
	call variants in targetted sites in same sorted bams with clair3 and output gvcf and vcf files. Must use 0 as the minimum allele frequency for both snps and indels to report coverage of reference calls in targetted sites (Ex. to report that there are no variants in pfk13 622I and we have X coverage to prove it.)
	"""
	input:
		samples_to_run=config["output_dir"]+'/mapping/bam_files/{sample}.sorted.bam',
		REF = config['genome_directory']+'/genome.fa', # TODO:make demo resources with Pf3D7 genome for built in tutorial demo data
		sif_path = config["clair3_sif_path"] # TODO: this can be handled by the singularity container eventually so users dont install clair3 themselves.
	output:
		vcf=config["output_dir"]+'/variant_calling/{sample}_targets'+'/merge_output.vcf.gz',
		gvcf=config["output_dir"]+'/variant_calling/{sample}_targets'+'/merge_output.gvcf.gz'
	params:
		INPUT_DIR = config["output_dir"]+'/mapping/bam_files',
		OUTPUT_DIR = config["output_dir"]+'/variant_calling/{sample}_targets',
		REF_DIR = config['genome_directory'],
		RESOURCE_DIR = config["project_resources"],
		THREADS = 4,
		MODEL_NAME = 'r941_prom_sup_g5014',
		targets_vcf = config["project_resources"]+'/targets.vcf.gz'
	shell:
		'''
		echo 'Processing sample: {wildcards.sample}'

		singularity exec \
			-B {params.INPUT_DIR} \
			-B {params.OUTPUT_DIR} \
			-B {params.REF_DIR} \
			-B {params.RESOURCE_DIR} \
			{input.sif_path} \
			/opt/bin/run_clair3.sh \
			--bam_fn={input.samples_to_run} \
			--ref_fn={input.REF} \
			--model_path=/opt/models/{params.MODEL_NAME} \
			--output={params.OUTPUT_DIR} \
			--threads={params.THREADS} \
			--platform=ont \
			--include_all_ctgs \
			--no_phasing_for_fa \
			--sample_name={wildcards.sample} \
			--gvcf \
			--include_all_ctgs \
			--threads=1 \
			--snp_min_af=0.0 \
			--indel_min_af=0.0 \
			--vcf_fn={params.targets_vcf} \
			--print_ref_calls
		'''

rule bcftools_index_gvcf_raw:
	input:
		original_gvcfs_to_index_raw=config["output_dir"]+'/variant_calling/{sample}_raw'+'/merge_output.gvcf.gz'
	output:
		index_file_gvcf_raw=config["output_dir"]+'/variant_calling/{sample}_raw'+'/merge_output.gvcf.gz.csi'
	shell:
		'''
		bcftools index {input.original_gvcfs_to_index_raw}
		'''

rule bcftools_index_gvcf_targets:
	input:
		original_gvcfs_to_index_targets=config["output_dir"]+'/variant_calling/{sample}_targets'+'/merge_output.gvcf.gz'
	output:
		index_file_gvcf_targets=config["output_dir"]+'/variant_calling/{sample}_targets'+'/merge_output.gvcf.gz.csi'
	shell:
		'''
		bcftools index {input.original_gvcfs_to_index_targets}
		'''

rule bcftools_index_vcf_raw:
	input:
		original_vcfs_to_index_raw=config["output_dir"]+'/variant_calling/{sample}_raw'+'/merge_output.vcf.gz'
	output:
		index_file_vcf_raw=config["output_dir"]+'/variant_calling/{sample}_raw'+'/merge_output.vcf.gz.csi'
	shell:
		'''
		bcftools index {input.original_vcfs_to_index_raw}
		'''

rule bcftools_index_vcf_targets:
	input:
		original_vcfs_to_index_targets=config["output_dir"]+'/variant_calling/{sample}_targets'+'/merge_output.vcf.gz'
	output:
		index_file_vcf_targets=config["output_dir"]+'/variant_calling/{sample}_targets'+'/merge_output.vcf.gz.csi'
	shell:
		'''
		bcftools index {input.original_vcfs_to_index_targets}
		'''

rule prepare_gvcf_lists:
	input: 
		#finished_index_file_gvcf_targets=config["output_dir"]+'/variant_calling/{sample}_targets'+'/merge_output.gvcf.gz.csi',
		finished_index_file_gvcf_targets=expand(config["output_dir"]+'/variant_calling/{sample}_targets'+'/merge_output.gvcf.gz.csi',sample=get_samples()),
		finished_index_file_gvcf_raw=expand(config["output_dir"]+'/variant_calling/{sample}_raw'+'/merge_output.gvcf.gz.csi',sample=get_samples()),
	output:
		raw_list = config["output_dir"] + '/gvcf_paths_raw.txt',
		targets_list = config["output_dir"] + '/gvcf_paths_targets.txt'
	params:
		py_script_path_gvcf = config["src_path"] + '/prepare_gvcf_lists.py',
		output_dir = config["output_dir"]
	shell:
		'''
		python {params.py_script_path_gvcf} {params.output_dir}
		'''

rule prepare_vcf_lists:
	input: 
		finished_index_file_vcf_targets=expand(config["output_dir"]+'/variant_calling/{sample}_targets'+'/merge_output.vcf.gz.csi',sample=get_samples()),
		finished_index_file_vcf_raw=expand(config["output_dir"]+'/variant_calling/{sample}_raw'+'/merge_output.vcf.gz.csi',sample=get_samples()),
	output:
		raw_list = config["output_dir"] + '/vcf_paths_raw.txt',
		targets_list = config["output_dir"] + '/vcf_paths_targets.txt'
	params:
		py_script_path_vcf = config["src_path"] + '/prepare_vcf_lists.py',
		output_dir = config["output_dir"]
	shell:
		'''
		python {params.py_script_path_vcf} {params.output_dir}
		'''

rule bcftools_merge_raw_gvcf:
	input:
		original_gvcfs_raw=expand(config["output_dir"]+'/variant_calling/{sample}_raw'+'/merge_output.gvcf.gz.csi',sample=get_samples()),
		raw_gvcf_list = config["output_dir"] + '/gvcf_paths_raw.txt'
	output:
		merged_gvcf_file_raw=config["output_dir"]+'/variant_calling/bcftools_merge_output_raw/merged_multisample.gvcf.gz'
	shell:
		'''
		bcftools merge --file-list {input.raw_gvcf_list} --force-samples -O z -o {output.merged_gvcf_file_raw}  
		'''

rule bcftools_merge_targets_gvcf:
	input:
		original_gvcfs_targets=expand(config["output_dir"]+'/variant_calling/{sample}_targets'+'/merge_output.gvcf.gz.csi',sample=get_samples()),
		targets_gvcf_list = config["output_dir"] + '/gvcf_paths_targets.txt'
	output:
		merged_gvcf_file_targets=config["output_dir"]+'/variant_calling/bcftools_merge_output_targets/merged_multisample.gvcf.gz'
	shell:
		'''
		bcftools merge --file-list {input.targets_gvcf_list} --force-samples -O z -o {output.merged_gvcf_file_targets} 
		'''

rule bcftools_merge_raw_vcf:
	input:
		original_vcfs=expand(config["output_dir"]+'/variant_calling/{sample}_raw'+'/merge_output.vcf.gz.csi',sample=get_samples()),
		raw_vcf_list = config["output_dir"] + '/vcf_paths_raw.txt'
	output:
		merged_vcf_file=config["output_dir"]+'/variant_calling/bcftools_merge_output_raw/merged_multisample.vcf.gz'
	shell:
		'''
		bcftools merge --file-list {input.raw_vcf_list} --force-samples -O z -o {output.merged_vcf_file} 
		'''

rule bcftools_merge_targets_vcf:
	input:
		original_vcfs=expand(config["output_dir"]+'/variant_calling/{sample}_targets'+'/merge_output.vcf.gz.csi',sample=get_samples()),
		targets_vcf_list = config["output_dir"] + '/vcf_paths_targets.txt'
	output:
		merged_vcf_file=config["output_dir"]+'/variant_calling/bcftools_merge_output_targets/merged_multisample.vcf.gz'
	shell:
		'''
		bcftools merge --file-list {input.targets_vcf_list} --force-samples -O z -o {output.merged_vcf_file}
		'''

rule bcftools_index_merged_raw_gvcf:
	input:
		merged_gvcfs_to_index_raw=config["output_dir"]+'/variant_calling/bcftools_merge_output_raw/merged_multisample.gvcf.gz'
	output:
		index_file_gvcf_targets=config["output_dir"]+'/variant_calling/bcftools_merge_output_raw/merged_multisample.gvcf.gz.csi'
	shell:
		'''
		bcftools index {input.merged_gvcfs_to_index_raw}
		'''

rule bcftools_index_merged_targets_gvcf:
	input:
		merged_gvcfs_to_index_targets=config["output_dir"]+'/variant_calling/bcftools_merge_output_targets/merged_multisample.gvcf.gz'
	output:
		index_file_gvcf_targets=config["output_dir"]+'/variant_calling/bcftools_merge_output_targets/merged_multisample.gvcf.gz.csi'
	shell:
		'''
		bcftools index {input.merged_gvcfs_to_index_targets}
		'''

rule bcftools_index_merged_raw_vcf:
	input:
		merged_vcfs_to_index_raw=config["output_dir"]+'/variant_calling/bcftools_merge_output_raw/merged_multisample.vcf.gz'
	output:
		index_file_vcf_targets=config["output_dir"]+'/variant_calling/bcftools_merge_output_raw/merged_multisample.vcf.gz.csi'
	shell:
		'''
		bcftools index {input.merged_vcfs_to_index_raw}
		'''

rule bcftools_index_merged_targets_vcf:
	input:
		merged_vcfs_to_index_targets=config["output_dir"]+'/variant_calling/bcftools_merge_output_targets/merged_multisample.vcf.gz'
	output:
		index_file_vcf_targets=config["output_dir"]+'/variant_calling/bcftools_merge_output_targets/merged_multisample.vcf.gz.csi'
	shell:
		'''
		bcftools index {input.merged_vcfs_to_index_targets}
		'''

rule bcftools_merge_raw_and_targets_gvcf:
	input:
		merged_raw_gvcf=config["output_dir"]+'/variant_calling/bcftools_merge_output_raw/merged_multisample.gvcf.gz',
		merged_targets_gvcf=config["output_dir"]+'/variant_calling/bcftools_merge_output_targets/merged_multisample.gvcf.gz',
		index_raw_gvcf=config["output_dir"]+'/variant_calling/bcftools_merge_output_raw/merged_multisample.gvcf.gz.csi',
		index_targets_gvcf=config["output_dir"]+'/variant_calling/bcftools_merge_output_targets/merged_multisample.gvcf.gz.csi'
	output:
		merged_raw_and_targets_gvcf_file=config["output_dir"]+'/variant_calling/nanomip_variants.gvcf.gz'
	shell:
		'''
		bcftools merge --force-samples -O z -o {output.merged_raw_and_targets_gvcf_file} {input.merged_raw_gvcf} {input.merged_targets_gvcf}
		'''

rule bcftools_merge_raw_and_targets_vcf:
	input:
		merged_raw_vcf=config["output_dir"]+'/variant_calling/bcftools_merge_output_raw/merged_multisample.vcf.gz',
		merged_targets_vcf=config["output_dir"]+'/variant_calling/bcftools_merge_output_targets/merged_multisample.vcf.gz',
		index_raw_vcf=config["output_dir"]+'/variant_calling/bcftools_merge_output_raw/merged_multisample.vcf.gz.csi',
		index_targets_vcf=config["output_dir"]+'/variant_calling/bcftools_merge_output_targets/merged_multisample.vcf.gz.csi'
	output:
		merged_raw_and_targets_vcf_file=config["output_dir"]+'/variant_calling/nanomip_variants.vcf.gz'
	shell:
		'''
		bcftools merge --force-samples -O z -o {output.merged_raw_and_targets_vcf_file} {input.merged_raw_vcf} {input.merged_targets_vcf}
		'''

