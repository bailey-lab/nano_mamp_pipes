configfile: 'mips_to_nanopore.yaml'

"""
Development Notebook:
_______________________________________________________________________________________________________________________________________
240207: Using nested config object defined in yaml allow safely reusing variable names like input_dir for different sections of the pipeline
		without having to worry about overwriting variables. Called in smk using syntax: config["all_section_paths"]["mapping_section"]["input_dir"] 
		and config["all_section_paths"]["variant_calling_section"]["input_dir"] (GT)

240207: I'm not sure that I have properly formatted the root directory in the yaml. I use the root: /nano_mamp_pipes/ in all paths 
		which are written (and subsequently read) by the tool. I think may need to use the `os` module to find the location of this smk file to
		properly format the working directory like this:

		# in beginning of smk file:
		root_dir = os.path.dirname(os.path.realpath(__file__)) # already implemented above

		# in yaml: #TODO
		mapping_section:
		  output_dir: "analysis/alignments"

		# in smk file - construct the path to the output directory"  #TODO
		mapping_output_dir = os.path.join(root_dir, config['mapping_section']['output_dir'])

		# in smk file - insert into rule:  #TODO
		rule example_rule:
    		output:
        		os.path.join(mapping_output_dir, "output_file.txt")
_______________________________________________________________________________________________________________________________________

240207: Instead of splitting this smk and yaml into sections,  we may want to look into using mutiple smk files and calling them in 
		this centralized smk using `include` or `subflow` to make the rule execution less complicated and code more 
		modular(readable/debuggable/etc) (GT)
_______________________________________________________________________________________________________________________________________
"""

import os
import glob

root_dir = os.path.dirname(os.path.realpath(__file__))

# Function to get sample names based on FASTQ file prefixes
def get_samples():
    return [os.path.basename(f).split("_")[0] for f in glob.glob(config['all_section_paths']['mapping_section']["path_to_fastq"] + "/*_R1*.fastq.gz")]

# Create directories if they don't exist
os.makedirs(config["all_section_paths"]["mapping_section"]["temp_dir"], exist_ok=True)
os.makedirs(config["all_section_paths"]["mapping_section"]["output_dir"], exist_ok=True)

# Rule all to specify final output files
rule all:
    input:
        expand(os.path.join(config["all_section_paths"]["mapping_section"]["output_dir"], "{sample}.sorted.bam.bai"), sample=get_samples()),
		merged_gvcf_file=config["all_section_paths"]["variant_calling_section"]["output_dir"]+'/bcftools_merge_output_combined/nanomip_variants.gvcf.gz',
		merged_vcf_file=config["all_section_paths"]["variant_calling_section"]["output_dir"]+'/bcftools_merge_output_combined/nanomip_variants.vcf.gz'


"""
_______________________________________________________________________________________________________________________________________

Section 1) guppy_basecalling
_______________________________________________________________________________________________________________________________________
"""

"""
_______________________________________________________________________________________________________________________________________

Section 2) cat_fastqs:
_______________________________________________________________________________________________________________________________________
"""

"""
_______________________________________________________________________________________________________________________________________

Section 3) get_barcode_numbers
_______________________________________________________________________________________________________________________________________
"""

"""
_______________________________________________________________________________________________________________________________________

Section 4) get_barcode_seqs
_______________________________________________________________________________________________________________________________________
"""

"""
_______________________________________________________________________________________________________________________________________

Section 5) demux_samples
_______________________________________________________________________________________________________________________________________
"""

"""
_______________________________________________________________________________________________________________________________________

Section 6) mapping_section
_______________________________________________________________________________________________________________________________________

"""

rule map_to_main_genome:
	"""
	Rule for aligning reads with minimap2 to .paf files
    """
	input:
		indexed_genome=genome_location+'/genome.fa',
		sample_fastq=config['all_section_paths']['mapping_section']['fastq_dir'] + '/{sample}.fastq.gz'
	output:
		sample_paf=config["all_section_paths"]["mapping_section"]["output_dir"]+'/paf_main_files/{sample}.paf'
	shell:
		'minimap2 -ax sr {input.indexed_genome} {input.sample_fastq} -o {output.sample_paf}'

rule make_sam:
	"""
	Rule for aligning reads with minimap2 to sam files
    """
	input:
		indexed_genome=genome_location+'/genome.fa',
		sample_fastq=config['all_section_paths']['mapping_section']['fastq_dir']+'/{sample}.fastq.gz'
	output:
		sample_sam=config["all_section_paths"]["mapping_section"]["output_dir"]+'/sam_main_files/{sample}.sam'
	shell:
		'minimap2 -ax sr {input.indexed_genome} {input.sample_fastq} -o {output.sample_sam}'

rule make_bam:
	"""
	Convert sam to bam with samtools view -b
    """
	input:
		sample_sam=config["all_section_paths"]["mapping_section"]["output_dir"]+'/sam_main_files/{sample}.sam'
	output:
		sample_bam=temp(config["all_section_paths"]["mapping_section"]["temp_dir"]+'/bam_main_files/{sample}.bam')
	shell:
		'samtools view -b -o {output.sample_bam} {input.sample_sam}'

rule sort_bam:
	"""
	Sort bams by coordinate with samtools sort
    """
	input:
		sample_bam=config["all_section_paths"]["mapping_section"]["output_dir"]+'/bam_main_files/{sample}.bam'
	output:
		sorted_bam=config["all_section_paths"]["mapping_section"]["output_dir"]+'/bam_main_files/{sample}_sorted.bam'
	shell:
		'samtools sort -o {output.sorted_bam} {input.sample_bam}'

rule index_bam:
	"""
	index bams with samtools index
    """
	input:
		sorted_bam=config["all_section_paths"]["mapping_section"]["output_dir"]+'/bam_main_files/{sample}_sorted.bam'
	output:
		bam_index=config["all_section_paths"]["mapping_section"]["output_dir"]+'/bam_main_files/{sample}_sorted.bam.bai'
	shell:
		'samtools index {input.sorted_bam}'

rule merge_bams:
	"""
	merge bams from all samples into one merged bam file for exploratory analysis by examining on UCSC cluster or for other summary purposes
    """
	input:
		bam_indices=expand(config["all_section_paths"]["mapping_section"]["output_dir"]+'/bam_main_files/{sample}_sorted.bam.bai', sample=config['samples']),
		sorted_bams=expand(config["all_section_paths"]["mapping_section"]["output_dir"]+'/bam_main_files/{sample}_sorted.bam', sample=config['samples'])
	output:
		merged_bams=temp(config["all_section_paths"]["mapping_section"]["temp_dir"]+'/overall_summary/merged.bam')
	shell:
		'samtools merge {output.merged_bams} {input.sorted_bams}'

rule sort_merged:
	"""
	sort the merged "all sample" bam file by coordinate
    """
	input:
		merged_bams=config["all_section_paths"]["mapping_section"]["output_dir"]+'/overall_summary/merged.bam'
	output:
		sorted_merged_bam=config["all_section_paths"]["mapping_section"]["output_dir"]+'/overall_summary/merged_sorted.bam'
	shell:
		'samtools sort -o {output.sorted_merged_bam} {input.merged_bams}'

rule index_merged:
	"""
	index the merged "all sample" bam file
    """
	input:
		sorted_bam=config["all_section_paths"]["mapping_section"]["output_dir"]+'/overall_summary/merged_sorted.bam'
	output:
		bam_index=config["all_section_paths"]["mapping_section"]["output_dir"]+'/overall_summary/merged_sorted.bam.bai'
	shell:
		'samtools index {input.sorted_bam}'

"""
_______________________________________________________________________________________________________________________________________

# Section 7) variant_calling_section
_______________________________________________________________________________________________________________________________________

"""

rule run_clair3_raw:
	"""
	call variants in all regions in sorted bams with clair3 and output gvcf and vcf files. Using 0.0001 as the minimum allele frequency for both snps and indels.
    """
	input:
		samples_to_run=config["all_section_paths"]["variant_calling_section"]["input_dir"]+'/{sample}_sorted.bam',
		REF = config["all_section_paths"]["variant_calling_section"]['reference_dir']+'/genome.fa', # make resources
		sif_path = config["all_section_paths"]["variant_calling_section"]['sif_path']
	output:
		vcf=config["all_section_paths"]["variant_calling_section"]['output_dir']+'/{sample}_raw'+'/merge_output.vcf.gz',
		gvcf=config["all_section_paths"]["variant_calling_section"]['output_dir']+'/{sample}_raw'+'/merge_output.gvcf.gz'
	params:
		INPUT_DIR = config["all_section_paths"]["variant_calling_section"]['input_dir'],
		OUTPUT_DIR = config["all_section_paths"]["variant_calling_section"]['output_dir']+'/{sample}_raw',
		REF_DIR = config["all_section_paths"]["variant_calling_section"]['reference_dir'],
		RESOURCE_DIR = config["all_section_paths"]["variant_calling_section"]['resource_dir'],
		THREADS = 4,
		MODEL_NAME = 'r941_prom_sup_g5014',
		targets_vcf = config["all_section_paths"]["variant_calling_section"]['resource_dir']+'/targets.vcf.gz'
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
			--snp_min_af=0.0001 \
			--indel_min_af=0.0001 \
		  	--include_all_ctgs \
			--threads=1 \
			--print_ref_calls
		  '''

rule run_clair3_targets:
	"""
	call variants in targetted sites in same sorted bams with clair3 and output gvcf and vcf files. Must use 0 as the minimum allele frequency for both snps and indels to report coverage of reference calls in targetted sites (Ex. to report that there are no variants in pfk13 622I and we have X coverage to prove it.)
    """
	input:
		samples_to_run=config["all_section_paths"]["variant_calling_section"]['input_dir']+'/{sample}_sorted.bam',
		REF = config["all_section_paths"]["variant_calling_section"]['reference_dir']+'/genome.fa', # make resources
		sif_path = config["all_section_paths"]["variant_calling_section"]['sif_path']
	output:
		vcf=config["all_section_paths"]["variant_calling_section"]['output_dir']+'/{sample}_targets'+'/merge_output.vcf.gz',
		gvcf=config["all_section_paths"]["variant_calling_section"]['output_dir']+'/{sample}_targets'+'/merge_output.gvcf.gz'
	params:
		INPUT_DIR = config["all_section_paths"]["variant_calling_section"]['input_dir'],
		OUTPUT_DIR = config["all_section_paths"]["variant_calling_section"]['output_dir']+'/{sample}_targets',
		REF_DIR = config["all_section_paths"]["variant_calling_section"]['reference_dir'],
		RESOURCE_DIR = config["all_section_paths"]["variant_calling_section"]['resource_dir'],
		THREADS = 4,
		MODEL_NAME = 'r941_prom_sup_g5014',
		targets_vcf = config["all_section_paths"]["variant_calling_section"]['resource_dir']+'/targets.vcf.gz'
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
		original_gvcfs_to_index_raw=config["all_section_paths"]["variant_calling_section"]['output_dir']+'/{sample}_raw'+'/merge_output.gvcf.gz'
	output:
		index_file_gvcf_raw=config["all_section_paths"]["variant_calling_section"]['output_dir']+'/{sample}_raw'+'/merge_output.gvcf.gz.csi'
	shell:
		'''
		bcftools index {input.original_gvcfs_to_index_raw}
		'''

rule bcftools_index_gvcf_targets:
	input:
		original_gvcfs_to_index_targets=config["all_section_paths"]["variant_calling_section"]['output_dir']+'/{sample}_targets'+'/merge_output.gvcf.gz'
	output:
		index_file_gvcf_targets=config["all_section_paths"]["variant_calling_section"]['output_dir']+'/{sample}_targets'+'/merge_output.gvcf.gz.csi'
	shell:
		'''
		bcftools index {input.original_gvcfs_to_index_targets}
		'''

rule bcftools_index_vcf_raw:
	input:
		original_vcfs_to_index_raw=config["all_section_paths"]["variant_calling_section"]['output_dir']+'/{sample}_raw'+'/merge_output.vcf.gz'
	output:
		index_file_vcf_raw=config["all_section_paths"]["variant_calling_section"]['output_dir']+'/{sample}_raw'+'/merge_output.vcf.gz.csi'
	shell:
		'''
		bcftools index {input.original_vcfs_to_index_raw}
		'''

rule bcftools_index_vcf_targets:
	input:
		original_vcfs_to_index_targets=config["all_section_paths"]["variant_calling_section"]['output_dir']+'/{sample}_targets'+'/merge_output.vcf.gz'
	output:
		index_file_vcf_targets=config["all_section_paths"]["variant_calling_section"]['output_dir']+'/{sample}_targets'+'/merge_output.vcf.gz.csi'
	shell:
		'''
		bcftools index {input.original_vcfs_to_index_targets}
		'''

rule bcftools_merge_raw_gvcf:
	input:
		original_gvcfs_raw=expand(config["all_section_paths"]["variant_calling_section"]['output_dir']+'/{sample}_raw'+'/merge_output.gvcf.gz.csi',sample=config["all_section_paths"]["variant_calling_section"]['samples'])
	output:
		merged_gvcf_file_raw=config["all_section_paths"]["variant_calling_section"]['output_dir']+'/bcftools_merge_output_raw/merged_multisample.gvcf.gz'
	shell:
		'''
		bcftools merge --file-list /nfs/jbailey5/baileyweb/gtollefs/mips_to_nanopore_gt/bams_to_vcf/resources/gvcf_list_raw.txt --force-samples -O z -o {output.merged_gvcf_file_raw}
		'''

rule bcftools_merge_targets_gvcf:
	input:
		original_gvcfs_targets=expand(config["all_section_paths"]["variant_calling_section"]['output_dir']+'/{sample}_targets'+'/merge_output.gvcf.gz.csi',sample=config["all_section_paths"]["variant_calling_section"]['samples'])
	output:
		merged_gvcf_file_targets=config["all_section_paths"]["variant_calling_section"]['output_dir']+'/bcftools_merge_output_targets/merged_multisample.gvcf.gz'
	shell:
		'''
		bcftools merge --file-list /nfs/jbailey5/baileyweb/gtollefs/mips_to_nanopore_gt/bams_to_vcf/resources/gvcf_list_targets.txt --force-samples -O z -o {output.merged_gvcf_file_targets}
		'''

rule bcftools_merge_raw_vcf:
	input:
		original_vcfs=expand(config["all_section_paths"]["variant_calling_section"]['output_dir']+'/{sample}_raw'+'/merge_output.vcf.gz.csi',sample=config["all_section_paths"]["variant_calling_section"]['samples'])
	output:
		merged_vcf_file=config["all_section_paths"]["variant_calling_section"]['output_dir']+'/bcftools_merge_output_raw/merged_multisample.vcf.gz'
	shell:
		'''
		bcftools merge --file-list /nfs/jbailey5/baileyweb/gtollefs/mips_to_nanopore_gt/bams_to_vcf/resources/vcf_list_raw.txt --force-samples -O z -o {output.merged_vcf_file}
		'''

rule bcftools_merge_targets_vcf:
	input:
		original_vcfs=expand(config["all_section_paths"]["variant_calling_section"]['output_dir']+'/{sample}_targets'+'/merge_output.vcf.gz.csi',sample=config["all_section_paths"]["variant_calling_section"]['samples'])
	output:
		merged_vcf_file=config["all_section_paths"]["variant_calling_section"]['output_dir']+'/bcftools_merge_output_targets/merged_multisample.vcf.gz'
	shell:
		'''
		bcftools merge --file-list /nfs/jbailey5/baileyweb/gtollefs/mips_to_nanopore_gt/bams_to_vcf/resources/vcf_list_targets.txt --force-samples -O z -o {output.merged_vcf_file}
		'''

rule bcftools_index_merged_raw_gvcf:
	input:
		merged_gvcfs_to_index_raw=config["all_section_paths"]["variant_calling_section"]['output_dir']+'/bcftools_merge_output_raw/merged_multisample.gvcf.gz'
	output:
		index_file_gvcf_targets=config["all_section_paths"]["variant_calling_section"]['output_dir']+'/bcftools_merge_output_raw/merged_multisample.gvcf.gz.csi'
	shell:
		'''
		bcftools index {input.merged_gvcfs_to_index_raw}
		'''

rule bcftools_index_merged_targets_gvcf:
	input:
		merged_gvcfs_to_index_targets=config["all_section_paths"]["variant_calling_section"]['output_dir']+'/bcftools_merge_output_targets/merged_multisample.gvcf.gz'
	output:
		index_file_gvcf_targets=config["all_section_paths"]["variant_calling_section"]['output_dir']+'/bcftools_merge_output_targets/merged_multisample.gvcf.gz.csi'
	shell:
		'''
		bcftools index {input.merged_gvcfs_to_index_targets}
		'''

rule bcftools_index_merged_raw_vcf:
	input:
		merged_vcfs_to_index_raw=config["all_section_paths"]["variant_calling_section"]['output_dir']+'/bcftools_merge_output_raw/merged_multisample.vcf.gz'
	output:
		index_file_vcf_targets=config["all_section_paths"]["variant_calling_section"]['output_dir']+'/bcftools_merge_output_raw/merged_multisample.vcf.gz.csi'
	shell:
		'''
		bcftools index {input.merged_vcfs_to_index_raw}
		'''

rule bcftools_index_merged_targets_vcf:
	input:
		merged_vcfs_to_index_targets=config["all_section_paths"]["variant_calling_section"]['output_dir']+'/bcftools_merge_output_targets/merged_multisample.vcf.gz'
	output:
		index_file_vcf_targets=config["all_section_paths"]["variant_calling_section"]['output_dir']+'/bcftools_merge_output_targets/merged_multisample.vcf.gz.csi'
	shell:
		'''
		bcftools index {input.merged_vcfs_to_index_targets}
		'''

rule bcftools_merge_raw_and_targets_gvcf:
	input:
		merged_raw_gvcf=config["all_section_paths"]["variant_calling_section"]['output_dir']+'/bcftools_merge_output_raw/merged_multisample.gvcf.gz',
		merged_targets_gvcf=config["all_section_paths"]["variant_calling_section"]['output_dir']+'/bcftools_merge_output_targets/merged_multisample.gvcf.gz',
		index_raw_gvcf=config["all_section_paths"]["variant_calling_section"]['output_dir']+'/bcftools_merge_output_raw/merged_multisample.gvcf.gz.csi',
		index_targets_gvcf=config["all_section_paths"]["variant_calling_section"]['output_dir']+'/bcftools_merge_output_targets/merged_multisample.gvcf.gz.csi'
	output:
		merged_raw_and_targets_gvcf_file=config["all_section_paths"]["variant_calling_section"]['output_dir']+'/bcftools_merge_output_combined/nanomip_variants.gvcf.gz'
	shell:
		'''
		bcftools merge --force-samples -O z -o {output.merged_raw_and_targets_gvcf_file} {input.merged_raw_gvcf} {input.merged_targets_gvcf}
		'''

rule bcftools_merge_raw_and_targets_vcf:
	input:
		merged_raw_vcf=config["all_section_paths"]["variant_calling_section"]['output_dir']+'/bcftools_merge_output_raw/merged_multisample.vcf.gz',
		merged_targets_vcf=config["all_section_paths"]["variant_calling_section"]['output_dir']+'/bcftools_merge_output_targets/merged_multisample.vcf.gz',
		index_raw_vcf=config["all_section_paths"]["variant_calling_section"]['output_dir']+'/bcftools_merge_output_raw/merged_multisample.vcf.gz.csi',
		index_targets_vcf=config["all_section_paths"]["variant_calling_section"]['output_dir']+'/bcftools_merge_output_targets/merged_multisample.vcf.gz.csi'
	output:
		merged_raw_and_targets_vcf_file=config["all_section_paths"]["variant_calling_section"]['output_dir']+'/bcftools_merge_output_combined/nanomip_variants.vcf.gz'
	shell:
		'''
		bcftools merge --force-samples -O z -o {output.merged_raw_and_targets_vcf_file} {input.merged_raw_vcf} {input.merged_targets_vcf}
		'''
"""
_______________________________________________________________________________________________________________________________________

# Section 8) vcf_to_tables_section
_______________________________________________________________________________________________________________________________________
"""



"""
_______________________________________________________________________________________________________________________________________

# Section 9) miplithon_and_mamplicoRn_section
_______________________________________________________________________________________________________________________________________
"""



"""
_______________________________________________________________________________________________________________________________________

# Section 9) generate_count_summaries_section
_______________________________________________________________________________________________________________________________________
"""

# prior contents - commented out to make this a fully functional smk file for now until we can get the rest of the pipeline working to make this section work:
# rule all:
# 	input:
# 		count_summary='count_summary.txt'

# rule generate_count_summaries:
# 	input:
# 		mip_info=config['project_resources']+'/mip_ids/mip_info.json',
# 		sam_files=expand(config['sam_folder']+'/{sample}.sam', sample=all_samples)
# 	output:
# 		count_summary='count_summary.txt'
# 	script:
# 		'generate_count_summaries.py'

# rule map_og_reads:
