configfile: 'mips_to_nanopore.yaml'
#root_folder='/nfs/jbailey5/baileyweb/asimkin/mips_to_nanopore'

def get_all_samples(sam_folder):
	import os
	sample_names=[]
	sam_paths=os.listdir(sam_folder)
	for sample in sam_paths:
		sample_names.append(sample[:-4])
	return sample_names

all_samples=get_all_samples(config['sam_folder'])

'''
rule guppy_basecalling:

rule cat_fastqs:

rule get_barcode_numbers:

rule get_barcode_seqs:

rule demux_samples:


rule map_demuxed_reads:
	input:
		fastq_folder=root_folder+'/run_230508_mipstoNanopore_P2/extractedByIlluminaBarcodes/fastq'
	output:
		mapped_reads=root_folder+'/snakemake_output/mapped_reads'
	script:
		'scripts/map_demuxed_reads.py'
'''

rule all:
	input:
		count_summary='count_summary.txt'

rule generate_count_summaries:
	input:
		mip_info=config['project_resources']+'/mip_ids/mip_info.json',
		sam_files=expand(config['sam_folder']+'/{sample}.sam', sample=all_samples)
	output:
		count_summary='count_summary.txt'
	script:
		'generate_count_summaries.py'

rule map_og_reads:
