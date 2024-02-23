configfile: 'demux_nanopore.yaml'

rule all:
	input:
		zipped_file=config['catted_fastq']+'.gz'

rule run_basecaller:
	'''
	TODO: make it so this script only 'completes' when every input fast5 has a
	corresponding fastq
	'''
	input:
		fast5_directory=config['fast5_dir']
	params:
		basecalled_fastq=config['fastq_dir'],
		config='dna_r10.4.1_e8.2_400bps_sup.cfg'
	output:
		complete_status='basecalling_completed.txt',
		basecalled_fastq=config['fastq_dir']+'/pass'
	resources:
		platform='gpu',
		mem_mb=48000,
		nodes=8,
		time_min=2880
	shell:
		'''
		bash scripts/guppy_basecaller2.sh {params.config} {input.fast5_directory} {params.basecalled_fastq}
		touch basecalling_completed.txt
		'''

rule cat_files:
	input:
		basecalled_fastq=config['fastq_dir']+'/pass'
	output:
		catted_unzipped_file=temp(config['catted_fastq']),
		zipped_file=config['catted_fastq']+'.gz'
	script:
		'''
		scripts/alternative_cat.py
		'''
'''
rule demux_files:
	input:
		zipped_file=config['catted_fastq']+'.gz',
		sample_barcodes=config['sample_barcodes']
	output:
		demuxed_fastq_folder=config['demuxed_folder']
'''