configfile: 'demux_nanopore.yaml'

rule all:
	input:
		demuxed_fastq_folder=config['output_folder']+'/demuxed_fastq'

rule run_basecaller:
	'''
	TODO: make a version that runs each fast5 file separately, and make it so
	this script only 'completes' when every input fast5 has status 'completed'
	(current version is all or nothing - succeeds or needs to restart
	basecalling from beginning
	'''
	input:
		fast5_directory=config['fast5_dir']
	params:
		basecalled_fastq=config['output_folder']+'/basecalling',
		config=config['nanopore_config_params']
	output:
		finished_telemetry=config['output_folder']+'/basecalling/sequencing_telemetry.js'
	resources:
		partition='gpu',
		mem_mb=48000,
		cluster_threads=8,
		time_min=2879,
		extra_flags='--gres=gpu:1 --gres-flags=enforce-binding -N 1'
	shell:
		'''
		bash scripts/guppy_basecaller2.sh {params.config} {input.fast5_directory} {params.basecalled_fastq}
		'''

rule cat_files:
	input:
		finished_telemetry=config['output_folder']+'/basecalling/sequencing_telemetry.js'
	params:
		basecalled_fastq=config['output_folder']+'/basecalling/pass',
		catted_unzipped_file=config['output_folder']+'/catted_pass.fastq'
	output:
		zipped_file=temp(config['output_folder']+'/catted_pass.fastq.gz')
	resources:
		time_min=2879
	script:
		'scripts/alternative_cat.py'

rule demux_files:
	input:
		zipped_file=config['output_folder']+'/catted_pass.fastq.gz',
		elucidator_sif=config['elucidator_sif'],
		barcodes_dir=config['barcodes_dir']
	params:
		fastq_dir=config['output_folder'],
		sample_barcodes=config['sample_barcodes']
	output:
		demuxed_fastq_folder=directory(config['output_folder']+'/demuxed_fastq')
	resources:
		time_min=2879,
		mem_mb=32000
	shell:
		'''
		singularity exec \
		-B {input.barcodes_dir}:/opt/resources \
		-B {params.fastq_dir}:/opt/demuxing {input.elucidator_sif} \
		elucidator extractByIlluminaAaptors --fastqgz /opt/demuxing/catted_pass.fastq.gz \
		--overWriteDir --dout /opt/demuxing/demuxed_fastq \
		--illuminaBarcodeSampleSheet /opt/resources/{params.sample_barcodes}
		'''