configfile: 'demux_nanopore.yaml'

rule all:
	input:
		demuxed_fastq_folder=config['output_folder']+'/demuxed_fastq'

rule run_basecaller:
	'''
	TODO: make it so this script only 'completes' (so that basecalling_
	completed.txt is only touched) when every input fast5 has a	corresponding
	fastq
	'''
	input:
		fast5_directory=config['fast5_dir']
	params:
		basecalled_fastq=config['output_folder']+'/basecalling',
		config='dna_r10.4.1_e8.2_400bps_sup.cfg'
	output:
		complete_status='basecalling_completed.txt',
		basecalled_fastq=directory(config['output_folder']+'/basecalling/pass')
	resources:
		platform='gpu',
		mem_mb=48000,
		nodes=8,
		time_min=2879,
		extra_flags='--gres=gpu:1 --gres-flags=enforce-binding -N 1'
	shell:
		'''
		bash scripts/guppy_basecaller2.sh {params.config} {input.fast5_directory} {params.basecalled_fastq}
		touch basecalling_completed.txt
		'''

rule cat_files:
	input:
		basecalled_fastq=config['output_folder']+'/basecalling/pass'
	params:
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
		time_min=2879
	shell:
		'''
		singularity exec \
		-B {input.barcodes_dir}:/opt/resources \
		-B {params.fastq_dir}:/opt/demuxing {input.elucidator_sif} \
		elucidator extractByIlluminaAaptors --fastqgz /opt/demuxing/catted_pass.fastq.gz \
		--overWriteDir --dout {output.demuxed_fastq_folder} \
		--illuminaBarcodeSampleSheet /opt/resources/{params.sample_barcodes}
		'''