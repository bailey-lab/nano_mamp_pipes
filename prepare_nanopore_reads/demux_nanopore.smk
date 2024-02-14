configfile: 'mips_to_nanopore.yaml'

rule run_basecaller:
	'''
	TODO: make it so this script only 'completes' when every input fast5 has a
	corresponding fastq
	'''
	input:
		fast5_directory=config['fast5_dir']
	params:
		basecalled_fastq=config['fastq_dir']
	output:
		complete_status='basecalling_completed.txt'
	shell:
		'''
		sbatch scripts/guppy_basecaller.sh {input.fast5_directory} {params.basecalled_fastq}
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

rule demux_files:
	input:
		zipped_file=config['catted_fastq']+'.gz',
		sample_barcodes=config['sample_barcodes']
	output:
		demuxed_fastq_folder=config['demuxed_folder']
