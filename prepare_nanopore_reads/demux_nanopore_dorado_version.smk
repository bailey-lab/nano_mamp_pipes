configfile: 'demux_nanopore_dorado_version.yaml'

rule all:
	input:
		demuxed_fastq_folder=config['output_folder']+'/demuxed_fastq'

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