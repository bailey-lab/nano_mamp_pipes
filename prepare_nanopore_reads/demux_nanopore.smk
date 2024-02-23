configfile: 'demux_nanopore.yaml'

rule run_basecaller:
	'''
	TODO: make it so this script only 'completes' when every input fast5 has a
	corresponding fastq.
	GT: Do we want to add a sample prefix pull from fast5 and check in fastq for completeness and auto rerunning of guppy on failed samples, or does the guppy_basecaller.sh handle this already? 
	#  ^ This way we could make fastq files with prefix the expected output from the rule and input for the next rule.
	#  ^ may be able to remove the squeue and scontrol checks if expeected fastq output samples are known
	GT: Currently rule should run as long as sbatch job is running now with the implemented scontrol/squeue check check
	'''
	input:
		fast5_directory=config['fast5_dir']
	params:
		basecalled_fastq=config['fastq_dir']
	output:
		complete_status=config['fastq_dir']+'basecalling_completed.txt' # maybe don't want to write this to fastq directory since prefixes are pulled from fastq directory. Also should we put these in a temp or separate non-demuxed fastq directory and make a fastq directory for demuxed fastq with real sample prefixes for reading for mapping rules?
	shell:
		'''
		jobid=$(sbatch --parsable scripts/guppy_basecaller.sh {input.fast5_directory} {params.basecalled_fastq})
		while : ; do 
			if ! squeue -j $jobid &>/dev/null; then
			# use scontrol to see if the job is no longer present to avoid unwanted error signal to snakemake
				if scontrol show job $jobid &>/dev/null; then
					# if job is still on slurm keep checking
					sleep 30
				else
					# if job no longer seen on SLURM, break to prevent error signal
					break
				fi
			fi
			sleep 30
		done
		touch {output.complete_status}
		'''

rule cat_files:
	input:
		basecalled_fastq=config['fastq_dir']+'/pass',
		complete_status=config['fastq_dir']+'basecalling_completed.txt' # add guppy finish check before running this rule - unless we implement {sample}fastq prefixes as input
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
