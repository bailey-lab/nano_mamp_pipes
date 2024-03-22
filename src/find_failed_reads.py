'''
counts the number of reads that failed during:
1. basecalling (reads that didn't pass quality score cutoffs)
2. demultiplexing (reads that had no sample barcode assigned to them)
3. amplicon assignment (reads that didn't have correct regional primers associated
4. seekdeep quality thresholds (too short reads, too long reads, low quality, missing one or both primers)
'''
import os
import gzip

fail_path='/nfs/jbailey5/baileyweb/asimkin/mips_to_nanopore/240215_MtoN_P2/guppy_version/prepped_reads/basecalling/fail'
success_path='/nfs/jbailey5/baileyweb/asimkin/mips_to_nanopore/240215_MtoN_P2/guppy_version/prepped_reads/basecalling/pass'
undetermined='/nfs/jbailey5/baileyweb/asimkin/mips_to_nanopore/240215_MtoN_P2/guppy_version/prepped_reads/demuxed_fastq/fastq/Undetermined.fastq.gz'
all_extracted='/nfs/jbailey5/baileyweb/asimkin/mips_to_nanopore/240215_MtoN_P2/guppy_version/prepped_reads/demuxed_fastq/fastq'

def get_demux_fail(undetermined):
	read_counter=0
	for line_number, line in enumerate(gzip.open(undetermined)):
		if line_number%4==0:
			read_counter+=1
	return read_counter

def get_total_reads(read_folder):
	files=os.listdir(read_folder)
	read_counter=0
	for file_name in files:
		for line_number, line in enumerate(gzip.open(read_folder+'/'+file_name)):
			if line_number%4==0:
				read_counter+=1
	return read_counter


print('basecalling fails is', get_total_reads(fail_path))
print('basecalling successes is', get_total_reads(success_path))
all_extracted=get_total_reads(all_extracted)
failed_demux=get_demux_fail(undetermined)
print('Undetermined is', failed_demux)
print('all extracted is', all_extracted)
print('successfully demuxed is', all_extracted-failed_demux)
