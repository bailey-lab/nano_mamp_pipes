'''
this is a troubleshooting script to figure out how much (if any) space is saved
by unzipping files, catting them, and gzipping all at once rather than catting
gzipped files. Theoretically, unzipping, catting, and gzipping should be a bit
more efficient because similar regions between different input files can be
compressed together instead of repeatedly compressed individually.
'''
import os
import gzip
import subprocess

input_dir=snakemake.params.basecalled_fastq
output_file_name=snakemake.params.catted_unzipped_file

output_file=open(output_file_name, 'w')
for file_name in os.listdir(input_dir):
	if file_name.endswith('.fastq.gz'):
		for line in gzip.open(input_dir+'/'+file_name, mode='rt'):
			output_file.write(line)
output_file.close()
subprocess.call(f'gzip {output_file_name}', shell=True)
