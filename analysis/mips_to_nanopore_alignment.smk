configfile: 'mips_to_nanopore.yaml'

import os
import glob

# get sample prefixes
def get_samples():
    fastq_files = glob.glob(os.path.join(config['fastq_dir'], '*.fastq.gz'))
    sample_names = [
        os.path.basename(f).replace('.fastq.gz', '') for f in fastq_files
        if 'Undetermined' not in f and 'NTC' not in os.path.basename(f)
    ]
    return sample_names

# Rule all to specify final output files
rule all:
    input:
        sorted_bams=expand(config["output_dir"]+'/mapping/bam_files/{sample}.sorted.bam', sample=get_samples()),
        bam_index=expand(config["output_dir"]+'/mapping/bam_files/{sample}.sorted.bam.bai', sample=get_samples()),
        coverage_table=config['output_dir'] + "/coverage_analysis/amplicon_coverage_table.tsv"

rule make_sam:
	"""
	Rule for aligning reads with minimap2 to sam files
	"""
	input:
		indexed_genome=config['genome_directory']+'/genome.fa',
		sample_fastq=config['fastq_dir']+'/{sample}.fastq.gz'
	output:
		sample_sam=temp(config["output_dir"]+'/mapping/sam_files/{sample}.sam')
	shell:
		'minimap2 -ax sr {input.indexed_genome} {input.sample_fastq} -o {output.sample_sam}'

rule make_bam:
	"""
	Convert sam to bam with samtools view -b
	"""
	input:
		sample_sam=config["output_dir"]+'/mapping/sam_files/{sample}.sam'
	output:
		sample_bam=temp(config["output_dir"]+'/mapping/bam_files/{sample}.bam')
	shell:
		'module load samtools && samtools view -b -o {output.sample_bam} {input.sample_sam}'

rule sort_bam:
	"""
	Sort bams by coordinate with samtools sort
	"""
	input:
		sample_bam_to_sort=config["output_dir"]+'/mapping/bam_files/{sample}.bam'
	output:
		sorted_bam=config["output_dir"]+'/mapping/bam_files/{sample}.sorted.bam'
	shell:
		'module load samtools && samtools sort -o {output.sorted_bam} {input.sample_bam_to_sort}'

rule index_bam:
	"""
	index bams with samtools index
	"""
	input:
		sorted_bam=config["output_dir"]+'/mapping/bam_files/{sample}.sorted.bam'
	output:
		bam_index=config["output_dir"]+'/mapping/bam_files/{sample}.sorted.bam.bai'
	shell:
		'module load samtools && samtools index {input.sorted_bam}'

rule check_amplicon_coverage:
    input:
        bam_files=expand(config["output_dir"] + '/mapping/bam_files/{sample}.sorted.bam', sample=get_samples()),
        bam_indices=expand(config["output_dir"] + '/mapping/bam_files/{sample}.sorted.bam.bai', sample=get_samples()),
        amplicon_targets=config["amplicon_targets"]
    output:
        coverage_table=config['output_dir'] + "/coverage_analysis/amplicon_coverage_table.tsv"
    params:
        output_dir=config["output_dir"] + "/coverage_analysis"
    shell:
        """
        module load samtools && \
        mkdir -p {params.output_dir} && \
        echo -e "Sample\tAmplicon\tChrom\tStart\tEnd\tReads_Mapped\tMean_Depth\tMedian_Depth" > {output.coverage_table} && \
        for bam in {input.bam_files}; do
            sample=$(basename $bam .sorted.bam)
            while read chrom start end amplicon; do
                # Create temporary filtered BAM with reads >= 2300bp
                temp_bam="/tmp/${{sample}}_${{amplicon}}_filtered.bam"
                samtools view -h $bam "$chrom:$start-$end" | awk 'substr($1,1,1)=="@" || length($10) >= 2300' | samtools view -b > $temp_bam
                samtools index $temp_bam
                
                # Count reads and calculate mean and median depth from filtered BAM
                reads=$(samtools view -c $temp_bam "$chrom:$start-$end")
                depth_stats=$(samtools depth -r "$chrom:$start-$end" $temp_bam | awk '{{depths[NR]=$3; sum+=$3}} END {{
                    if(NR==0) {{
                        print "0\t0"
                    }} else {{
                        mean = sum/NR
                        asort(depths)
                        if(NR%2==1) {{
                            median = depths[(NR+1)/2]
                        }} else {{
                            median = (depths[NR/2] + depths[NR/2+1])/2
                        }}
                        print mean "\t" median
                    }}
                }}')
                
                echo -e "$sample\t$amplicon\t$chrom\t$start\t$end\t$reads\t$depth_stats" >> {output.coverage_table}
                
                # Clean up temporary files
                rm -f $temp_bam $temp_bam.bai
            done < {input.amplicon_targets}
        done
        """
