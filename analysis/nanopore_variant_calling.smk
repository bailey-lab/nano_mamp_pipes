configfile: 'mips_to_nanopore.yaml'

import os
import glob

# Get all sample names
def get_samples():
    fastq_files = glob.glob(os.path.join(config['fastq_dir'], '*.fastq.gz'))
    sample_names = [
        os.path.basename(f).replace('.fastq.gz', '') for f in fastq_files
        if 'Undetermined' not in f and 'NTC' not in os.path.basename(f)
    ]
    return sample_names

# Rule all to check final files
rule all:
    input:
        finished_gvcf=config['output_dir']+'/variant_calling/merged_multisample.gvcf.gz',
        finished_vcf=config['output_dir']+'/variant_calling/merged_multisample.vcf.gz'

# Run Clair3 variant calling for each sample
rule run_clair3_raw:
    input:
        bam=config["output_dir"]+'/mapping/bam_files/{sample}.sorted.bam',
        REF=config['genome_directory']+'/genome.fa',
        sif_path=config["clair3_sif_path"]
    output:
        vcf=config["output_dir"]+'/variant_calling/{sample}_raw/merge_output.vcf.gz',
        gvcf=config["output_dir"]+'/variant_calling/{sample}_raw/merge_output.gvcf.gz'
    params:
        INPUT_DIR=config["output_dir"]+'/mapping/bam_files',
        OUTPUT_DIR=config["output_dir"]+'/variant_calling/{sample}_raw',
        GENOME_DIR=config['genome_directory'],
        MODEL_NAME='r1041_e82_400bps_sup_v430_bacteria_finetuned',
        THREADS=4
    shell:
        '''
        echo 'Processing sample: {wildcards.sample}'
        singularity exec \
            -B {params.INPUT_DIR} \
            -B {params.OUTPUT_DIR} \
            -B {params.GENOME_DIR} \
            {input.sif_path} \
            /opt/bin/run_clair3.sh \
            --bam_fn={input.bam} \
            --ref_fn={input.REF} \
            --model_path=/opt/models/{params.MODEL_NAME} \
            --output={params.OUTPUT_DIR} \
            --threads={params.THREADS} \
            --platform=ont \
            --include_all_ctgs \
            --no_phasing_for_fa \
            --sample_name={wildcards.sample} \
            --gvcf \
            --snp_min_af=0.001 \
            --indel_min_af=0.001 \
            --print_ref_calls
        '''

# Index GVCFs and VCFs
rule index_gvcf:
    input:
        gvcf=config["output_dir"]+'/variant_calling/{sample}_raw/merge_output.gvcf.gz'
    output:
        gvcf_index=config["output_dir"]+'/variant_calling/{sample}_raw/merge_output.gvcf.gz.csi'
    shell:
        '''
        module load bcftools;
        bcftools index {input.gvcf}
        '''

rule index_vcf:
    input:
        vcf=config["output_dir"]+'/variant_calling/{sample}_raw/merge_output.vcf.gz'
    output:
        vcf_index=config["output_dir"]+'/variant_calling/{sample}_raw/merge_output.vcf.gz.csi'
    shell:
        '''
        module load bcftools;
        bcftools index {input.vcf}
        '''

# Merge GVCFs
rule merge_gvcfs:
    input:
        gvcf_indices=expand(config["output_dir"]+'/variant_calling/{sample}_raw/merge_output.gvcf.gz.csi', sample=get_samples())
    output:
        merged_gvcf=config["output_dir"]+'/variant_calling/merged_multisample.gvcf.gz'
    shell:
        '''
        module load bcftools;
        bcftools merge \
            {config[output_dir]}/variant_calling/*_raw/merge_output.gvcf.gz \
            --force-samples -O z -o {output.merged_gvcf}
        '''

# Merge VCFs
rule merge_vcfs:
    input:
        vcf_indices=expand(config["output_dir"]+'/variant_calling/{sample}_raw/merge_output.vcf.gz.csi', sample=get_samples())
    output:
        merged_vcf=config["output_dir"]+'/variant_calling/merged_multisample.vcf.gz'
    shell:
        '''
        module load bcftools;
        bcftools merge \
            {config[output_dir]}/variant_calling/*_raw/merge_output.vcf.gz \
            --force-samples -O z -o {output.merged_vcf}
        '''
