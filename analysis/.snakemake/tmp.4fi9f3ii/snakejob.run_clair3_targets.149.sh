#!/bin/sh
# properties = {"type": "single", "rule": "run_clair3_targets", "local": false, "input": ["/nfs/jbailey5/baileyweb/gtollefs/mips_to_nanopore_gt/240215_MtoN_P2/second_output_rerun_for_completeness_check/mapping/bam_files/3d7-4000.sorted.bam", "/nfs/jbailey5/baileyweb/bailey_share/resources/MIP_species_resources/pf/Pf_3D7/genomes/genome.fa", "/nfs/jbailey5/baileyweb/asimkin/aapter_test/run_230508_mipstoNanopore_P2/bam_to_vcf/clair3_latest.sif"], "output": ["/nfs/jbailey5/baileyweb/gtollefs/mips_to_nanopore_gt/240215_MtoN_P2/second_output_rerun_for_completeness_check/variant_calling/3d7-4000_targets/merge_output.vcf.gz", "/nfs/jbailey5/baileyweb/gtollefs/mips_to_nanopore_gt/240215_MtoN_P2/second_output_rerun_for_completeness_check/variant_calling/3d7-4000_targets/merge_output.gvcf.gz"], "wildcards": {"samples_for_targeted_calling": "3d7-4000"}, "params": {"INPUT_DIR": "/nfs/jbailey5/baileyweb/gtollefs/mips_to_nanopore_gt/240215_MtoN_P2/second_output_rerun_for_completeness_check/mapping/bam_files", "OUTPUT_DIR": "/nfs/jbailey5/baileyweb/gtollefs/mips_to_nanopore_gt/240215_MtoN_P2/second_output_rerun_for_completeness_check/variant_calling/3d7-4000_targets", "REF_DIR": "/nfs/jbailey5/baileyweb/bailey_share/resources/MIP_species_resources/pf/Pf_3D7/genomes", "RESOURCE_DIR": "/nfs/jbailey5/baileyweb/gtollefs/mips_to_nanopore_gt/bams_to_vcf/resources", "THREADS": 4, "MODEL_NAME": "r941_prom_sup_g5014", "targets_vcf": "/nfs/jbailey5/baileyweb/gtollefs/mips_to_nanopore_gt/bams_to_vcf/resources/targets.vcf.gz"}, "log": [], "threads": 1, "resources": {"mem_mb": 4000, "mem_mib": 3815, "disk_mb": 1819, "disk_mib": 1735, "tmpdir": "<TBD>", "nodes": 1, "time_min": 120}, "jobid": 149, "cluster": {}}
cd '/nfs/jbailey5/baileyweb/gtollefs/mips_to_nanopore_gt/nano_mamp_pipes/analysis' && /oscar/home/gtollefs/mambaforge/envs/snakemake/bin/python3.11 -m snakemake --snakefile '/nfs/jbailey5/baileyweb/gtollefs/mips_to_nanopore_gt/nano_mamp_pipes/analysis/mips_to_nanopore_variant_calling_with_target_check.smk' --target-jobs 'run_clair3_targets:samples_for_targeted_calling=3d7-4000' --allowed-rules 'run_clair3_targets' --cores 'all' --attempt 1 --force-use-threads  --resources 'mem_mb=4000' 'mem_mib=3815' 'disk_mb=1819' 'disk_mib=1735' 'nodes=1' 'time_min=120' --wait-for-files '/nfs/jbailey5/baileyweb/gtollefs/mips_to_nanopore_gt/nano_mamp_pipes/analysis/.snakemake/tmp.4fi9f3ii' '/nfs/jbailey5/baileyweb/gtollefs/mips_to_nanopore_gt/240215_MtoN_P2/second_output_rerun_for_completeness_check/mapping/bam_files/3d7-4000.sorted.bam' '/nfs/jbailey5/baileyweb/bailey_share/resources/MIP_species_resources/pf/Pf_3D7/genomes/genome.fa' '/nfs/jbailey5/baileyweb/asimkin/aapter_test/run_230508_mipstoNanopore_P2/bam_to_vcf/clair3_latest.sif' --force --keep-target-files --keep-remote --max-inventory-time 0 --nocolor --notemp --no-hooks --nolock --ignore-incomplete --rerun-triggers 'software-env' 'input' 'mtime' 'params' 'code' --skip-script-cleanup  --conda-frontend 'mamba' --wrapper-prefix 'https://github.com/snakemake/snakemake-wrappers/raw/' --latency-wait 60 --scheduler 'ilp' --scheduler-solver-path '/oscar/home/gtollefs/mambaforge/envs/snakemake/bin' --default-resources 'mem_mb=4000' 'disk_mb=max(2*input.size_mb, 1000)' 'tmpdir=system_tmpdir' 'nodes=1' 'time_min=120' --mode 2 && touch '/nfs/jbailey5/baileyweb/gtollefs/mips_to_nanopore_gt/nano_mamp_pipes/analysis/.snakemake/tmp.4fi9f3ii/149.jobfinished' || (touch '/nfs/jbailey5/baileyweb/gtollefs/mips_to_nanopore_gt/nano_mamp_pipes/analysis/.snakemake/tmp.4fi9f3ii/149.jobfailed'; exit 1)

