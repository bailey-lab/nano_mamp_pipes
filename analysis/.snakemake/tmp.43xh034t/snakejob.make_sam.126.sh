#!/bin/sh
# properties = {"type": "single", "rule": "make_sam", "local": false, "input": ["/nfs/jbailey5/baileyweb/bailey_share/resources/MIP_species_resources/pf/Pf_3D7/genomes/genome.fa", "/nfs/jbailey5/baileyweb/gtollefs/mips_to_nanopore_gt/240215_MtoN_P2/output/dorado_output/demuxed_fastq/fastq/RW-26026-078078-ArtR23RW-1.fastq.gz"], "output": ["/nfs/jbailey5/baileyweb/gtollefs/mips_to_nanopore_gt/240215_MtoN_P2/dorado_pipeline_output/mapping/sam_files/RW-26026-078078-ArtR23RW-1.sam"], "wildcards": {"sample": "RW-26026-078078-ArtR23RW-1"}, "params": {}, "log": [], "threads": 1, "resources": {"mem_mb": 4000, "mem_mib": 3815, "disk_mb": 1000, "disk_mib": 954, "tmpdir": "<TBD>", "nodes": 1, "time_min": 120}, "jobid": 126, "cluster": {}}
cd '/nfs/jbailey5/baileyweb/gtollefs/mips_to_nanopore_gt/nano_mamp_pipes/analysis' && /oscar/home/gtollefs/mambaforge/envs/snakemake/bin/python3.11 -m snakemake --snakefile '/nfs/jbailey5/baileyweb/gtollefs/mips_to_nanopore_gt/nano_mamp_pipes/analysis/mips_to_nanopore_alignment.smk' --target-jobs 'make_sam:sample=RW-26026-078078-ArtR23RW-1' --allowed-rules 'make_sam' --cores 'all' --attempt 1 --force-use-threads  --resources 'mem_mb=4000' 'mem_mib=3815' 'disk_mb=1000' 'disk_mib=954' 'nodes=1' 'time_min=120' --wait-for-files '/nfs/jbailey5/baileyweb/gtollefs/mips_to_nanopore_gt/nano_mamp_pipes/analysis/.snakemake/tmp.43xh034t' '/nfs/jbailey5/baileyweb/bailey_share/resources/MIP_species_resources/pf/Pf_3D7/genomes/genome.fa' '/nfs/jbailey5/baileyweb/gtollefs/mips_to_nanopore_gt/240215_MtoN_P2/output/dorado_output/demuxed_fastq/fastq/RW-26026-078078-ArtR23RW-1.fastq.gz' --force --keep-target-files --keep-remote --max-inventory-time 0 --nocolor --notemp --no-hooks --nolock --ignore-incomplete --rerun-triggers 'input' 'params' 'code' 'mtime' 'software-env' --skip-script-cleanup  --conda-frontend 'mamba' --wrapper-prefix 'https://github.com/snakemake/snakemake-wrappers/raw/' --latency-wait 60 --scheduler 'ilp' --scheduler-solver-path '/oscar/home/gtollefs/mambaforge/envs/snakemake/bin' --default-resources 'mem_mb=4000' 'disk_mb=max(2*input.size_mb, 1000)' 'tmpdir=system_tmpdir' 'nodes=1' 'time_min=120' --mode 2 && touch '/nfs/jbailey5/baileyweb/gtollefs/mips_to_nanopore_gt/nano_mamp_pipes/analysis/.snakemake/tmp.43xh034t/126.jobfinished' || (touch '/nfs/jbailey5/baileyweb/gtollefs/mips_to_nanopore_gt/nano_mamp_pipes/analysis/.snakemake/tmp.43xh034t/126.jobfailed'; exit 1)

