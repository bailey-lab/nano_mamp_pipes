#!/bin/sh
# properties = {"type": "single", "rule": "make_sam", "local": false, "input": ["/nfs/jbailey5/baileyweb/bailey_share/resources/MIP_species_resources/pf/Pf_3D7/genomes/genome.fa", "/nfs/jbailey5/baileyweb/asimkin/mips_to_nanopore/240215_MtoN_P2/prepped_reads/demuxed_fastq/fastq/RW013013-263263-ArtR23RW-1.fastq.gz"], "output": ["/nfs/jbailey5/baileyweb/gtollefs/mips_to_nanopore_gt/240215_MtoN_P2/output_rerun_for_completeness_check/mapping/sam_files/RW013013-263263-ArtR23RW-1.sam"], "wildcards": {"sample": "RW013013-263263-ArtR23RW-1"}, "params": {}, "log": [], "threads": 1, "resources": {"mem_mb": 4000, "mem_mib": 3815, "disk_mb": 1000, "disk_mib": 954, "tmpdir": "<TBD>", "nodes": 1, "time_min": 120}, "jobid": 204, "cluster": {}}
cd '/nfs/jbailey5/baileyweb/gtollefs/mips_to_nanopore_gt/nano_mamp_pipes/analysis' && /oscar/home/gtollefs/mambaforge/envs/snakemake/bin/python3.11 -m snakemake --snakefile '/nfs/jbailey5/baileyweb/gtollefs/mips_to_nanopore_gt/nano_mamp_pipes/analysis/mips_to_nanopore.smk' --target-jobs 'make_sam:sample=RW013013-263263-ArtR23RW-1' --allowed-rules 'make_sam' --cores 'all' --attempt 1 --force-use-threads  --resources 'mem_mb=4000' 'mem_mib=3815' 'disk_mb=1000' 'disk_mib=954' 'nodes=1' 'time_min=120' --wait-for-files '/nfs/jbailey5/baileyweb/gtollefs/mips_to_nanopore_gt/nano_mamp_pipes/analysis/.snakemake/tmp.sqysjyfu' '/nfs/jbailey5/baileyweb/bailey_share/resources/MIP_species_resources/pf/Pf_3D7/genomes/genome.fa' '/nfs/jbailey5/baileyweb/asimkin/mips_to_nanopore/240215_MtoN_P2/prepped_reads/demuxed_fastq/fastq/RW013013-263263-ArtR23RW-1.fastq.gz' --force --keep-target-files --keep-remote --max-inventory-time 0 --nocolor --notemp --no-hooks --nolock --ignore-incomplete --rerun-triggers 'input' 'mtime' 'code' 'software-env' 'params' --skip-script-cleanup  --conda-frontend 'mamba' --wrapper-prefix 'https://github.com/snakemake/snakemake-wrappers/raw/' --latency-wait 60 --scheduler 'ilp' --scheduler-solver-path '/oscar/home/gtollefs/mambaforge/envs/snakemake/bin' --default-resources 'mem_mb=4000' 'disk_mb=max(2*input.size_mb, 1000)' 'tmpdir=system_tmpdir' 'nodes=1' 'time_min=120' --mode 2 && touch '/nfs/jbailey5/baileyweb/gtollefs/mips_to_nanopore_gt/nano_mamp_pipes/analysis/.snakemake/tmp.sqysjyfu/204.jobfinished' || (touch '/nfs/jbailey5/baileyweb/gtollefs/mips_to_nanopore_gt/nano_mamp_pipes/analysis/.snakemake/tmp.sqysjyfu/204.jobfailed'; exit 1)

