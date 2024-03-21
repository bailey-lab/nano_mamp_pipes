#!/bin/sh
# properties = {"type": "single", "rule": "index_bam", "local": false, "input": ["/nfs/jbailey5/baileyweb/gtollefs/mips_to_nanopore_gt/240215_MtoN_P2/dorado_pipeline_output/mapping/bam_files/RW-26026-012012-ArtR23RW-1.sorted.bam"], "output": ["/nfs/jbailey5/baileyweb/gtollefs/mips_to_nanopore_gt/240215_MtoN_P2/dorado_pipeline_output/mapping/bam_files/RW-26026-012012-ArtR23RW-1.sorted.bam.bai"], "wildcards": {"sample": "RW-26026-012012-ArtR23RW-1"}, "params": {}, "log": [], "threads": 1, "resources": {"mem_mb": 4000, "mem_mib": 3815, "disk_mb": 1000, "disk_mib": 954, "tmpdir": "<TBD>", "nodes": 1, "time_min": 120}, "jobid": 250, "cluster": {}}
cd '/nfs/jbailey5/baileyweb/gtollefs/mips_to_nanopore_gt/nano_mamp_pipes/analysis' && /oscar/home/gtollefs/mambaforge/envs/snakemake/bin/python3.11 -m snakemake --snakefile '/nfs/jbailey5/baileyweb/gtollefs/mips_to_nanopore_gt/nano_mamp_pipes/analysis/mips_to_nanopore_alignment.smk' --target-jobs 'index_bam:sample=RW-26026-012012-ArtR23RW-1' --allowed-rules 'index_bam' --cores 'all' --attempt 1 --force-use-threads  --resources 'mem_mb=4000' 'mem_mib=3815' 'disk_mb=1000' 'disk_mib=954' 'nodes=1' 'time_min=120' --wait-for-files '/nfs/jbailey5/baileyweb/gtollefs/mips_to_nanopore_gt/nano_mamp_pipes/analysis/.snakemake/tmp.43xh034t' '/nfs/jbailey5/baileyweb/gtollefs/mips_to_nanopore_gt/240215_MtoN_P2/dorado_pipeline_output/mapping/bam_files/RW-26026-012012-ArtR23RW-1.sorted.bam' --force --keep-target-files --keep-remote --max-inventory-time 0 --nocolor --notemp --no-hooks --nolock --ignore-incomplete --rerun-triggers 'input' 'params' 'code' 'mtime' 'software-env' --skip-script-cleanup  --conda-frontend 'mamba' --wrapper-prefix 'https://github.com/snakemake/snakemake-wrappers/raw/' --latency-wait 60 --scheduler 'ilp' --scheduler-solver-path '/oscar/home/gtollefs/mambaforge/envs/snakemake/bin' --default-resources 'mem_mb=4000' 'disk_mb=max(2*input.size_mb, 1000)' 'tmpdir=system_tmpdir' 'nodes=1' 'time_min=120' --mode 2 && touch '/nfs/jbailey5/baileyweb/gtollefs/mips_to_nanopore_gt/nano_mamp_pipes/analysis/.snakemake/tmp.43xh034t/250.jobfinished' || (touch '/nfs/jbailey5/baileyweb/gtollefs/mips_to_nanopore_gt/nano_mamp_pipes/analysis/.snakemake/tmp.43xh034t/250.jobfailed'; exit 1)

