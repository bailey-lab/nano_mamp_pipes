#!/bin/sh
# properties = {"type": "single", "rule": "sort_bam", "local": false, "input": ["/nfs/jbailey5/baileyweb/gtollefs/mips_to_nanopore_gt/240215_MtoN_P2/dorado_pipeline_output/mapping/bam_files/RW-17017-455455-ArtR23RW-1.bam"], "output": ["/nfs/jbailey5/baileyweb/gtollefs/mips_to_nanopore_gt/240215_MtoN_P2/dorado_pipeline_output/mapping/bam_files/RW-17017-455455-ArtR23RW-1.sorted.bam"], "wildcards": {"sample": "RW-17017-455455-ArtR23RW-1"}, "params": {}, "log": [], "threads": 1, "resources": {"mem_mb": 4000, "mem_mib": 3815, "disk_mb": 1000, "disk_mib": 954, "tmpdir": "<TBD>", "nodes": 1, "time_min": 120}, "jobid": 55, "cluster": {}}
cd '/nfs/jbailey5/baileyweb/gtollefs/mips_to_nanopore_gt/nano_mamp_pipes/analysis' && /oscar/home/gtollefs/mambaforge/envs/snakemake/bin/python3.11 -m snakemake --snakefile '/nfs/jbailey5/baileyweb/gtollefs/mips_to_nanopore_gt/nano_mamp_pipes/analysis/mips_to_nanopore_alignment.smk' --target-jobs 'sort_bam:sample=RW-17017-455455-ArtR23RW-1' --allowed-rules 'sort_bam' --cores 'all' --attempt 1 --force-use-threads  --resources 'mem_mb=4000' 'mem_mib=3815' 'disk_mb=1000' 'disk_mib=954' 'nodes=1' 'time_min=120' --wait-for-files '/nfs/jbailey5/baileyweb/gtollefs/mips_to_nanopore_gt/nano_mamp_pipes/analysis/.snakemake/tmp.43xh034t' '/nfs/jbailey5/baileyweb/gtollefs/mips_to_nanopore_gt/240215_MtoN_P2/dorado_pipeline_output/mapping/bam_files/RW-17017-455455-ArtR23RW-1.bam' --force --keep-target-files --keep-remote --max-inventory-time 0 --nocolor --notemp --no-hooks --nolock --ignore-incomplete --rerun-triggers 'input' 'params' 'code' 'mtime' 'software-env' --skip-script-cleanup  --conda-frontend 'mamba' --wrapper-prefix 'https://github.com/snakemake/snakemake-wrappers/raw/' --latency-wait 60 --scheduler 'ilp' --scheduler-solver-path '/oscar/home/gtollefs/mambaforge/envs/snakemake/bin' --default-resources 'mem_mb=4000' 'disk_mb=max(2*input.size_mb, 1000)' 'tmpdir=system_tmpdir' 'nodes=1' 'time_min=120' --mode 2 && touch '/nfs/jbailey5/baileyweb/gtollefs/mips_to_nanopore_gt/nano_mamp_pipes/analysis/.snakemake/tmp.43xh034t/55.jobfinished' || (touch '/nfs/jbailey5/baileyweb/gtollefs/mips_to_nanopore_gt/nano_mamp_pipes/analysis/.snakemake/tmp.43xh034t/55.jobfailed'; exit 1)

