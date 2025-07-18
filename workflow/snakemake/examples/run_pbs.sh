#! /usr/bin/bash
snakemake --snakefile /data03/lead/userdata/huanghuichang/Work/Pipeline/Space-sketcher-dev/workflow/snakemake/Snakefile \
    --configfile input.yaml \
    --executor slurm \
    --jobs 4 -n
    # --cluster-generic-submit-cmd "qsub -V -l walltime=1000:00:00" \
    
