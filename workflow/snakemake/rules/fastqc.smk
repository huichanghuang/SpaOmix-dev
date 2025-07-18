import os

rule fastq_qc:
    input: 
        RNA=lambda wildcards: config["samples"][wildcards.sample]["RNA"].split(";"),
        Oligo=lambda wildcards: config["samples"][wildcards.sample]["Oligo"].split(";")
    output:
        touch(os.path.join(config["outdir"], "{sample}/.fastqc.done")),
    params:
        fastqc = "/data03/lead/userdata/huanghuichang/Software/miniconda3/bin/fastqc",
        out_dir = os.path.join(config["outdir"], "{sample}"),
    threads: 3,
    resources:
        mem_gib=3,
    log: os.path.join(config["outdir"], "log/{sample}.fastq_qc.log"),
    shell:
        """
        {params.fastqc} --outdir {params.out_dir} --threads {threads} {input.RNA} {input.Oligo} >& {log}
        """
