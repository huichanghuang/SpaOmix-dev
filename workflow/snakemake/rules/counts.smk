import os

def get_references(species):
    try:
        genomedir = config["references"][species]
    except Exception:
        raise ValueError(f"配置文件中，references 需要配置对应的物种 {species}")
    return genomedir


def _set_count_params(wildcards):
    """
    
    """
    space_sketcher = "/data03/lead/userdata/huanghuichang/Software/miniconda3/envs/space-sketcker/bin/space-sketcher"
    sample = wildcards.sample
    outdir = config["outdir"]
    value = config["samples"][wildcards.sample]
    RNA = config["samples"][wildcards.sample]["RNA"].split(";")
    Oligo = config["samples"][wildcards.sample]["Oligo"].split(";")
    coord = config["samples"][wildcards.sample]["coord"]
    species = value.get("species",config["species"])
    genome_dir = get_references(species)

    ##可选参数
    rnachemistry = value.get("rnachemistry", config["rnachemistry"])
    oligochip = value.get("oligochip", config["oligochip"])
    
    cbwhitelist = value.get("cbwhitelist", config["cbwhitelist"])
    cbwhitelist = '' if not cbwhitelist else f"--cbwhitelist {cbwhitelist}"
    
    sbwhitelist = value.get("sbwhitelist", config["sbwhitelist"])
    sbwhitelist = '' if not sbwhitelist else f"--sbwhitelist {sbwhitelist}"
    
    mapparams = value.get("mapparams", config["mapparams"])
    mapparams = '' if not mapparams else f"--mapparams {mapparams}"
    
    nobam = value.get("nobam", config["nobam"])
    nobampar = '--nobam ' if nobam else ""
    dev = value.get("dev", config["dev"])
    devpar = '--dev ' if dev else ""
    velo = value.get("dev", config["velo"])
    velopar = '--velo ' if velo else ""
    threads = value.get("threads", config["threads"])
    
    cmd = (
        f'{space_sketcher} rna run '
        f'--genomeDir {genome_dir} '
        f'--rna1 {RNA[0]} --rna2 {RNA[1]} '
        f'--oligor1 {Oligo[0]} --oligor2 {Oligo[1]} '
        f'--rnachemistry {rnachemistry} '
        f'--oligochip {oligochip} '
        f'--coordfile {coord} '
        f'{cbwhitelist} '
        f'{sbwhitelist} '
        f'{mapparams} '
        f'--outdir {outdir} --name {sample} '
        f'--reference {species} '
        f'-t {threads} '
        f'{devpar} {nobampar} {velopar}'
    )
    return cmd


rule counts:
    input:
        RNA=lambda wildcards: config["samples"][wildcards.sample]["RNA"].split(";"),
        Oligo=lambda wildcards: config["samples"][wildcards.sample]["Oligo"].split(";"),
    output:
        touch(os.path.join(config["outdir"], "{sample}/.all_analysis.done")),
    params:
        cmd = lambda wildcards: _set_count_params(wildcards),
    threads: config["threads"],
    resources:
        mem_gib=60
    log: os.path.join(config["outdir"], "log/{sample}.analysis.log"),
    shell:
        """
        echo "{params.cmd}" >& {log}
        {params.cmd} >>{log} 2>&1
        """
