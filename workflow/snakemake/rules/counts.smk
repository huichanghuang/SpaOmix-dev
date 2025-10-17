import os

def get_references(species):
    try:
        genomedir = config["references"][species]
    except Exception:
        raise ValueError(f"配置文件中，references 需要配置对应的物种 {species}")
    return genomedir


def _set_count_params(wildcards):

    space_sketcher = "/data03/lead/userdata/huanghuichang/Software/miniconda3/bin/space-sketcher"
    sample = wildcards.sample
    outdir = config["outdir"]
    value = config["samples"][wildcards.sample]
    RNA = config["samples"][wildcards.sample]["RNA"].split(";")
    Oligo = config["samples"][wildcards.sample]["Oligo"].split(";")
    coord = config["samples"][wildcards.sample]["coord"]
    species = value.get("species",config["species"])
    # forcecells = value.get("forcecells",config["forcecells"])
    # minrnaumi = value.get("forcecells",config["minrnaumi"])
    genome_dir = get_references(species)

    ##可选参数
    rnachemistry = value.get("rnachemistry", config["rnachemistry"])
    oligochip = value.get("oligochip", config["oligochip"])
    
    cbwhitelist = value.get("cbwhitelist", config["cbwhitelist"])
    cbwhitelist = '' if not cbwhitelist else f"--cbwhitelist {cbwhitelist}"
    
    sbwhitelist = value.get("sbwhitelist", config["sbwhitelist"])
    sbwhitelist = '' if not sbwhitelist else f"--sbwhitelist {sbwhitelist}"
    
    mapparams = value.get("mapparams", config["mapparams"])
    mapparams = '' if not mapparams else f"--mapparams \'{mapparams}\'"
    
    nobam = value.get("nobam", config["nobam"])
    nobampar = '--nobam ' if nobam else ""
    dev = value.get("dev", config["dev"])
    devpar = '--dev ' if dev else ""
    velo = value.get("dev", config["velo"])
    velopar = '--velo ' if velo else ""
    thread = value.get("thread", config["thread"])
    
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
        f'-t {thread} '
        # f'--minrnaumi {minrnaumi} '
        # f'--forcecells {forcecells} '
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
    threads: config["thread"],
    resources:
        mem_gib=60
    log: os.path.join(config["outdir"], "log/{sample}.analysis.log"),
    shell:
        """
        echo "{params.cmd}" >& {log}
        {params.cmd} >>{log} 2>&1
        """

rule check_contamination:
    input:
        os.path.join(config["outdir"], "{sample}/.all_analysis.done"),
    output:
        touch(os.path.join(config["outdir"], "{sample}/.check_contamination.done")),
    params:
        Rpath = os.path.join(config["conda_env_path"], "bin/Rscript"),
        Rscript = os.path.join(SMK_PATH, "src/check_contamination.R"),
        matrixdir = os.path.join(config["outdir"], "{sample}/04.report/SCST"),
        outdir = os.path.join(config["outdir"], "{sample}"),
    threads: config["thread"],
    log:
        os.path.join(config["outdir"], "log/{sample}.check_contamination.log"),
    shell:
        "time {params.Rpath} {params.Rscript} --matrixdir {params.matrixdir} "
        "--outdir {params.outdir} --samplename {wildcards.sample} >{log} 2>&1\n"

rule cb_distance:
    input:
        os.path.join(config["outdir"], "{sample}/.all_analysis.done"),
    output:
        os.path.join(config["outdir"], "{sample}/{sample}.cb_distance.png"),
    params:
        python = os.path.join(config["conda_env_path"], "bin/python"),
        pysrc = os.path.join(SMK_PATH, "src/calculate_distance_distribution.py"),
        infile = os.path.join(config["outdir"], "{sample}/04.report/SCST/spatial_location_information.txt.gz"),
        outdir = os.path.join(config["outdir"], "{sample}"),
    threads: config["thread"],
    log:
        os.path.join(config["outdir"], "log/{sample}.cb_distance.log"),
    shell:
        "time {params.python} {params.pysrc} --infile {params.infile} --outfile {output} >{log} 2>&1\n"


rule sunburst_plot:
    input:
        os.path.join(config["outdir"], "{sample}/.all_analysis.done"),
    output:
        os.path.join(config["outdir"], "{sample}/{sample}.sunburst_plot.html"),
    params:
        python = os.path.join(config["conda_env_path"], "bin/python"),
        pysrc = os.path.join(SMK_PATH, "src/sunburst_plot.py"),
        cellreads = os.path.join(config["outdir"], "{sample}/01.count/Solo.out/GeneFull_Ex50pAS/CellReads.stats"),
        filtered_barcodes = os.path.join(config["outdir"], "{sample}/01.count/Solo.out/GeneFull_Ex50pAS/filtered/barcodes.tsv"),
        clustered_barcodes = os.path.join(config["outdir"], "{sample}/02.oligo/filtered_matrix/barcodes.tsv.gz"),
        summaryfile = os.path.join(config["outdir"], "{sample}/04.report/summary.csv"),
        sb_umis = os.path.join(config["outdir"], "{sample}/02.oligo/spatial_umis.csv.gz"),
        coord = lambda wildcards: config["samples"][wildcards.sample]["coord"],
        oligochip = lambda wildcards: config["samples"][wildcards.sample].get("oligochip", config["oligochip"]),
    threads: config["thread"],
    log:
        os.path.join(config["outdir"], "log/{sample}.sunburst_plot.log"),
    shell:
        "time {params.python} {params.pysrc} --cellreads {params.cellreads} --filtered_barcodes {params.filtered_barcodes} "
        "--clustered_barcodes {params.clustered_barcodes} --summaryfile {params.summaryfile} --sb_umis {params.sb_umis} "
        "--coordfile {params.coord} --oligochip {params.oligochip} --output {output} >{log} 2>&1\n"


rule plot_chip_dot:
    input:
        os.path.join(config["outdir"], "{sample}/.all_analysis.done"),
    output:
        touch(os.path.join(config["outdir"], "{sample}/{sample}.plot_chip_dot.done")),
    params:
        python = os.path.join(config["conda_env_path"], "bin/python"),
        pysrc = os.path.join(SMK_PATH, "src/plot_dot.py"),
        indir = os.path.join(config["outdir"], "{sample}"),
        coord = lambda wildcards: config["samples"][wildcards.sample]["coord"],
    threads: config["thread"],
    log:
        os.path.join(config["outdir"], "log/{sample}.plot_chip_dot.log"),
    shell:
        "time {params.python} {params.pysrc} {params.coord} {params.indir} >{log} 2>&1\n"
