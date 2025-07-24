import os
import subprocess
from pathlib import Path
import typer
from typing_extensions import Annotated
from typing import Optional

class Count:
    def __init__(self,args):
        self.rna1 = args.rna1
        self.rna2 = args.rna2
        self.threads = args.threads
        self.name = args.name
        self.rnachemistry = args.rnachemistry
        self.cbwhitelist = args.cbwhitelist
        self.mapparams = args.mapparams
        self.genomeDir = args.genomeDir
        self.calling_method = args.calling_method
        self.velo = args.velo
        self.forcecells = args.forcecells
        self.minumi = args.minumi
        self.outdir = os.path.abspath(os.path.join(args.outdir,args.name))

    def Prepare_mapping_params(self) -> str:
        """
        Check the rnachemistry and prepare mapping parameters.
        Returns:
            str: The mapping parameters.
        """
        ###load function
        from space_sketcher.tools.utils import gunzip
        ###to avoid STAR memery error
        if self.cbwhitelist is not None and self.cbwhitelist.endswith(".gz"):
            self.cbwhitelist = gunzip(self.cbwhitelist)

        mapping_params = ""
        if self.velo:
            mapping_params += "--soloFeatures GeneFull_Ex50pAS Velocyto "
        else:
            mapping_params += "--soloFeatures GeneFull_Ex50pAS "
                    
        if self.rnachemistry == "10X":
            mapping_params += "--soloType CB_UMI_Simple "
            mapping_params += "--soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 12 "
            mapping_params += "--soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --soloBarcodeReadLength 0 "
            mapping_params += f"--soloCBwhitelist {self.cbwhitelist} "
            if self.mapparams is not None:
                mapping_params += self.mapparams
        elif self.rnachemistry == "leader_v1":
            mapping_params += "--soloType CB_UMI_Complex "
            mapping_params += "--soloCBposition 0_0_0_9 0_10_0_19 --soloUMIposition 0_20_0_29 "
            mapping_params += "--soloCBmatchWLtype EditDist_2 "
            mapping_params += f"--soloCBwhitelist {self.cbwhitelist} {self.cbwhitelist} "
            if self.mapparams is not None:
                mapping_params += self.mapparams
        elif self.rnachemistry == "other" and self.mapparams is not None:
            mapping_params += self.mapparams ###whitelist 在mapparams中提供
        else:
            print("Not available rnachemistry or mapping params")

        return mapping_params
        
    def run(self):
        ### import lib
        from space_sketcher.tools.utils import str_mkdir, judgeFilexits
        from space_sketcher.__init__ import __root_dir__
        from space_sketcher.rna.src.saturation import count_saturation
        ### run
        judgeFilexits(
            self.rna1,
            self.rna2,
            self.genomeDir,
            )

        rnadir = os.path.join(self.outdir, "01.count")
        str_mkdir(rnadir)
        
        star_version = subprocess.check_output(f"{__root_dir__}/software/STAR --version", shell=True)
        print(f"STAR 版本号：{star_version.decode('utf8')}")

        mapping_pars = self.Prepare_mapping_params()
        STAR_cmd = (
            f"{__root_dir__}/software/STAR "
            f"--runMode alignReads "
            "--quantMode GeneCounts "
            f"--soloCellFilter {self.calling_method} "
            "--outFilterScoreMin 30 "
            "--soloStrand Unstranded "
            "--readFilesCommand zcat "
            "--soloCellReadStats Standard "
            "--soloMultiMappers EM "
            "--soloUMIdedup 1MM_CR "
            "--soloUMIfiltering MultiGeneUMI_CR "
            "--clipAdapterType CellRanger4 "
            "--outSAMtype BAM SortedByCoordinate "
            "--outSAMattributes CR CY UR UY NH HI nM AS GX GN gx gn CB UB sS sQ sM "
            "--soloOutFileNames Solo.out/ genes.tsv barcodes.tsv matrix.mtx "
            f"{mapping_pars} "
            f"--readFilesIn {self.rna2} {self.rna1} "
            f"--genomeDir {self.genomeDir} "
            f"--outFileNamePrefix {rnadir}/ "
            f"--runThreadN {self.threads} "
            "--outBAMsortingThreadN 6 "
        )

        if not os.path.exists(Path(self.outdir) / ".STAR.done"):
            subprocess.check_call(STAR_cmd, shell=True)
            ###change output directory permissions
            chmod_cmd = f"chmod -R a+r {rnadir}/Solo.out && chmod a+x $(find {rnadir}/Solo.out -type d)"
            subprocess.check_call(chmod_cmd, shell=True)
            (Path(self.outdir) / ".STAR.done").touch()
        else:
            print(".STAR.done exits, skip running STAR.")

        ###recall cell by umi=200 or forcell
        ###Extract spatial barcodes
        # judge file exits
        matrixdir = os.path.join(rnadir, "Solo.out/GeneFull_Ex50pAS/raw")
        outmatrix = os.path.join(rnadir, "Solo.out/GeneFull_Ex50pAS/callcell")
        judgeFilexits(matrixdir)
        callcell_cmd = (
            f"{__root_dir__}/software/Rscript {__root_dir__}/rna/src/cell_identify.R "
            f"--inputdir {matrixdir} "
            f"--outputdir {outmatrix} "
            f"--minumi {self.minumi} "
            f"--force_cells {self.forcecells} "
        )
        print('Performing cell calling.')
        subprocess.check_call(callcell_cmd, shell=True)

        ###calculate saturation
        count_saturation(rnadir, self.threads)


def count_app(
    # 必需参数
    name: Annotated[str, typer.Option("--name", "-n", help="Sample name", prompt=True, show_default=False)],
    rna1: Annotated[str, typer.Option("--rna1", "-r1", help="RNA R1 fastq file", prompt=True, show_default=False)],
    rna2: Annotated[str, typer.Option("--rna2", "-r2", help="RNA R2 fastq file", prompt=True, show_default=False)],
    genomeDir: Annotated[str, typer.Option("--genomeDir", "-g", help="Path to genome directory", prompt=True, show_default=False)],
    # 可选参数
    outdir: Annotated[str, typer.Option("--outdir", "-o", help=f"Output directory")] = os.getcwd(),
    rnachemistry: Annotated[str, typer.Option("--rnachemistry", "-rc", help="Chemistry version: 10X/leader_v1/other")] = "leader_v1",
    cbwhitelist: Annotated[Optional[str], typer.Option("--cbwhitelist", "-cb", help="Cell barcode whitelist file")] = None,
    mapparams: Annotated[Optional[str], typer.Option("--mapparams", "-mp", help="Additional STAR mapping parameters (required for 'other' chemistry)")] = None,
    threads: Annotated[int, typer.Option("--threads", "-t", help="Number of threads")] = 4,
    calling_method: Annotated[str, typer.Option("--calling_method", help="Cell calling method: CellRanger2.2/EmptyDrops_CR")] = "EmptyDrops_CR",
    forcecells: Annotated[int, typer.Option("--forcecells", help="Force pipeline to use this number of beads, bypassing cell calling algorithm.")] = 0,
    minumi: Annotated[int, typer.Option("--minumi", help="The min umi for use emptydrops")] = 200,
    velo: Annotated[bool, typer.Option("--velo", help="Enable STARsolo Velocyto mode (default: False)")] = False
):
    """
    Run RNA mapping and counting pipeline.\n
    
    Example:\n
        space-sketcher rna count --name S1 --rna1 R1.fq --rna2 R2.fq --genomeDir /path/to/genome --outdir /path/to/outdir
    """
    # 将参数打包成类似argparse的Namespace对象
    class Args:
        pass
    args = Args()
    
    # 动态设置属性
    for k, v in locals().items():
        if k != 'args':
            setattr(args, k, v)
    
    # 执行处理流程
    processor = Count(args)
    processor.run()

# 导出函数
__all__ = ["count_app"]
