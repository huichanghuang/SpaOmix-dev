import os
import collections
import time
import typer
from typing_extensions import Annotated
from typing import Optional
from space_sketcher.__init__ import __root_dir__


class Runpipe:
    def __init__(self, args):
        self.name = args.name
        self.outdir = os.path.abspath(args.outdir)
        self.rna1 = args.rna1
        self.rna2 = args.rna2
        self.oligor1 = args.oligor1
        self.oligor2 = args.oligor2
        self.genomeDir = os.path.abspath(args.genomeDir)
        self.coordfile = args.coordfile
        self.rnachemistry = args.rnachemistry
        self.oligochip = args.oligochip        
        self.mapparams = args.mapparams
        self.threads = args.threads        
        self.calling_method = args.calling_method
        self.maxoligoumi = args.maxoligoumi
        self.minoligoumi = args.minoligoumi
        self.eps = args.eps
        self.min_samples = args.min_samples
        self.minfeatures = args.minfeatures
        self.mincells = args.mincells
        self.ndims = args.ndims
        self.nvariables = args.nvariables
        self.resolution = args.resolution
        self.forcecells = args.forcecells
        self.minrnaumi = args.minrnaumi
        self.dev = args.dev
        self.nobam = args.nobam
        self.velo = args.velo

        if args.cbwhitelist is not None:
            self.cbwhitelist = args.cbwhitelist
        elif self.rnachemistry != "other" and (self.mapparams is None or "--soloCBwhitelist" not in self.mapparams):
            self.cbwhitelist = os.path.join(__root_dir__, 
                              f"data/cbwhitelist/{self.rnachemistry}/{self.rnachemistry}.cbwhitelist.txt")   
        else:
            self.cbwhitelist = None

        if args.sbwhitelist is not None:
            self.sbwhitelist = args.sbwhitelist
        elif self.oligochip == "LD":
            self.sbwhitelist = os.path.join(__root_dir__, 
                              "data/sbwhitelist/sbwhitelist.txt")   
        else:
            self.sbwhitelist = None

        self.reference = (
            args.reference
            if args.reference is not None
            else os.path.basename(os.path.basename(self.genomeDir))
        )

    def runpipe(self):

        ### import lib
        from space_sketcher.tools.utils import (
            str_mkdir, 
            judgeFilexits,
            execute_and_log,
            bin_path,
            rm_temp,
        )
        
        print("test run")
        ### run       
        judgeFilexits(
            self.rna1,
            self.rna2,
            self.oligor1,
            self.oligor2,
            self.genomeDir,
            self.coordfile,
            )

        print(bin_path())

        if self.cbwhitelist is None and self.mapparams is None:
            print("Both cbwhitelist and mapparams are None, please check!")
            exit(1)
        # Base command components
        count_cmd = [
            f"{bin_path()}/space-sketcher rna count",
            f"--name {self.name}",
            f"--outdir {self.outdir}",
            f"--rna1 {self.rna1}",
            f"--rna2 {self.rna2}",
            f"--rnachemistry {self.rnachemistry}",
            f"--threads {self.threads}",
            f"--genomeDir {self.genomeDir}",
            f"--calling_method {self.calling_method}",
            f"--forcecells {self.forcecells}",
            f"--minumi {self.minrnaumi}",
        ]
        # Add optional parameters if they exist
        if self.cbwhitelist is not None:
            count_cmd.append(f"--cbwhitelist {self.cbwhitelist}")
        if self.mapparams is not None:
            count_cmd.append(f"--mapparams \'{self.mapparams}\'")
        if self.velo:
            count_cmd.append(f"--velo")
        count_cmd  = ' '.join(count_cmd) 


        oligo_cmd = [
            f"{bin_path()}/space-sketcher rna oligo",
            f"--name {self.name}",
            f"--outdir {self.outdir}",
            f"--oligor1 {self.oligor1}",
            f"--oligor2 {self.oligor2}",
            f"--rnachemistry {self.rnachemistry}",
            f"--oligochip {self.oligochip}",
            f"--threads {self.threads}",
            f"--coordfile {self.coordfile}",
            f"--maxumi {self.maxoligoumi}",
            f"--minumi {self.minoligoumi}",
            f"--eps {self.eps}",
            f"--min_samples {self.min_samples}",
        ]
        # Add optional parameters if they exist
        if self.cbwhitelist is not None:
            oligo_cmd.append(f"--cbwhitelist {self.cbwhitelist}")
        if self.sbwhitelist is not None:
            oligo_cmd.append(f"--sbwhitelist {self.sbwhitelist}")
        if self.mapparams is not None:
            oligo_cmd.append(f"--mapparams \'{self.mapparams}\'")
        oligo_cmd  = ' '.join(oligo_cmd) 


        analysis_cmd = [
            f'{bin_path()}/space-sketcher rna analysis',
            f"--name {self.name}",
            f"--outdir {self.outdir}",
            f"--minfeatures {self.minfeatures}",
            f"--mincells {self.mincells}",
            f"--ndims {self.ndims}",
            f"--nvariables {self.nvariables}",
            f"--resolution {self.resolution}",
        ]
        analysis_cmd  = ' '.join(analysis_cmd)
        
        
        report_cmd = [
            f'{bin_path()}/space-sketcher rna report',
            f'--name {self.name}',
            f"--outdir {self.outdir}",
            f"--reference {self.reference}",
            f"--rnachemistry {self.rnachemistry}",
            f"--oligochip {self.oligochip}",
        ]
        if self.dev:
            report_cmd.append(f"--dev")
        report_cmd = ' '.join(report_cmd)
        
        cmdlist = collections.OrderedDict()
        cmdlist['count'] = count_cmd
        cmdlist['oligo'] = oligo_cmd
        cmdlist['analysis'] = analysis_cmd
        cmdlist['report'] = report_cmd

        logdir = os.path.join(self.outdir,self.name)
        str_mkdir(logdir)
        start_time = time.time()
        for pipe, pipecmd in cmdlist.items():
            execute_and_log(pipecmd, pipe, logdir)
        
        end_time = time.time()
        analysis_time = end_time - start_time
        analysis_time_minutes, analysis_time_seconds = divmod(analysis_time, 60)
        analysis_time_hours, analysis_time_minutes = divmod(analysis_time_minutes, 60)

        print(f'\nAnalysis Finished')
        print(f'Elapsed Time: {int(analysis_time_hours)} hours {int(analysis_time_minutes)} minutes {int(analysis_time_seconds)} seconds')

        ###remove bamfile if no bam
        bamfile = os.path.join(self.outdir,self.name,"01.count/Aligned.sortedByCoord.out.bam")
        if self.nobam and os.path.exists(bamfile):
            rm_temp(bamfile)

        ###rm temp file
        if not self.dev:
            rm_temp(f"{self.outdir}/{self.name}/*/*.temp.*")
            rm_temp(f"{self.outdir}/{self.name}/02.oligo/*.png")
            rm_temp(f"{self.outdir}/{self.name}/02.oligo/dbscan*")
            rm_temp(f"{self.outdir}/{self.name}/02.oligo/CB_cluster.txt")



def run_app(
    # 必需参数
    name: Annotated[str, typer.Option("--name", "-n", help="Sample name", prompt=True, show_default=False)],
    rna1: Annotated[str, typer.Option("--rna1", "-r1", help="RNA R1 fastq file", prompt=True, show_default=False)],
    rna2: Annotated[str, typer.Option("--rna2", "-r2", help="RNA R2 fastq file", prompt=True, show_default=False)],
    oligor1: Annotated[str, typer.Option("--oligor1", "-o1", help="Oligo R1 fastq file", prompt=True, show_default=False)],
    oligor2: Annotated[str, typer.Option("--oligor2", "-o2", help="Oligo R2 fastq file", prompt=True, show_default=False)],
    genomeDir: Annotated[str, typer.Option("--genomeDir", "-g", help="Path to genome directory", prompt=True, show_default=False)],
    coordfile: Annotated[str, typer.Option("--coordfile", "-cf", help="Coordinate file with spatial barcodes", prompt=True, show_default=False)],
    
    # 可选参数
    outdir: Annotated[str, typer.Option("--outdir", "-o", help=f"Output directory")] = os.getcwd(),
    cbwhitelist: Annotated[Optional[str], typer.Option("--cbwhitelist", "-cb", help="Cell barcode whitelist file")] = None,
    sbwhitelist: Annotated[Optional[str], typer.Option("--sbwhitelist", "-sb", help="Spatial barcode whitelist file")] = None,
    rnachemistry: Annotated[str, typer.Option("--rnachemistry", "-rc", help="Chemistry version: 10X/leader_v1/other")] = "leader_v1",
    oligochip: Annotated[str, typer.Option("--oligochip", "-oc", help="Spatial chip version: LD/GM")] = "LD",
    mapparams: Annotated[Optional[str], typer.Option("--mapparams", "-mp", help="Additional STAR mapping parameters (required for 'other' chemistry)")] = None,
    threads: Annotated[int, typer.Option("--threads", "-t", help="Number of threads")] = 4,
    calling_method: Annotated[str, typer.Option("--calling_method", help="Cell calling method: CellRanger2.2/EmptyDrops_CR")] = "EmptyDrops_CR",
    maxoligoumi: Annotated[int, typer.Option("--maxumi", help="Max oligo UMI threshold for filtering")] = 5000,
    minoligoumi: Annotated[int, typer.Option("--minoligoumi", help="Min oligo UMI threshold for filtering")] = 2,
    eps: Annotated[float, typer.Option("--eps", help="DBSCAN epsilon (max distance)")] = 150.0,
    min_samples: Annotated[int, typer.Option("--min_samples", help="DBSCAN min points per cluster")] = 6,
    minfeatures: Annotated[int, typer.Option("--minfeatures", help="Min features per cell")] = 5,
    mincells: Annotated[int, typer.Option("--mincells", help="Min cells per gene")] = 3,
    ndims: Annotated[int, typer.Option("--ndims", help="PCA dimensions")] = 30,
    nvariables: Annotated[int, typer.Option("--nvariables", help="Number of variable genes")] = 2000,
    resolution: Annotated[float, typer.Option("--resolution", help="Clustering resolution")] = 0.5,
    reference: Annotated[Optional[str], typer.Option("--reference", help="Reference name (default: genomeDir basename)", show_default=False)] = None,
    forcecells: Annotated[int, typer.Option("--forcecells", help="Force pipeline to use this number of beads, bypassing cell calling algorithm.")] = 0,
    minrnaumi: Annotated[int, typer.Option("--minrnaumi", help="The min rna umi for use emptydrops")] = 200,
    dev: Annotated[bool, typer.Option("--dev", help="Enable development mode (default: False)")] = False,
    nobam: Annotated[bool, typer.Option("--nobam", help="Remove BAM files after running (default: False)")] = False,
    velo: Annotated[bool, typer.Option("--velo", help="Enable STARsolo Velocyto mode (default: False)")] = False
):
    """
    Run spatial RNA analysis pipeline.\n
    
    Example: \n
        space-sketcher rna run --name S1 --rna1 R1.fq --rna2 R2.fq --oligor1 O1.fq --oligor2 O2.fq --genomeDir /path/to/genome --coordfile coords.txt --outdir /path/to/outdir
    """
    # 将参数打包成类似argparse的Namespace对象
    class Args:
        pass
    args = Args()
    
    # 动态设置属性
    for k, v in locals().items():
        if k != 'args':
            setattr(args, k, v)
    
    # 如果没有指定reference，使用genomeDir的basename
    if args.reference is None:
        args.reference = os.path.basename(os.path.normpath(args.genomeDir))
    
    # 执行处理流程
    processor = Runpipe(args)
    processor.runpipe()

# 导出函数
__all__ = ["run_app"]
