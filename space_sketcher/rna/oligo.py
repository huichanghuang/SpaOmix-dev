import os
from pathlib import Path
import typer
from typing_extensions import Annotated
from typing import Optional
from space_sketcher.__init__ import __root_dir__

class Oligo:
    def __init__(self, args):
        self.name = args.name
        self.outdir = os.path.abspath(os.path.join(args.outdir, args.name))
        self.oligor1 = args.oligor1
        self.oligor2 = args.oligor2
        self.rnachemistry = args.rnachemistry
        self.mapparams = args.mapparams
        self.oligochip = args.oligochip
        self.threads = args.threads
        self.coordfile = args.coordfile
        self.maxumi = args.maxumi
        self.minumi = args.minumi
        self.eps = args.eps
        self.min_samples = args.min_samples

        if args.cbwhitelist is not None:
            self.cbwhitelist = args.cbwhitelist
        elif self.rnachemistry != "other" and (self.mapparams is None or "--soloCBwhitelist" not in self.mapparams):
            self.cbwhitelist = os.path.join(__root_dir__, 
                              f"data/cbwhitelist/{self.rnachemistry}/{self.rnachemistry}.cbwhitelist.txt")   
        else:
            self.cbwhitelist = None

        self.sbwhitelist = (
            args.sbwhitelist
            if args.sbwhitelist is not None
            else os.path.join(__root_dir__, 
                              "data/sbwhitelist/sbwhitelist.txt")
        )               

    def run(self):
        ### import lib
        from space_sketcher.tools.utils import str_mkdir, judgeFilexits, get_formatted_time
        from space_sketcher.__init__ import __root_dir__
        from space_sketcher.rna.src.spatial_barcode_extraction import stat_spatial_barcodes
        from space_sketcher.rna.src.assign_coordinate import assign_coordinate
        from space_sketcher.rna.src.dbscan_filter import dbscan_filter
        print("test oligo")
        
        oligodir = os.path.join(self.outdir, "02.oligo")
        str_mkdir(oligodir)

        ###Extract spatial barcodes
        # judge file exits
        judgeFilexits(
            self.oligor1,
            self.oligor2,
            self.cbwhitelist,
            self.sbwhitelist
            )
        # run stat_spatial_barcodes
        if not os.path.exists(f"{oligodir}/spatial_umis.csv.gz"):
            print(f'\n{get_formatted_time()}\n'
                f'Extracting spatial barcode information.')
            stat_spatial_barcodes(self.oligor1, self.oligor2, 
                                self.oligochip, 
                                self.rnachemistry, self.mapparams,
                                self.cbwhitelist, self.sbwhitelist,
                                oligodir, self.threads)
        else:
            print(f"{oligodir}/spatial_umis.csv.gz already exits, skip Extracting spatial barcode information.")
            
        ###Assign spatial barcodes coordinate
        # judge file exits
        judgeFilexits(
            self.coordfile,
            f"{self.outdir}/01.count/Solo.out/GeneFull_Ex50pAS/filtered/barcodes.tsv",
            f"{oligodir}/spatial_umis.csv.gz",
            f"{oligodir}/sb_library_summary.temp.csv",
            )

        #run assign_coordinate
        print(f'\n{get_formatted_time()}\n'
            f'Assigning coordinate for each spatial barcode.')
        assign_coordinate(self.coordfile, self.oligochip, f"{oligodir}/spatial_umis.csv.gz",
                          f"{self.outdir}/01.count/Solo.out/GeneFull_Ex50pAS/filtered/barcodes.tsv",
                          f"{oligodir}/sb_library_summary.temp.csv",
                          self.sbwhitelist, oligodir)
        
        ###Filter cell barcode by dbscan clustering
        # judge file exits
        judgeFilexits(
            f"{self.outdir}/01.count/Solo.out/GeneFull_Ex50pAS/filtered",
            f"{self.outdir}/01.count/Solo.out/GeneFull_Ex50pAS/CellReads.stats",
            f"{oligodir}/cb_sb_coord.txt",
            )

        #run dbscan_filter
        print(f'\n{get_formatted_time()}\n'
            f'Performing dbscan filtering.')
        dbscan_filter(f"{oligodir}/cb_sb_coord.txt",
                      oligodir, 
                      self.maxumi, self.minumi, 
                      f"{self.outdir}/01.count/Solo.out/GeneFull_Ex50pAS/filtered", 
                      f"{self.outdir}/01.count/Solo.out/GeneFull_Ex50pAS/CellReads.stats", 
                      self.eps, self.min_samples, self.threads)        

        (Path(self.outdir) / ".oligo.done").touch()


def oligo_app(
    # 必需参数
    name: Annotated[str, typer.Option("--name", "-n", help="Sample name", prompt=True, show_default=False)],
    oligor1: Annotated[str, typer.Option("--oligor1", "-o1", help="Oligo R1 fastq file", prompt=True, show_default=False)],
    oligor2: Annotated[str, typer.Option("--oligor2", "-o2", help="Oligo R2 fastq file", prompt=True, show_default=False)],
    coordfile: Annotated[str, typer.Option("--coordfile", "-cf", help="Coordinate file with spatial barcodes", prompt=True, show_default=False)],
    # 可选参数
    outdir: Annotated[str, typer.Option("--outdir", "-o", help="Output directory ")] = os.getcwd(),
    cbwhitelist: Annotated[Optional[str], typer.Option("--cbwhitelist", "-cb", help="Cell barcode whitelist file")] = None,
    sbwhitelist: Annotated[Optional[str], typer.Option("--sbwhitelist", "-sb", help="Spatial barcode whitelist file")] = None,
    rnachemistry: Annotated[str, typer.Option("--rnachemistry", "-rc", help="Chemistry version: 10X/leader_v1/other")] = "leader_v1",
    mapparams: Annotated[Optional[str], typer.Option("--mapparams", "-mp", help="Additional STAR mapping parameters (required for 'other' chemistry)")] = None,
    oligochip: Annotated[str, typer.Option("--oligochip", "-oc", help="Spatial chip version: LD/GM")] = "LD",
    threads: Annotated[int, typer.Option("--threads", "-t", help="Number of threads")] = 4,
    maxumi: Annotated[int, typer.Option("--maxumi", "-max", help="Max UMI threshold for filtering")] = 5000,
    minumi: Annotated[int, typer.Option("--minumi", "-min", help="Min UMI threshold for filtering")] = 2,
    eps: Annotated[float, typer.Option("--eps", "-e", help="DBSCAN epsilon (max distance)")] = 150.0,
    min_samples: Annotated[int, typer.Option("--min_samples", "-ms", help="DBSCAN min points per cluster")] = 6
):
    """
    Process spatial oligo sequencing data.\n
    
    Example: \n
        space-sketcher rna run --name S1 --oligor1 O1.fq --oligor2 O2.fq --coordfile coords.txt --outdir /path/to/outdir
    """
    # 将参数转换为类似argparse的Namespace对象
    class Args:
        pass
    args = Args()
    for k, v in locals().items():
        if k != 'args':
            setattr(args, k, v)
    
    # 执行处理流程
    processor = Oligo(args)
    processor.run()

# 导出函数
__all__ = ["oligo_app"]
