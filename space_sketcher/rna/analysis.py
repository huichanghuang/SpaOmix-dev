import os
from pathlib import Path
from typing_extensions import Annotated
import typer

class Analysis:
    def __init__(self,args):
        self.name = args.name
        self.outdir = os.path.abspath(os.path.join(args.outdir,args.name))
        self.minfeatures = args.minfeatures
        self.mincells = args.mincells
        self.ndims = args.ndims
        self.nvariables = args.nvariables
        self.resolution = args.resolution

    
    def run(self):
        ### import lib
        from space_sketcher.tools.utils import str_mkdir, judgeFilexits
        from space_sketcher.rna.src.matrix_QC import perform_matrix_qc
        from space_sketcher.__init__ import __root_dir__
        # import subprocess
        
        anadir = os.path.join(self.outdir, "03.analysis")
        str_mkdir(anadir)

        ###Extract spatial barcodes
        # judge file exits
        matrixdir = os.path.join(self.outdir, "02.oligo/filtered_matrix")
        judgeFilexits(matrixdir)
        # # run perform_matrix_qc by R seurat
        # matrix_cmd = (
        #     f"{__root_dir__}/software/Rscript {__root_dir__}/rna/src/matrix_QC.R "
        #     f"--matrixdir {matrixdir} "
        #     f"--outdir {anadir} "
        #     f"--mincells {self.mincells} "
        #     f"--minfeatures {self.minfeatures} "
        #     f"--nvariables {self.nvariables} "
        #     f"--ndims {self.ndims} "
        #     f"--resolution {self.resolution} "
        # )
        # print(f'\n{get_formatted_time()}\n'
        #     f'Performing matrix QC.')
        # subprocess.check_call(matrix_cmd, shell=True)
        # run perform_matrix_qc
        perform_matrix_qc(matrixdir, anadir, self.minfeatures, 
                        self.mincells, self.nvariables,
                        self.ndims, self.resolution)

        (Path(self.outdir) / ".analysis.done").touch()


def analysis_app(
    # 必需参数
    name: Annotated[str, typer.Option("--name", "-n", help="Sample name", prompt=True, show_default=False)],
    # 可选参数
    outdir: Annotated[str, typer.Option("--outdir", "-o", help=f"Output directory")] = os.getcwd(),
    minfeatures: Annotated[int, typer.Option("--minfeatures","-f", help="Min features per cell")] = 5,
    mincells: Annotated[int, typer.Option("--mincells", "-c", help="Min cells per gene")] = 3,
    ndims: Annotated[int, typer.Option("--ndims", "-d", help="PCA dimensions")] = 30,
    nvariables: Annotated[int, typer.Option("--nvariables", "-v", help="Number of variable genes")] = 2000,
    resolution: Annotated[float, typer.Option("--resolution", "-r", help="Clustering resolution")] = 0.5,
):
    """
    Perform matrix analysis.
    """
    # 将参数转换为类似argparse的Namespace对象
    class Args:
        pass
    args = Args()
    for k, v in locals().items():
        if k != 'args':
            setattr(args, k, v)
    
    # 执行处理流程
    processor = Analysis(args)
    processor.run()

# 导出函数
__all__ = ["analysis_app"]
