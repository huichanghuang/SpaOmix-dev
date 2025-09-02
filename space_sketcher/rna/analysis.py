import os
from pathlib import Path
from typing_extensions import Annotated
import typer
from loguru import logger

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

        ###perform_matrix_qc
        # judge file exits
        matrixdir = os.path.join(self.outdir, "02.oligo/filtered_matrix")
        judgeFilexits(matrixdir)
        ###检查过滤后的细胞数是否大于50
        barcodefile = os.path.join(self.outdir, "02.oligo/filtered_matrix/barcodes.tsv.gz")
        # 读取barcode文件获取细胞数量
        cell_count = 0
        try:
            # 如果是gzip压缩文件
            import gzip
            with gzip.open(barcodefile, 'rt') as f:
                cell_count = sum(1 for _ in f)
        except:
            # 如果不是压缩文件或读取失败，尝试普通方式读取
            try:
                with open(barcodefile, 'r') as f:
                    cell_count = sum(1 for _ in f)
            except Exception as e:
                logger.error(f"Failed to read barcode file: {e}")
                cell_count = 0

        logger.info(f"Filtered cell count: {cell_count}")
        # 判断细胞数是否大于50，否则不运行matrix分析
        if cell_count > 50:
            logger.info("Processing matrix QC...")
            perform_matrix_qc(matrixdir, anadir, self.minfeatures, 
                            self.mincells, self.nvariables,
                            self.ndims, self.resolution)
            
            (Path(self.outdir) / ".analysis.done").touch()
        else:
            logger.warning(f"Filtered cell count ({cell_count}) is less than or equal to 50. Skipping matrix analysis.")
            # 可以选择创建一个标记文件表示分析被跳过
            (Path(self.outdir) / ".analysis.skipped").touch()
        # logger.info("Processing matrix QC...")
        # perform_matrix_qc(matrixdir, anadir, self.minfeatures, 
        #                 self.mincells, self.nvariables,
        #                 self.ndims, self.resolution)

        # (Path(self.outdir) / ".analysis.done").touch()


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
