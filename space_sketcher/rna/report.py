import os,shutil
from pathlib import Path
import typer
from typing_extensions import Annotated
from typing import Optional

class Report:
    def __init__(self,args):
        self.outdir = os.path.abspath(os.path.join(args.outdir,args.name))
        self.name = args.name
        self.kit = args.rnachemistry
        self.reference = args.reference
        self.dev = args.dev
        self.oligochip = args.oligochip

    def run(self):
        ### import lib
        from space_sketcher.tools.utils import str_mkdir
        from space_sketcher.rna.src.generate_report import generate_report

        ###generate report
        reportdir = os.path.join(self.outdir, "04.report")
        str_mkdir(reportdir) 
        ### run
        generate_report(
            self.outdir,
            self.name,
            self.kit,
            self.reference,
            self.oligochip,
            self.dev
        )
        
        ###copy files
        os.system(f"cp -r {self.outdir}/02.oligo/filtered_matrix {reportdir}/SCST")
        (Path(self.outdir) / ".report.done").touch()


def report_app(
    # 必需参数
    name: Annotated[str, typer.Option("--name", "-n", help="Sample name", prompt=True, show_default=False)],       
    # 可选参数 
    outdir: Annotated[str, typer.Option("--outdir", "-o", help=f"Output directory")] = os.getcwd(),
    rnachemistry: Annotated[str, typer.Option("--rnachemistry", "-rc", help="Chemistry version: 10X/leader_v1/other")] = "leader_v1",
    reference: Annotated[Optional[str], typer.Option("--reference", help="Reference name (default: genomeDir basename)", show_default=False)] = None,    
    dev: Annotated[bool, typer.Option("--dev", help="Enable development mode (default: False)")] = False,
    oligochip: Annotated[str, typer.Option("--oligochip", "-oc", help="Spatial chip version: LD/GM")] = "LD",
):
    """
    Generate report.
    """
    # 将参数转换为类似argparse的Namespace对象
    class Args:
        pass
    args = Args()
    for k, v in locals().items():
        if k != 'args':
            setattr(args, k, v)
    
    # 执行处理流程
    processor = Report(args)
    processor.run()

# 导出函数
__all__ = ["report_app"]
