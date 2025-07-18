import typer
from .run import run_app
from .count import count_app
from .oligo import oligo_app
from .analysis import analysis_app
from .report import report_app
from .mkref import mkref_app

app = typer.Typer(name="rna", help="Perform quality control, alignment, and filtering using spatial RNA and oligo library sequencing data.")

# 添加子命令
app.command("run")(run_app)
app.command("count")(count_app)
app.command("oligo")(oligo_app)
app.command("analysis")(analysis_app)
app.command("report")(report_app)
app.command("mkref")(mkref_app)

__all__ = ["app"]