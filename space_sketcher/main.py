import sys
from space_sketcher.__init__ import __version__, __category__
import typer
import sys
from pathlib import Path

def main():
    # 确保能正确导入子模块
    sys.path.append(str(Path(__file__).parent))

    app = typer.Typer(name="main", help="Spatial Multi-omics Analysis Workflow Suite. This suite provides a comprehensive set of tools for spatial multi-omics analysis, including RNA sequencing workflows.")

    # 导入rna模块
    try:
        from rna import app as rna_app
        app.add_typer(rna_app, name="rna", help="Spatial RNA analysis workflows.")
    except ImportError as e:
        typer.echo(f"Error loading RNA module: {e}", err=True)
        raise typer.Exit(1)

    return app()

if __name__ == "__main__":
    main()