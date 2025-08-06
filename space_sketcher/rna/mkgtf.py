import polars as pl
import re
import os
import typer
from typing_extensions import Annotated
from typing import Optional

REQUIRED_COLUMNS = [
    "seqname",
    "source",
    "feature",
    "start",
    "end",
    "score",
    "strand",
    "frame",
    "attribute",
]

def read_header(gtf_file):
    """
    读取 gtf 头部注释
    """
    header = []
    with open(gtf_file, "r") as f:
        for line in f:
            if line.startswith("##"):
                header.append(line)
            else:
                break
    return header


def parse_attributes(attr_str: str) -> dict:
    """
    >>> attr_str = 'gene_id "ENSG00000290825.1"; gene_type "lncRNA"; gene_name "DDX11L2"; level 2; tag "overlaps_pseudogene";'
    >>> parse_attributes(attr_str)
    """
    # to do: 正则解析是否过慢，如何加快
    pattern = re.compile(r'(\w+)\s+["\']([^"]*)["\']')
    attributes = pattern.findall(attr_str)
    return dict(attributes)


class Gtf:
    def __init__(self, args):
        self.ingtf = os.path.abspath(args.ingtf)
        self.output = os.path.abspath(args.output)
        self.filter_attrs = args.filter_attrs

    def run(self):
        schema = pl.Schema(
            [
                ("seqname", pl.String),
                ("source", pl.String),
                ("feature", pl.String),
                ("start", pl.Int64),
                ("end", pl.Int64),
                ("score", pl.String),
                ("strand", pl.String),
                ("frame", pl.String),
                ("attribute", pl.String),
            ]
        )
        df = pl.read_csv(
            self.ingtf,
            separator="\t",
            comment_prefix="#",
            has_header=False,
            new_columns=REQUIRED_COLUMNS,
            schema=schema,
        )

        # 按照 gene 分组
        df = df.with_columns(
            [pl.when(pl.col("feature") == "gene").then(1).otherwise(0).alias("group_id")]
        ).with_columns([pl.col("group_id").cum_sum()])

        # 从命令行中构建 polars 过滤表达式
        filter_condition = pl.lit(True)

        if isinstance(self.filter_attrs, str):
            self.filter_attrs = [self.filter_attrs]
        for item in self.filter_attrs:
            item = item.split(":")
            k = item[0].strip()
            v = item[1].strip().split(",")
            filter_condition = filter_condition & pl.col(k).is_in(v)

        def _filter_group(group_df, filter_condition):
            """
            >>> filter_expr = {"gene_type": ["lncRNA","antisense"]}
            """
            sub_df = group_df.filter(pl.col("feature") == "gene")
            attribute = sub_df["attribute"].item()
            attribute = parse_attributes(attribute)
            attr_df = pl.DataFrame(attribute)

            attr_df = attr_df.filter(filter_condition)
            if attr_df.shape[0] != 0:
                return group_df
            else:
                return pl.DataFrame(schema=sub_df.schema)

        res = []
        import tqdm
        total = df['group_id'].n_unique()
        for group_df in tqdm.tqdm(df.group_by("group_id", maintain_order=True),total=total):
            group_df = group_df[1]
            tmp_df = _filter_group(group_df,filter_condition)
            res.append(tmp_df)
        res = pl.concat(res).drop("group_id")

        header = read_header(self.ingtf)
        with open(self.output, "w") as f:
            for line in header:
                f.write(line)
        with open(self.output, "a") as f:
            res.write_csv(f, separator="\t", include_header=False, quote_style="never")


def mkgtf_app(
    # 必需参数
    ingtf: Annotated[
        str, typer.Option(..., "-i", "--ingtf", help="Path to input gtf file.")
    ],
    output: Annotated[
        str, typer.Option(..., "-o", "--output", help="Path to output gtf file.")
    ],
    # 可选参数
    filter_attrs: Annotated[
        str,
        typer.Option(
            ...,
            "-f",
            "--filter_attrs",
            help="Filter attributes as key:value pairs, separated by commas.",
        ),
    ] = "gene_type:protein_coding,lncRNA,antisense,IG_C_gene,IG_D_gene,IG_J_gene"
    ",IG_LV_gene,IG_V_gene,IG_V_pseudogene,IG_J_pseudogene,IG_C_pseudogene"
    ",TR_C_gene,TR_D_gene,TR_J_gene,TR_V_gene,TR_V_pseudogene,TR_J_pseudogene",
):
    """
    Filter GTF file.\n
    
    Example:\n
        space-sketcher rna mkgtf --ingtf /path/to/ingtf --output /path/to/outputgtf
    """
    # 将参数转换为类似argparse的Namespace对象
    class Args:
        pass
    args = Args()
    for k, v in locals().items():
        if k != 'args':
            setattr(args, k, v)
    
    # 执行处理流程
    processor = Gtf(args)
    processor.run()

# 导出函数
__all__ = ["mkgtf_app"]