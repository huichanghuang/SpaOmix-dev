import os
import copy
import argparse
import pandas as pd
from pathlib import Path
from jinja2 import FileSystemLoader, Environment
from space_sketcher.tools.utils import csv2dict, judgeFilexits
from space_sketcher.__init__ import __root_dir__, __version__
from space_sketcher.rna.src.plotly_to_report import (
    plot_barcode_ranks,
    saturation_plot,
    distribution_violin,
    spatial_scatter,
)

def base64uri(filepath):
    """
    将文件转换为 base64 URI 格式   
    参数:
        filepath: 文件路径
    返回:
        base64 URI 字符串，格式为: data:[MIME类型];base64,[base64编码数据]
    """
    import mimetypes
    import base64
    # 猜测文件类型
    mime_type = mimetypes.guess_type(filepath)[0] or 'application/octet-stream'
    # 读取文件并编码为base64
    with open(filepath, "rb") as f:
        encoded_data = base64.b64encode(f.read()).decode('utf-8')
    
    return f"data:{mime_type};base64,{encoded_data}"

def help_icon():
    """
    帮助图标
    """
    icon = """<svg t="1687338260183" class="icon" viewBox="0 0 1024 1024" version="1.1" xmlns="http://www.w3.org/2000/svg" p-id="2757" width="18" height="18" data-spm-anchor-id="a313x.7781069.0.i9"><path d="M512 0C229.23 0 0 229.23 0 512s229.23 512 512 512 512-229.23 512-512S794.77 0 512 0zM512 928c-229.75 0-416-186.25-416-416S282.25 96 512 96s416 186.25 416 416S741.75 928 512 928z" p-id="2758" fill="#707070"></path><path d="M537.64 343.452c47.074 0 83.266-37.528 83.266-78.072 0-32.46-20.832-60.878-62.496-60.878-54.816 0-82.178 44.618-82.178 77.11C475.144 320.132 498.152 343.452 537.64 343.452z" p-id="2759" fill="#707070"></path><path d="M533.162 728.934c-7.648 0-10.914-10.136-3.264-39.55l43.25-166.406c16.386-60.848 10.944-100.398-21.92-100.398-39.456 0-131.458 39.83-211.458 107.798l16.416 27.392c25.246-17.256 67.906-34.762 77.792-34.762 7.648 0 6.56 10.168 0 35.508l-37.746 158.292c-23.008 89.266 1.088 109.538 33.984 109.538 32.864 0 117.808-30.47 195.57-109.632l-18.656-25.34C575.354 716.714 543.05 728.934 533.162 728.934z" p-id="2760" fill="#707070"></path></svg>"""
    return icon

def help_collapse(help_msg: list | str | None = None, _counter=[0]):
    """
    生成折叠的帮助文档，自动生成唯一ID
    参数:
        help_msg: 帮助信息，可以是字符串、列表（支持混合str和dict[name,value]）
        _counter: 内部使用的计数器（无需手动设置）
    
    返回:
        {"icon": 可点击的折叠图标HTML, "text": 折叠内容HTML}
    """
    if help_msg is None:
        return {"icon": "", "text": ""}
    # 生成唯一ID（闭包实现计数器）
    _counter[0] += 1
    collapse_id = f"collapse-{_counter[0]}"
    
    # 处理帮助信息（支持字符串或列表）
    if isinstance(help_msg, list):
        help_msg_str = ""
        for item in help_msg:
            if isinstance(item, str):
                help_msg_str += f"<p>{item}</p>"
            else:
                help_msg_str += f'<h6><strong>{item["name"]}</strong></h6><p>{item["value"]}</p>'
    else:
        help_msg_str = help_msg
    
    # 返回HTML结构
    return {
        "icon": f'<a data-bs-toggle="collapse" href="#{collapse_id}" aria-expanded="false" aria-controls="help">{help_icon()}</a>',
        "text": f'<div class="collapse" id="{collapse_id}">{help_msg_str}</div>'
    }

def to_bootstrap_table(values: dict):
    """
    表格转 html
    """
    table_str = ""
    for k, v in values.items():
        table_str += f'<tr><td class="gt_row gt_left">{k}</td><td class="gt_row gt_right">{v}</td></tr>'

    return f"""
    <table class="gt_table">
    <tbody class="gt_table_body">
        {table_str}
    </tbody>
    </table>
    """

def report_card(
    title: str = None, help_msg: list | None = None, content: str = "", style: str = ""
):
    """
    生成报告的卡片

    Args:
        title (str): 标题
        help_msg (list | None): 帮助文档
        content (str): html 的内容，可以是表格，图片等
    """
    help_msg = help_collapse(help_msg)
    if title is None:
        id = ""
        title = ""
    else:
        id = f'id="{title}"'
        title = f'<div class = "card-header fw-bold">{title}{help_msg["icon"]}</div>'

    html_str = f"""
<div class="card" {id} style="{style}">
    {title}
    <div class="card-body">
        {help_msg["text"]}
        {content}
    </div>
</div>
"""
    return html_str

def get_resource(cdn=False):
    """获取css 以及 js等资源"""
    if cdn:
        links = [
            "https://cdn.jsdelivr.net/npm/bootstrap@5.3.3/dist/css/bootstrap.min.css",
            "https://cdn.jsdelivr.net/npm/@jstable/jstable@1.6.5/dist/jstable.min.css"
        ]
        scripts = [
            "https://cdn.plot.ly/plotly-3.0.0.min.js",
            "https://cdn.jsdelivr.net/npm/bootstrap@5.3.3/dist/js/bootstrap.min.js",
            "https://cdn.jsdelivr.net/npm/@jstable/jstable@1.6.5/dist/jstable.min.js"
        ]
    else:
        path = os.path.join(__root_dir__, "data/template/libs")
        links = [
            base64uri(os.path.join(path, "bootstrap/bootstrap.min.css")),
            base64uri(os.path.join(path, "jstable.min.css")),
        ]

        import plotly

        plotly_js = os.path.join(Path(plotly.__file__).parent, "package_data/plotly.min.js")
        # 轻量
        scripts = [
            base64uri(plotly_js),
            base64uri(os.path.join(path, "bootstrap/bootstrap.min.js")),
            base64uri(os.path.join(path, "jstable.min.js")),
        ]
    return links, scripts

def to_table_html(values: dict):
    """
    表格转 html
    """
    table_str = ""
    for k, v in values.items():
        table_str += f'<tr><td class="gt_row gt_left">{k}</td><td class="gt_row gt_right">{v}</td></tr>'

    return f"""
    <table class="gt_table">
    <tbody class="gt_table_body">
        {table_str}
    </tbody>
    </table>
    """

def tbl(df: pd.DataFrame, id=None):
    """
    数据框转 html
    """
    if id is None:
        id = ""
    else:
        id = f'id="{id}"'
    t_body = ""
    for i, rows in df.iterrows():
        td_str = ""
        for value in rows:
            td_str += f'<td class="gt_row gt_left">{value}</td>'
        t_body += f"<tr>{td_str}</tr>"
    t_head = ""
    for column in df.columns:
        t_head += f'<th class="gt_left" scope="col">{column}</th>'
    html_str = f"""
        <table {id} class="gt_table">
        <thead>
        <tr class="header gt_col_headings">
        {t_head}
        </tr>
        </thead>
        <tbody class="gt_table_body">
        {t_body}
        </tbody>
        </table>
    """
    return html_str


def extract_qc_from_cellread_stat(cellread_stat_file):
    """
    从 STARsolo 输出文件 CellReads.stats 提取 QC 信息
    """
    stats = {}
    df = pd.read_csv(cellread_stat_file, sep="\t")
    stats["Reads Mapped to Genome"] = df["genomeU"].sum() + df["genomeM"].sum()
    stats["Reads Mapped Confidently to Genome"] = df["genomeU"].sum()
    stats["Reads Mapped Confidently to Transcriptome"] = df["featureU"].sum()
    stats["Reads Mapped Confidently to Exonic"] = df["exonic"].sum()
    stats["Reads Mapped Confidently to Intronic"] = df["intronic"].sum()
    # 比对到反义链
    stats["Reads Mapped to antisense"] = df["exonicAS"].sum() + df["intronicAS"].sum()
    stats["Reads Mapped Confidently to Intergenic"] = (
        stats["Reads Mapped to Genome"]
        - stats["Reads Mapped Confidently to Exonic"]
        - stats["Reads Mapped Confidently to Intronic"]
        - stats["Reads Mapped to antisense"]
    )
    stats["Reads Mapped to mito"] = df["mito"].sum()
    for k, v in stats.items():
        stats[k] = v / df["cbMatch"].sum()
    return stats

def get_metrics(indir: Path|str, sample, kit, reference):
    """
    Generate Total summary metrix
    """
    summary = {"Sample": sample, 
               "Reference": reference,
               "Kit": kit,
               "Pipeline version": __version__}
    ###old rna summary generated by starsolo
    old_rna_summary_file = os.path.join(indir, "01.count/Solo.out/GeneFull_Ex50pAS/Summary.csv")
    ###new rna summary after dbscan filtering
    new_rna_summary_file = os.path.join(indir, "02.oligo/dbscan_filtered_cells.summary.csv")
    ###spatial library summary generated by spatial_barcode_extraction and assign_coordinate
    spatial_summary_file = os.path.join(indir, "02.oligo/sb_library_summary.csv")
    
    judgeFilexits(
        old_rna_summary_file,
        new_rna_summary_file,
        spatial_summary_file
    )

    old_rna_summary_ori = csv2dict(old_rna_summary_file)
    old_rna_summary_ori["Total Gene Detected"] = old_rna_summary_ori.pop("Total GeneFull_Ex50pAS Detected")
    retain_keys = ["Number of Reads", "Reads With Valid Barcodes", 
                              "Sequencing Saturation", "Q30 Bases in CB+UMI", 
                              "Q30 Bases in RNA read", "Estimated Number of Cells",
                              "Total Gene Detected"]
    old_rna_summary = {k: old_rna_summary_ori[k] for k in retain_keys if k in old_rna_summary_ori}
    new_rna_summary = csv2dict(new_rna_summary_file)
    spatial_summary = csv2dict(spatial_summary_file)

    ###merge
    summary.update(old_rna_summary)
    summary.update(new_rna_summary)
    summary.update(spatial_summary)

    # Mapping
    cellread_stat_file = (
        os.path.join(os.path.dirname(old_rna_summary_file), "CellReads.stats")
    )
    stats = extract_qc_from_cellread_stat(cellread_stat_file)
    summary.update(stats)

    out_summary = os.path.join(indir, "04.report/summary.csv")

    # 转换为 DataFrame
    outdf = pd.DataFrame.from_dict(summary, orient='index', columns=['Value'])
    outdf.index.name = 'Metric'  # 设置索引列的名称（可选）
    outdf.reset_index(inplace=True)  # 将索引转为列
    outdf.to_csv(out_summary, header=False, index=False)


def format_summary(summary: dict):
    """
    格式化 summary，返回一个新的字典
    """
    import re
    formatted_summary = summary.copy()

    def convert_to_number(value):
        if isinstance(value, str):
            try:
                # 尝试去掉千分位逗号（如 "1,000" → "1000"）
                cleaned = value.replace(",", "")
                # 如果能转为整数，则返回 int，否则返回 float
                return int(cleaned) if cleaned.isdigit() else float(cleaned)
            except ValueError:
                return value  # 转换失败则保持原样
        return value  # 非字符串直接返回
    
    formatted_summary = {k: convert_to_number(v) for k, v in formatted_summary.items()}

    int_keys = [
        "Number of Reads",
        "Total Gene Detected",
        "Estimated Number of Cells",
        "Cells with Certain Location",
        "Unique Reads in Cells Mapped to Gene",
        "Total Spatial Reads",
        "Spatial Reads with Valid Cellbarcode",
        "Valid Spatial Reads",
        "Valid Spatial UMIs",
        "Valid Spatial Reads in Cells",
        "Valid Spatial UMIs in Cells",
        "Total Spatial Barcodes with Location on Chip",
        "Unique Valid Spatial Barcodes with Location on Chip",
        "Spatial Reads in Cells with Location on Chip",
        *[key for key in formatted_summary.keys() if re.search(r"per Cell$", key)],
    ]
    if "Median top100 Spatial UMI Mean per Cell" in int_keys:
        int_keys.remove("Median top100 Spatial UMI Mean per Cell")
    if "Mean top100 Spatial UMI Mean per Cell" in int_keys:
        int_keys.remove("Mean top100 Spatial UMI Mean per Cell")
        
    for k in int_keys:
        if not isinstance(formatted_summary[k], str):
            formatted_summary[k] = f"{int(formatted_summary[k]):,}"

    pct_keys = [
        "Reads With Valid Barcodes",
        "Sequencing Saturation",
        "Spatial Barcode Saturation",
        *[key for key in formatted_summary.keys() if re.search(r"^Q30|^Reads Mapped|^Fraction of|-cluster$", key)],
    ]
    for k in pct_keys:
        if not isinstance(formatted_summary[k], str):
            formatted_summary[k] = f"{formatted_summary[k]:.2%}"

    return formatted_summary

def generate_report(indir: Path|str, sample, kit, reference, oligochip, dev=True, cdn = False):

    ###prepare summary
    get_metrics(indir, sample, kit, reference)

    ###jugde file exits
    summaryfile = os.path.join(indir, "04.report/summary.csv")
    template_file = os.path.join(__root_dir__, "data/template/count_template.html")
    kneefile1 = os.path.join(indir, "01.count/cell_rna_umi.rank.txt")
    kneefile2 = os.path.join(indir, "02.oligo/cell_sb_umi.rank.txt")
    sequencing_saturation_file = os.path.join(indir, "01.count/saturation.xls")
    metadata = os.path.join(indir, "03.analysis/UMAPpos.txt")
    markergenefile = os.path.join(indir, "03.analysis/top30DEG.txt")

    judgeFilexits(
        summaryfile,
        template_file,
        kneefile1,
        kneefile2,
        sequencing_saturation_file,
        metadata,
        markergenefile
    )

    _summary = csv2dict(summaryfile)    
    template_dir = os.path.join(__root_dir__, "data/template")
    env = Environment(loader=FileSystemLoader(template_dir))
    template = env.get_template("count_template.html")

    n_cells = _summary["Cells with Certain Location"]
    summary = format_summary(_summary)
    ######Summary 页面
    ###顶部指标项
    keys = [
        "Cells with Certain Location",
        "Median Gene per Cell",
        "Mean Reads per Cell",
        "Median UMI per Cell",
    ]
    metric_str = "".join(
        [
            f"""
        <div class="col-6 col-md-3">
            <div class="metric-box d-flex">
                <div class="flex-grow-1">
                    <div class="metric-title">{key}</div>
                    <div class="metric-number">{summary[key]}</div>
                </div>
                <div class="metric-icon">
                    <img src={base64uri(os.path.join(template_dir,"bar-icon.png"))} alt="bar" width="60px">
                </div>
            </div>
        </div>
        """
            for key in keys
        ]
    )

    ###Cells
    keys = [
        "Estimated Number of Cells",
        "Cells with Certain Location",
        "Fraction of Unique Reads in Cells",
        "Fraction of UMIs in Cells",
        "Mean Reads per Cell",
        "Median Reads per Cell",
        "Mean UMI per Cell",
        "Median UMI per Cell",
        "Mean Gene per Cell",
        "Median Gene per Cell",
        "Total Gene Detected",
    ]
    cell_help_msg = {
        "Estimated Number of Cells": "The number of barcodes associated with at least one cell.",
        "Cells with Certain Location": "The number of cell barcodes that can be assigned to a certain location on a chip.",
        "Fraction of Unique Reads in Cells": "The proportion of unique mapped reads assigned to valid cell associated barcodes.",
        "Fraction of UMIs in Cells": "The proportion of UMIs assigned to valid cell associated barcodes.",
        "Mean Reads per Cell": "The average number of mapped reads per cell detected.",
        "Median Reads per Cell": "The median number of unique mapped reads per cell detected.",
        "Mean UMI per Cell": "The average number of UMIs per cell detected.",
        "Median UMI per Cell": "The median number of UMIs per cell detected.",
        "Mean Gene per Cell": "The average number of Genes per cell detected.",
        "Median Gene per Cell": "The median number of Genes per cell detected.",
        "Total Gene Detected": "The number of genes with at least one UMI count in cell detected.",
    }
    if not dev:
        to_remove = ["Fraction of Unique Reads in Cells", "UMIs in Cells"]
        keys = list(filter(lambda x: x not in to_remove, keys))
        temp_cell_help_msg = {key:val for key, val in cell_help_msg.items() if key not in to_remove }
        cell_help_msg = copy.deepcopy(temp_cell_help_msg)

    content = (
        f'<div class="row"><div class = "col-6">'
        f'{to_bootstrap_table({key: summary[key] for key in keys[0:5]})}'
        f'</div><div class = "col-6">'
        f'{to_bootstrap_table({key: summary[key] for key in keys[5:]})}'
        f'</div></div>'
    )
    cell_str = report_card(
        title="Cells",
        help_msg=[{"name": k, "value": v} for (k, v) in cell_help_msg.items()],
        content=content,
    )
 
    ###Spatial library
    keys = [
        "Total Spatial Reads",
        "Spatial Reads with Valid Cellbarcode",
        "Fraction of Valid Spatial Reads",
        "Valid Spatial UMIs",
        "Spatial Barcode Saturation",
        "Fraction of Valid Spatial Reads in Cells",
        "Fraction of Valid Spatial UMIs in Cells",
        "Total Spatial Barcodes with Location on Chip",
        "Fraction of Unique Valid Spatial Barcodes with Location on Chip",
        "Fraction of Spatial Reads in Cells with Location on Chip",
        "Median top100 Spatial UMI Mean per Cell",
        "Mean top100 Spatial UMI Mean per Cell",
        "single-cluster",
        "multi-cluster",
        "no-cluster",
    ]
    spatial_help_msg = {
        "Total Spatial Reads": "The total number of reads in spatial library.",
        "Spatial Reads with Valid Cellbarcode": "The number of spatial reads with a valid cell barcode.",
        "Fraction of Valid Spatial Reads": "The fraction of spatial reads with a valid cell barcode and valid spatial barcode.",
        "Valid Spatial UMIs": "The number of total UMIs with a valid cell barcode and valid spatial barcode.",
        "Spatial Barcode Saturation": "The percentage of spatial reads originating from a duplicate UMI.",
        "Fraction of Valid Spatial Reads in Cells": "The fraction of valid spatial reads in cell associated barcodes.",
        "Fraction of Valid Spatial UMIs in Cells": "The fraction of valid spatial UMIs in cell associated barcodes.",
        "Total Spatial Barcodes with Location on Chip": "The number of spatial barcodes with location on chip.",
        "Fraction of Unique Valid Spatial Barcodes with Location on Chip": "The fraction of spatial barcodes on the chip that are both unique and match the whitelist.",
        "Fraction of Spatial Reads in Cells with Location on Chip": "The fraction of spatial reads in cell associated barcodes and with location on chip.",
        "Median top100 Spatial UMI Mean per Cell": "Median across cells of the mean UMI count in the top 100 spatial barcodes per cell.",
        "Mean top100 Spatial UMI Mean per Cell": "Mean across cells of the mean UMI count in the top 100 spatial barcodes per cell",
        "single-cluster": "The proportion of cell associated barcodes with a single cluster in dbscan clustering.",
        "multi-cluster": "The proportion of cell associated barcodes with multiple cluster in dbscan clustering.",
        "no-cluster": "The proportion of cell associated barcodes without any cluster in dbscan clustering.",
    }
    if not dev:
        to_remove = ["Median top100 Spatial UMI Mean per Cell",
                    "Mean top100 Spatial UMI Mean per Cell",
                    "single-cluster",
                    "multi-cluster",
                    "no-cluster"
        ]
        keys = list(filter(lambda x: x not in to_remove, keys))
        temp_spatial_help_msg = {key:val for key, val in spatial_help_msg.items() if key not in to_remove }
        spatial_help_msg = copy.deepcopy(temp_spatial_help_msg)

    content = (
        f'<div class="row"><div class = "col-6">'
        f'{to_bootstrap_table({key: summary[key] for key in keys[0:7]})}'
        f'</div><div class = "col-6">'
        f'{to_bootstrap_table({key: summary[key] for key in keys[7:]})}'
        f'</div></div>'
    )
    spatial_str = report_card(
        title="Spatial library",
        help_msg=[{"name": k, "value": v} for (k, v) in spatial_help_msg.items()],
        content=content,
    )

    # mapping
    mapping_help_msg = {
        "Reads Mapped to Genome": "The number and fraction of reads mapped to the genome.",
        "Reads Mapped Confidently to Genome": "The number and fraction of reads uniquely mapped to the genome.",
        "Reads Mapped Confidently to Transcriptome": "The number and fraction of reads uniquely mapped to the transcriptome. The read must be consistent with annotated splice junctions.",
        "Reads Mapped Confidently to Exonic": "The number and fraction of reads that uniquely mapped to an exonic region of the genome.",
        "Reads Mapped Confidently to Intronic": "The number and fraction of reads that uniquely mapped to an intronic region of the genome.",
        "Reads Mapped Confidently to Intergenic": "The number and fraction of reads that uniquely mapped to the intergenic region of the genome.",
        "Reads Mapped to antisense": "The number and fraction of reads uniquely mapped to the antisense strand of gene.",
        "Reads Mapped to mito": "Fraction of reads that mapped uniquely to mitochondria region of the genome.",
    }

    keys = [
        "Reads Mapped to Genome",
        "Reads Mapped Confidently to Genome",
        "Reads Mapped Confidently to Transcriptome",
        "Reads Mapped Confidently to Exonic",
        "Reads Mapped Confidently to Intronic",
        "Reads Mapped Confidently to Intergenic",
        "Reads Mapped to antisense",
        "Reads Mapped to mito",
    ]
    content = (
        f'<div class="row"><div class = "col-6">{to_bootstrap_table({key: summary[key] for key in keys[0:4]})}'
        f'</div><div class = "col-6">'
        f"{to_bootstrap_table({key: summary[key] for key in keys[4:8]})}</div></div>"
    )
    mapping_str = report_card(
        title="Mapping",
        help_msg=[{"name": k, "value": v} for (k, v) in mapping_help_msg.items()],
        content=f"{content}",
    )

    # summary
    summary_help_msg = {
        "Number of Reads": "The number of total raw reads.",
        "Reads With Valid Barcodes": "The number and fraction of reads with valid barcodes.",
        "Sequencing Saturation": "Sequencing saturation refers to the percentage of reads originating from a duplicate UMI (in other words, an mRNA molecule that was sequenced more than 1 time).",
        "Q30 Bases in CB+UMI": "The fraction of cell-barcode and UMI bases with Q-scores >= 30.",
        "Q30 Bases in RNA read": "The fraction of RNA read bases with Q-scores >= 30.",
        "Sequencing Saturation Plot": "This plot shows the Sequencing Saturation metric as a function of downsampled sequencing depth, up to the observed sequencing depth. A high saturation indicates we are detecting the vast majority of mRNA molecules in the samples, and thus don't need to sequence the libraries any deeper.",
    }
    keys = [
        "Sample",
        "Reference",
        "Kit",
        "Pipeline version",
        "Number of Reads",
        "Reads With Valid Barcodes",
        "Sequencing Saturation",
        "Q30 Bases in CB+UMI",
        "Q30 Bases in RNA read",
    ]
    content = (
        f'<div class="row"><div class = "col-6">{to_bootstrap_table({key: summary[key] for key in keys[0:5]})}'
        f'</div><div class = "col-6">'
        f"{to_bootstrap_table({key: summary[key] for key in keys[5:10]})}</div></div>"
    )
    summary_str = report_card(
        title="Summary",
        help_msg=[{"name": k, "value": v} for (k, v) in summary_help_msg.items()],
        content=content,
    )

    ###骤降曲线
    knee_fig1 = plot_barcode_ranks(kneefile1, "RNA")
    knee_fig2 = plot_barcode_ranks(kneefile2, "Spatial")
    UMI_rank_str = report_card(
        title="UMI rank",
        help_msg=[
            {
                "name": "RNA UMI(left)",
                "value": "Total RNA UMI count for each cell barcode plotted against its rank. The distribution of total counts exhibits a sharp transition between barcodes with large and small total counts, probably corresponding to cell-containing and empty droplets respectively. The blue part of the curve indicates cell-associated barcodes, while the gray part may be the background-associated barcodes.",
            },
            {
                "name": "Spatial UMI(left)",
                "value": "Total spatial UMI count for each cell barcode plotted against its rank. The blue part of the curve indicates cell-associated barcodes, while the gray part may be the background-associated barcodes.",
            },
        ],
        content=f"""
            <div class="row">
                <div class="col-6">
                    {
                        knee_fig1.to_html(
                            full_html=False,
                            default_height="380px",
                            include_plotlyjs=False,
                        )
                    }
                </div>
                <div class="col-6">
                    {
                        knee_fig2.to_html(
                            full_html=False,
                            default_height="380px",
                            include_plotlyjs=False,
                        )
                    }
                </div>
            </div>
            """,
    )

    ###饱和度曲线
    sat_fig1, sat_fig2 = saturation_plot(sequencing_saturation_file)
    sequencing_saturation_str = report_card(
        title="Sequencing_saturation",
        help_msg=[
            {
                "name": "Left",
                "value": "This plot shows the Sequencing Saturation metric as a function of downsampled sequencing depth (measured in mean reads per cell), up to the observed sequencing depth. Sequencing Saturation is a measure of the observed library complexity, and approaches 1.0 (100%) when all converted mRNA transcripts have been sequenced. The slope of the curve near the endpoint can be interpreted as an upper bound to the benefit to be gained from increasing the sequencing depth beyond this point. The dotted line is drawn at a value reasonably approximating the saturation point.",
            },
            {
                "name": "Right",
                "value": "This plot shows the Median Genes per Cell as a function of downsampled sequencing depth in mean reads per cell, up to the observed sequencing depth. The slope of the curve near the endpoint can be interpreted as an upper bound to the benefit to be gained from increasing the sequencing depth beyond this point.",
            },
        ],
        content=f"""
            <div class="row">
                <div class="col-6">
                    {
                        sat_fig1.to_html(
                            full_html=False,
                            default_height="380px",
                            include_plotlyjs=False,
                        )
                    }
                </div>
                <div class="col-6">
                    {
                        sat_fig2.to_html(
                            full_html=False,
                            default_height="380px",
                            include_plotlyjs=False,
                        )
                    }
                </div>
            </div>
            """,
    )
    
    summary_page = f"""
    <div class="row mb-4">
        <div class="col-12">
            <div class="card">
                <div class="card-body">
                    <div class="row">
                        {metric_str}
                    </div>
                </div>
            </div>
        </div>
    </div>
    <div class="row">
        <div class="col-12 mb-4 col-lg-6">
            {cell_str}
        </div>
        <div class="col-12 mb-4 col-lg-6">
            {spatial_str}
        </div>
        <div class="col-12 mb-4 col-lg-6">
            {UMI_rank_str}
        </div>
        <div class="col-12 mb-4 col-lg-6">
            {sequencing_saturation_str}
        </div>
        <div class="col-12 mb-4 col-lg-6">
            {summary_str}
        </div>
        <div class="col-12 mb-4 col-lg-6">
            {mapping_str}
        </div>
    </div>
    """

    ######Analysis 页面
    # 小提琴图
    fig = distribution_violin(metadata, samplename=summary["Sample"])
    content = fig.to_html(
        full_html=False,
        default_height="450px",
        default_width="950px",
        include_plotlyjs=False,
    )
    violin_str = report_card(
        title="Distribution",
        help_msg=[
            {
                "name": "Gene counts",
                "value": "The distribution of effective gene counts detected in each cell.",
            },
            {
                "name": "UMI counts",
                "value": "The distribution of UMI counts detected in each cell.",
            },
            {
                "name": "Mito percentage",
                "value": "The distribution of mitochondrial fraction detected in each cell.",
            },
        ],
        content=f'<div class="d-flex justify-content-end">{content}</div>',
        style="width: 1000px",
    )

    #UMI distribution on UMAP and Spatial location
    umi_fig1, umi_fig2 = spatial_scatter(metadata, "UMI", oligochip)
    def _to_html(fig):
        return fig.to_html(
            full_html=False,
            default_height="550px",
            default_width="480px",
            include_plotlyjs=False,
        )

    content = (
        f'<div class="row"><div class = "col-6">{_to_html(umi_fig1)}'
        f'</div><div class = "col-6">'
        f"{_to_html(umi_fig2)}</div></div>"
    )
    umi_str = report_card(
        title="UMI Counts",
        help_msg=[
            "The display is limited to a random subset of cells.",
            {
                "name": "left",
                "value": "This plot displays the log2(UMI counts) for each spot barcode, each associated with an exact coordinate on the chip, and is colored by the log2 value of UMI counts.",
            },
            {
                "name": "right",
                "value": "This plot shows the log2(UMI counts) for each spot barcode. Each dot associated with a spot barcode and is colored by the log2 value of UMI counts. The coordinate axes represents the 2-dimensional embedding produced by the UMAP(Uniform Manifold Approximation and Projection)algorithm.",
            },
        ],
        content=content,
        style="width: 1000px",
    )

    #Cluster distribution on UMAP and Spatial location
    cluster_fig1, cluster_fig2 = spatial_scatter(metadata, "Cluster", oligochip)
    content = (
        f'<div class="row"><div class = "col-6">{_to_html(cluster_fig1)}'
        f'</div><div class = "col-6">'
        f"{_to_html(cluster_fig2)}</div></div>"
    )
    cluster_str = report_card(
        title="Cluster",
        help_msg=[
            "The display is limited to a random subset of cells.",
            {
                "name": "left",
                "value": "This plot shows the automated clustering result for each cell-barcode by UMAP algorithm. each associated with an exact coordinate on the chip, and is colored according to different cluster.",
            },
            {
                "name": "right",
                "value": "This plot shows the automated clustering result for each cell-barcode by UMAP algorithm. Each dot associated with a cell barcode and is colored according to different cluster. The coordinate axes represents the 2-dimensional embedding produced by the UMAP(Uniform Manifold Approximation and Projection)algorithm.",
            },
        ],
        content=content,
        style="width: 1000px",
    )

    # marker gene
    marker_genes = pd.read_csv(markergenefile, sep = "\t")
    marker_genes.columns = ["cluster", "gene", "scores", "log2FC", "pvals", "pvals_adj", "pct_nz_group", "pct_nz_reference"]
    marker_genes = (
        marker_genes.groupby("cluster", observed=True)
        .apply(
            lambda x: x.sort_values(by="pvals_adj", ascending=True).head(30),
            include_groups=False,
        )
        .reset_index(level=0)
    )
    # 格式化
    marker_genes["scores"] = marker_genes["scores"].map(lambda x: f"{x:.4f}")
    marker_genes["log2FC"] = marker_genes["log2FC"].map(lambda x: f"{x:.4f}")
    marker_genes["pvals"] = marker_genes["pvals"].map(lambda x: f"{x:.4f}")
    marker_genes["pvals_adj"] = marker_genes["pvals_adj"].map(lambda x: f"{x:.4f}")
    marker_genes["pct_nz_group"] = marker_genes["pct_nz_group"].map(lambda x: f"{x:.4f}")
    marker_genes["pct_nz_reference"] = marker_genes["pct_nz_reference"].map(lambda x: f"{x:.4f}")

    marker_str = report_card(
        title="Genes",
        help_msg=[
            "The table shows the top 30 differentially expressed genes for each cluster. \
            Here a differential expression test was performed between each cluster and the rest of the sample for each feature. \
            The avg_log2FC is an estimate of the log2 ratio of expression in a cluster to that in all other cells. \
            The p_val is a measure of the statistical significance of the expression difference and \
            the p_val_adj is adjusted p-value, based on bonferroni correction using all features in the dataset. \
            pct_nz_group is the percentage of cells within the current cluster/group where the gene is detected (expression value > 0), \
            pct_nz_reference is the percentage of cells outside the target group (i.e., all other cells) where the gene is detected."
        ],
        content=tbl(marker_genes, id="table"),
    )
    js = """<script>
        let myTable = new JSTable("#table", {
            sortable: true,
            searchable: true,
            perPage: 10,
        });
        </script>"""

    analysis_page = f"""
        <div class="row">
            <div class="col mb-4">
                {violin_str}
                <div class="mt-4">
                {umi_str}
                </div>
                <div class="mt-4">
                {cluster_str}
                </div>
            </div>
            <div class="col mb-4">
                {marker_str}
            </div>
        </div>
        """
    analysis_page += js

    links, scripts = get_resource(cdn=cdn)

    content = template.render(
        title="Count",
        logo_uri=base64uri(os.path.join(template_dir, "logo-high.png")),
        favicon=base64uri(os.path.join(template_dir, "favicon.ico")),
        links=links,
        scripts=scripts,
        summary=summary_page,
        analysis=analysis_page,
    )

    outfile = os.path.join(indir, "04.report/summary.html")
    with open(outfile, "w", encoding="utf-8") as f:
        f.write(content)



def parse_args():
    parser = argparse.ArgumentParser(description='Generate report.')
    parser.add_argument('-i', '--indir', 
        metavar='PATH', 
        type=str,
        help='The input directory.'
        )
    parser.add_argument('-d', '--dev', 
        metavar='BOOL', 
        type=bool,
        default=True,
        help='Development mode, selected from [False, True], default=True.'
        )
    parser.add_argument('-s', '--sample', 
        metavar='FILE', 
        type=str,
        help='Sample name.'
        )
    parser.add_argument('-k', '--kit', 
        metavar='FILE', 
        type=str,
        help='The chemistry kit version.'
        )
    parser.add_argument('-r', '--reference', 
        metavar='FILE', 
        type=str,
        help='The mapping reference.'
        )
    parser.add_argument(
        '-o', '--oligochip',
        metavar='STR',
        choices=["LD","GM"],
        help='spatial chip version, only can be "LD" or "GM". [default: LD].',
        default='LD'
        )
    args = parse_args()
    return args
    

if __name__=='__main__':
    args = parse_args()
    generate_report(args.indir, args.sample, args.kit, args.reference, args.oligochip, args.dev)
