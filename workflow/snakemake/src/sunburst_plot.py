import pandas as pd
import gzip
import plotly.express as px
import plotly.graph_objects as go
from collections import defaultdict
import argparse

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--cellreads", type=str, required=True)
    parser.add_argument("--filtered_barcodes", type=str, required=True)
    parser.add_argument("--clustered_barcodes", type=str, required=True)
    parser.add_argument("--summaryfile", type=str, required=True)
    parser.add_argument("--sb_umis", type=str, required=True)
    parser.add_argument("--coordfile", type=str, required=True)
    parser.add_argument("--oligochip", type=str, default="LD")
    parser.add_argument("--output", type=str, default="sunburst_plots.html")
    return parser.parse_args()

comp_tab = str.maketrans('ATCG', 'TAGC')
def reverse_complement(seq):
    return seq.translate(comp_tab)[::-1]

def get_barcode(_file):
    barcodes = set()
    if _file.endswith("gz"):
        with gzip.open(_file) as inf:
            for line in inf:
                barcodes.add(line.strip().decode("utf-8"))
    else:
        with open(_file) as inf:
            for line in inf:
                barcodes.add(line.strip())        

    return barcodes

def extract_rna_statistic(_rnainfile, _filtered_barcodes, _clustered_barcodes):
    indata = pd.read_csv(_rnainfile, header=0, sep = "\t")
    totalreads = indata['cbMatch'].sum()

    validdata = indata[indata['CB']!='CBnotInPasslist']
    valid_reads = validdata['cbMatch'].sum()
    valid_mapreads_genome = validdata['genomeU'].sum()+validdata['genomeM'].sum()
    valid_unimap_genome = validdata['genomeU'].sum()
    valid_unimap_transcriptome = validdata['featureU'].sum()
    valid_reads_in_filtered_cells = validdata[validdata['CB'].isin(_filtered_barcodes)]['countedU'].sum()
    valid_reads_in_clustered_cells = validdata[validdata['CB'].isin(_clustered_barcodes)]['countedU'].sum()
    unique_rna_reads = validdata[validdata['CB'].isin(_clustered_barcodes)]['nUMIunique'].sum()

    labels = ["Total RNA Reads", "Valid Barcode Reads", "Invalid Barcode Reads",
            "Unique Mapped Reads", "Multi-mapped Reads", "Unmapped Reads",
            "Transcriptome Mapped Reads", "Unannotated Reads",
            "Reads in Filtered Cells", "Reads in Background",
            "Reads in Clustered Cells", "Unclustered Reads",
            "Unique Reads", "Saturation"]
    parents = ["", "Total RNA Reads", "Total RNA Reads",
            "Valid Barcode Reads", "Valid Barcode Reads", "Valid Barcode Reads",
            "Unique Mapped Reads", "Unique Mapped Reads",
            "Transcriptome Mapped Reads", "Transcriptome Mapped Reads",
            "Reads in Filtered Cells", "Reads in Filtered Cells",
            "Reads in Clustered Cells", "Reads in Clustered Cells"]
    values = [totalreads, valid_reads, totalreads-valid_reads,
            valid_unimap_genome, valid_mapreads_genome-valid_unimap_genome, valid_reads-valid_mapreads_genome, 
            valid_unimap_transcriptome, valid_unimap_genome-valid_unimap_transcriptome,
            valid_reads_in_filtered_cells, valid_unimap_transcriptome-valid_reads_in_filtered_cells,
            valid_reads_in_clustered_cells, valid_reads_in_filtered_cells-valid_reads_in_clustered_cells,
            unique_rna_reads, valid_reads_in_clustered_cells-unique_rna_reads]

    return labels, parents, values


def extract_oligo_statistic(_summaryfile, _sb_umis, _coordfile, _clustered_barcode, _oligochip):
    summary = pd.read_csv(_summaryfile, header=None)
    summary.columns = ["key", "value"]
    total_spatial_reads = int(float(summary[summary["key"]=="Total Spatial Reads"]['value'].iloc[0]))
    reads_with_valid_cellbarcode =  int(float(summary[summary["key"]=="Spatial Reads with Valid Cellbarcode"]['value'].iloc[0]))
    valid_spatial_reads = int(float(summary[summary["key"]=="Valid Spatial Reads"]['value'].iloc[0]))
    valid_spatial_reads_in_filtered_cells = int(float(summary[summary["key"]=="Valid Spatial Reads in Cells"]['value'].iloc[0]))
    reads_in_filtered_cells_with_location = int(float(summary[summary["key"]=="Spatial Reads in Cells with Location on Chip"]['value'].iloc[0]))

    df_umis = pd.read_csv(f"{_sb_umis}", compression="gzip")
    df_umis = df_umis.dropna(subset=["Spatial_Barcode"])
    truecell_umis = df_umis[df_umis["Cell_Barcode"].isin(_clustered_barcode)]

    if _oligochip == "LD":
        ###spatial barcode截取, 为了匹配puckfile中的barcode
        truecell_umis["SUB_SB"] = truecell_umis["Spatial_Barcode"].str[:10]+truecell_umis["Spatial_Barcode"].str[12:18]
    elif _oligochip == "GM":
        truecell_umis["SUB_SB"] = truecell_umis['Spatial_Barcode'].apply(reverse_complement)
    else:
        truecell_umis["SUB_SB"] = truecell_umis["Spatial_Barcode"]

    coordf = pd.read_csv(_coordfile, header=0)
    coordf.columns = ["SUB_SB", "x", "y"]
    mergedf = pd.merge(truecell_umis, coordf, on="SUB_SB", how="inner")
    reads_in_clustered_cells = int(mergedf["Read_Count"].sum())
    mergedf_dropped = mergedf.drop("Read_Count", axis=1)
    mergedf_dropped["UMI_count"] = mergedf_dropped.groupby(['Cell_Barcode', 'SUB_SB'])['UMI'].transform('nunique')
    df_sorted = mergedf_dropped.sort_values(by=['Cell_Barcode', 'UMI_count'], ascending=[True, False])
    df_sorted = df_sorted.drop("UMI", axis=1)
    df_sorted_dedup = df_sorted.drop_duplicates()
    df_sorted_dedup.columns = ["cb", "sb", "subsb", "xcoord", "ycoord", "umi_count"]
    unique_reads = int(df_sorted_dedup["umi_count"].sum())

    labels = ["Total Spatial Reads", "Reads with Valid Cellbarcode", "Invalid Cell Barcode Reads",
            "Valid Spatial Reads", "Invalid Spatial Reads",
            "Reads in Filtered Cells", "Reads in Background",
            "Reads with Location", "Reads without location",
            "Reads in Clustered Cells", "Unclustered Reads",
            "Unique Reads", "Saturation"]
    parents = ["", "Total Spatial Reads", "Total Spatial Reads",
            "Reads with Valid Cellbarcode", "Reads with Valid Cellbarcode", 
            "Valid Spatial Reads","Valid Spatial Reads",
            "Reads in Filtered Cells", "Reads in Filtered Cells",
            "Reads with Location", "Reads with Location",
            "Reads in Clustered Cells", "Reads in Clustered Cells"]
    values = [total_spatial_reads, reads_with_valid_cellbarcode, total_spatial_reads-reads_with_valid_cellbarcode,
            valid_spatial_reads, reads_with_valid_cellbarcode-valid_spatial_reads,
            valid_spatial_reads_in_filtered_cells, valid_spatial_reads-valid_spatial_reads_in_filtered_cells,
            reads_in_filtered_cells_with_location, valid_spatial_reads_in_filtered_cells-reads_in_filtered_cells_with_location,
            reads_in_clustered_cells, reads_in_filtered_cells_with_location-reads_in_clustered_cells,
            unique_reads, reads_in_clustered_cells-unique_reads]

    return labels, parents, values

def sunburst_plot(labels, parents, values, title):
    # 计算层次深度和同一层次的扇区计数
    levels = []
    sibling_count = defaultdict(int)  # 记录每个母节点的子节点数量
    sibling_index = {}  # 记录每个标签在其姐妹节点中的索引

    for i, (label, parent) in enumerate(zip(labels, parents)):
        level = 0
        current = parent
        while current != "":
            level += 1
            if current in labels:
                current = parents[labels.index(current)]
            else:
                break
        levels.append(level)
        # 记录姐妹节点信息
        sibling_count[parent] += 1
        sibling_index[label] = sibling_count[parent] - 1

    max_level = max(levels)

    def adjust_brightness_simple(color, factor):
        """简化版亮度调整"""
        color = color.lstrip('#')
        rgb = tuple(int(color[i:i+2], 16) for i in (0, 2, 4))
        new_rgb = [min(255, max(0, int(channel * factor))) for channel in rgb]
        return '#{:02x}{:02x}{:02x}'.format(*new_rgb)

    base_colors = [
        '#FF9AA2',  # 粉红色
        '#FFB7B2',  # 浅珊瑚
        '#FFDAC1',  # 桃色
        '#E2F0CB',  # 淡绿
        '#B5EAD7',  # 薄荷色
        '#C7CEEA',  # 淡紫
        '#F8B195',  # 杏色
        '#F67280',  # 玫瑰红
        '#C06C84',  # 紫红色
        '#6C5B7B',  # 淡紫色
        '#355C7D',  # 钢蓝色
        '#99B898'   # 淡青绿
    ]
    # 为每个标签分配颜色
    colors = []
    for i, (label, level) in enumerate(zip(labels, levels)):
        base_color = base_colors[level % len(base_colors)]
        # 从内到外的基础亮度
        inner_position = (max_level - level) / max_level if max_level > 0 else 0
        base_brightness = 0.9 - (0.1 * inner_position)
        # 同一层次内不同扇区的亮度变化
        parent = parents[i]
        total_siblings = sibling_count[parent]
        current_index = sibling_index[label]
    
        if total_siblings > 1:
            # 为姐妹节点添加亮度变化
            sibling_variation = (current_index / (total_siblings - 1)) * 0.2  # 最大0.3的亮度变化
            brightness_factor = base_brightness + sibling_variation
        else:
            brightness_factor = base_brightness
        
        color = adjust_brightness_simple(base_color, brightness_factor)
        colors.append(color)

    # 创建太阳图
    fig = go.Figure(go.Sunburst(
        labels=labels,
        parents=parents,
        values=values,
        branchvalues="total",
        marker=dict(
            colors=colors,
            line=dict(width=2, color='white')
        ),
        hovertemplate='<b>%{label}</b><br>数量: %{value:,}<br>占比: %{percentParent:.2%}<extra></extra>'
    ))

    fig.update_layout(
        margin=dict(t=50, l=0, r=0, b=0),
        title={
            'text': title,
            'x': 0.5,
            'xanchor': 'center',
            'font': {'size': 16}
        },
        width=600,
        height=600
    )
    
    return fig

def create_combined_html(fig1, fig2, output_file):
    """将两个图表组合到一个HTML文件中，左右排布 - 精确边框版本"""
    from plotly.subplots import make_subplots
    import plotly.io as pio
    
    # 创建组合图表
    combined_fig = make_subplots(
        rows=1, cols=2,
        specs=[[{"type": "sunburst"}, {"type": "sunburst"}]],
        subplot_titles=(
            "<b>📊 RNA Sequencing Data</b>", 
            "<b>📍 Spatial Sequencing Data</b>"
        ),
        horizontal_spacing=0.08
    )
    
    # 添加第一个太阳图
    for trace in fig1.data:
        combined_fig.add_trace(trace, row=1, col=1)
    
    # 添加第二个太阳图
    for trace in fig2.data:
        combined_fig.add_trace(trace, row=1, col=2)
    
    # 更新布局
    combined_fig.update_layout(
        title_text="<b>RNA and Spatial Sequencing Data Quality Analysis</b>",
        title_x=0.5,
        title_font=dict(size=20, family="Arial, sans-serif"),
        height=750,
        width=1300,
        showlegend=False,
        paper_bgcolor='#f8f9fa',
        # 为整个图表区域添加边框
        margin=dict(l=50, r=50, t=80, b=50),
    )
    
    # 使用正确的坐标引用为每个子图添加边框
    # 左边图表边框
    combined_fig.add_shape(
        type="rect",
        xref="paper", yref="paper",
        x0=0.02, x1=0.48, y0=0.05, y1=0.95,
        line=dict(width=2, color="#3498db"),
        fillcolor="rgba(52, 152, 219, 0.05)",
    )
    
    # 右边图表边框
    combined_fig.add_shape(
        type="rect",
        xref="paper", yref="paper",
        x0=0.52, x1=0.98, y0=0.05, y1=0.95,
        line=dict(width=2, color="#e74c3c"),
        fillcolor="rgba(231, 76, 60, 0.05)",
    )
    
    # 更新子图标题样式
    combined_fig.update_annotations(
        font=dict(size=14, color='white', family="Arial, sans-serif", weight="bold"),
        bordercolor='#2c3e50',
        borderwidth=2,
        borderpad=6,
        bgcolor='#2c3e50',
        xanchor='center'
    )
    
    # 添加分隔线
    combined_fig.add_shape(
        type="line",
        xref="paper", yref="paper",
        x0=0.5, x1=0.5, y0=0, y1=1,
        line=dict(width=2, color="#95a5a6", dash="dot")
    )
    
    # 保存为HTML文件
    pio.write_html(combined_fig, file=output_file, auto_open=False)
    print(f"图表已保存到: {output_file}")

if __name__ == "__main__":
    args = parse_args()
    
    # 获取barcode数据
    filtered_barcodes = get_barcode(args.filtered_barcodes)
    clustered_barcodes = get_barcode(args.clustered_barcodes)
    
    # 提取RNA统计数据
    rna_labels, rna_parents, rna_values = extract_rna_statistic(
        args.cellreads, filtered_barcodes, clustered_barcodes
    )
    
    # 提取Spatial统计数据
    spatial_labels, spatial_parents, spatial_values = extract_oligo_statistic(
        args.summaryfile, args.sb_umis, args.coordfile, clustered_barcodes, args.oligochip
    )
    
    # 创建两个太阳图
    rna_fig = sunburst_plot(rna_labels, rna_parents, rna_values, "RNA Sequencing Data Quality")
    spatial_fig = sunburst_plot(spatial_labels, spatial_parents, spatial_values, "Spatial Data Quality")
    
    # 将两个图表组合并保存到HTML文件
    create_combined_html(rna_fig, spatial_fig, args.output)
    