import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import numpy as np
import math, collections
from plotly.subplots import make_subplots
from scipy.interpolate import make_interp_spline

def cross_plot(crossfile):
    crossdf = pd.read_csv(crossfile, header=0, sep="\t")
    fig = px.scatter(crossdf, x = "column_sum_in_matrix1", y = "column_sum_in_matrix2", color="species")
    fig.update_traces({'opacity': 1.0, 'marker':{'color':'lightskyblue','size':3.5}}, selector={'name': 'GRCh38'})
    fig.update_traces({'opacity': 1.0, 'marker':{'color':'mediumseagreen','size':3.5}}, selector={'name': 'mm10'})
    fig.update_traces({'opacity': 0.5, 'marker':{'color':'lightgray','size':3.5}}, selector={'name': 'Multiplet'})

    fig.update_layout(
        # width=600, height=500,
        plot_bgcolor='white', ###set the background color
        margin=dict(l=20, r=20, t=30, b=20), ###set the margin of the plot
        title=dict(text="Cell UMI Counts", font=dict(size=15),x=0.5,y=0.99),
        xaxis_title="mm10", yaxis_title="GRCh38")
    ##update the axis style
    fig.update_xaxes(
        mirror=True,
        ticks='outside',
        showline=True,
        linecolor='black',
        gridcolor='whitesmoke')
    fig.update_yaxes(
        mirror=True,
        ticks='outside',
        showline=True,
        linecolor='black',
        gridcolor='whitesmoke')
    return fig

# barcoderanks plot density
def segment_log_plot(y_data, x_start, x_end):
    log_max_x = np.log(len(y_data))
    log_max_y = np.log(max(y_data))
    segment_len = 0.0
    segment_idx = [x_start]
    for i in range(x_start, x_end):
        last_i = max(x_start, i-1)
        dx = (np.log(i) - np.log(last_i)) / log_max_x
        dy = (np.log(y_data[i]) - np.log(y_data[last_i])) / log_max_y
        segment_len += np.linalg.norm([dx, dy])
        if segment_len >= 0.02 and i > (segment_idx[-1] + 20):
            segment_idx.append(i+1)
            segment_len = 0.0
    if segment_idx[-1] != x_end:
        segment_idx.append(x_end)
    return segment_idx

# barcoderanks plot density
def plot_cmap(density):

    plot_colors =  [
        "#DDDDDD","#D6D9DC","#CFD6DB","#C8D3DA","#C1D0D9",
        "#BACDD9","#B3C9D8","#ACC6D7","#A5C3D6","#9EC0D6",
        "#97BDD5","#90BAD4","#89B6D3","#82B3D3","#7BB0D2",
        "#74ADD1","#6DAAD0","#66A6CF","#5FA3CF","#58A0CE",
        "#539DCC","#4F99CA","#4C95C8","#4992C6","#458EC3",
        "#428AC1","#3F87BF","#3B83BD","#3880BA","#347CB8",
        "#3178B6","#2E75B4","#2A71B1","#276DAF","#236AAD",
        "#2066AB","#1D62A8","#195FA6","#165BA4","#1358A2"]
    levels = len(plot_colors)
    ind = min(levels - 1, int(math.floor(levels * density)))
    return plot_colors[ind]

def plot_barcode_ranks(kneefile, _type="RNA"):

    dataframe_df = pd.read_csv(kneefile, header=0, sep = "\t")
    dataframe_df = dataframe_df.sort_values(by="UMI" , ascending=False)
    dataframe_df = dataframe_df.reset_index(drop=True)
    dataframe_df['New']=dataframe_df.index
    cell_bc = np.array(dataframe_df[dataframe_df['is_cell_barcode'] == 1].index)
    sorted_bc = np.array(dataframe_df.index)
    sorted_counts = np.array(dataframe_df['UMI'])
    total_bc = len(sorted_bc)
    ix1 = dataframe_df.drop_duplicates('is_cell_barcode',keep='first').index[1]-1
    ix2 = dataframe_df.drop_duplicates('is_cell_barcode',keep='last').index[0]
    plot_segments = []
    barcodeSegment = collections.namedtuple(
        'barcodeSegment', 
        ['start', 'end', 'density', 'legend'])

    plot_segments.append(barcodeSegment(start=0, end=ix1, density=1.0, legend=True))
    plot_segments.append(barcodeSegment(start=ix2+1, end=total_bc, density=0.0, legend=True))
    mixed_segments = segment_log_plot(sorted_counts, ix1, ix2)

    for i in range(len(mixed_segments) - 1):
        num_cells = sum([1 for i in range(mixed_segments[i], mixed_segments[i + 1]) if sorted_bc[i] in cell_bc])
        density = float(num_cells)/float(mixed_segments[i + 1]-mixed_segments[i])
        plot_segments.append(barcodeSegment(
                start=mixed_segments[i], end = mixed_segments[i + 1], density=density, legend=False))

    plot_data = []
    for plot_segment in plot_segments:
        start = max(0, plot_segment.start - 1)
        end = plot_segment.end
        selct_count = dataframe_df[start:end]
        dp_first = set(selct_count[selct_count[["UMI"]].duplicated(keep="first")].index)
        dp_last = set(selct_count[selct_count[["UMI"]].duplicated(keep="last")].index)
        dp_inter = dp_first & dp_last
        selct_count=selct_count.drop(list(dp_inter),axis=0)
        x = list(selct_count['New'])
        y = list(selct_count['UMI'])
        name = 'Cell' if plot_segment.density > 0 else 'Background'
        if plot_segment.density > 0:
            n_barcodes = plot_segment.end - plot_segment.start
            n_cells = int(round(plot_segment.density * n_barcodes))
            hover = "{:.0f}% Cell<br>({}/{})".format(100 * plot_segment.density, n_cells, n_barcodes)
        else:
            hover = "NOISE"

        data_dict = {
            "x": x,"y": y,"name": name, "hoverinfo": "text",
            "text": hover,"type": "scattergl","mode": "lines",
            "line": {"width": 3,"color": plot_cmap(plot_segment.density),},
            "showlegend": plot_segment.legend}
        plot_data.append(data_dict)

    plotly_data = [
        go.Scatter(
            x=dat['x'], y=dat['y'], name=dat['name'], mode=dat['mode'], 
            showlegend=dat['showlegend'],
            marker={'color': dat['line']['color']}, 
            line=dat['line'], text=dat['text']
            ) for dat in plot_data]
    
    layout = go.Layout(
        xaxis = dict(
            type="log", gridcolor="whitesmoke", title="Barcode in Rank-descending Order",
            color="black", showline=True, zeroline=True, linewidth=1, fixedrange= True,
            linecolor="black"),
        yaxis = dict(
            type="log", title=f"{_type} UMI counts", gridcolor="whitesmoke",
            linewidth=1, fixedrange= True, color="black", linecolor="black"),
        # height= 500, width= 500,
        plot_bgcolor='rgba(0,0,0,0)',hovermode='closest',paper_bgcolor='white',
        legend = dict(
            x=0.75,y=0.99,traceorder="normal",
            font = dict(
                family="Arial",size=12,color="black"),
            bordercolor="Black",borderwidth=0),
        margin = dict(l=0,r=0,b=0,t=0,pad=0.1),
        title=dict(text=f"{_type} UMI", font=dict(size=15),x=0.5,y=0.99),
        font = dict(size=10))
    
    fig = go.Figure(
        data=plotly_data, layout=layout)

    return fig


def saturation_plot(saturationfile):
    """序列饱和度绘制

    Args:
        saturationfile: 序列饱和度文件
    """
    df = pd.read_csv(saturationfile, sep = "\t")
    # 添加一行 0 值
    df = pd.concat(
        [pd.DataFrame([[0.0] * len(df.columns)], columns=df.columns), df],
        ignore_index=True,
    )
    fig1 = px.line(df, x="Mean Reads per Cell", y="Sequencing Saturation")
    fig1.update_traces(line_color="cornflowerblue", line_width=3)
    fig1.add_hline(y=0.9, line_width=1, line_dash="dash", line_color="black")
    fig1.update_layout(
        plot_bgcolor="white",
        margin=dict(l=20, r=20, t=30, b=20),
        showlegend=False,
    )
    # Update xaxis properties
    fig1.update_xaxes(
        mirror=True,
        ticks="outside",
        showline=True,
        linecolor="black",
        gridcolor="whitesmoke",
        title_text="Mean Reads per Cell",
    )
    fig1.update_yaxes(
        mirror=True,
        ticks="outside",
        showline=True,
        linecolor="black",
        gridcolor="whitesmoke",
        title_text="Sequencing Saturation",
        range=[0, 1],
    )

    fig2 = px.line(df, x="Mean Reads per Cell", y="Median Genes per Cell")
    fig2.update_traces(line_color="cornflowerblue", line_width=3)
    fig2.update_layout(
        plot_bgcolor="white",
        margin=dict(l=20, r=20, t=30, b=20),
        showlegend=False,
    )
    fig2.update_xaxes(
        mirror=True,
        ticks="outside",
        showline=True,
        linecolor="black",
        gridcolor="whitesmoke",
    )
    fig2.update_yaxes(
        mirror=True,
        ticks="outside",
        showline=True,
        linecolor="black",
        gridcolor="whitesmoke",
    )
    return fig1, fig2


my36colors = ['#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', 
                '#476D87', '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', 
                '#8C549C', '#585658', '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', 
                '#58A4C3', '#E4C755', '#F7F398', '#AA9A59', '#E63863', '#E39A35', 
                '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B', '#712820', '#DCC1DD', 
                '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963', '#968175']


def distribution_violin(clusterfile, samplename):

    clusterdf = pd.read_csv(clusterfile, header=0, sep="\t")
    cols = ["CB", "UMI counts", "Gene counts", "Mito Percentage", "xcoord", "ycoord", "UMAP1", "UMAP2", "Cluster"]
    clusterdf.columns = cols
    clusterdf["sample"] = samplename
    fig = make_subplots(rows=1, cols=3, 
                        subplot_titles=("Gene counts", "UMI counts", "Mito Percentage(%)"),
                        horizontal_spacing=0.1)
    fig.add_trace(go.Figure(data=go.Violin(y=clusterdf["Gene counts"], x=clusterdf["sample"], box_visible=True, line_color='black', meanline_visible=True, fillcolor='#E5D2DD', opacity=0.6, name='Gene counts')).data[0], row=1, col=1)
    fig.add_trace(go.Figure(data=go.Violin(y=clusterdf["UMI counts"], x=clusterdf["sample"], box_visible=True, line_color='black', meanline_visible=True, fillcolor='#53A85F', opacity=0.6, name='UMI counts')).data[0], row=1, col=2)
    fig.add_trace(go.Figure(data=go.Violin(y=clusterdf["Mito Percentage"], x=clusterdf["sample"], box_visible=True, line_color='black', meanline_visible=True, fillcolor='#F1BB72', opacity=0.6, name='Mito Percentage(%)')).data[0], row=1, col=3)
    
    fig.update_layout(
        # width=1200, height=500,
        plot_bgcolor='white', ###set the background color
        margin=dict(l=20, r=20, t=30, b=20))
    ##update the axis style
    fig.update_xaxes(
        mirror=True,
        ticks='outside',
        showline=True,
        linecolor='black',
        gridcolor='whitesmoke')
    
    fig.update_yaxes(
        mirror=True,
        ticks='outside',
        showline=True,
        linecolor='black',
        gridcolor='whitesmoke')    
    return fig

from bokeh.io import output_file, show
from bokeh.models import ColumnDataSource, Slider, HoverTool
from bokeh.plotting import figure
from bokeh.layouts import column
import pandas as pd
import numpy as np

def _umap_theme(fig):
    """
    白色背景
    """
    fig.update_layout(
        xaxis=dict(
            mirror=True,
            gridcolor="whitesmoke",
            color="black",
            showline=True,
            zeroline=True,
            linewidth=1,
            linecolor="black",
            zerolinecolor="whitesmoke",
        ),
        yaxis=dict(
            mirror=True,
            gridcolor="whitesmoke",
            linewidth=1,
            color="black",
            linecolor="black",
            zerolinecolor="whitesmoke",
        ),
        plot_bgcolor="rgba(0,0,0,0)",
        paper_bgcolor="white",
    )
    return fig


def check_layout(clusterdf, fig, _type="LD"):
    """
    调整绘图布局，根据类型设置坐标轴范围和刻度
    参数:
    clusterdf - 包含坐标数据的DataFrame
    fig - 要调整的图形对象
    _type - 类型，可选"LD"或"GM"或"UMAP"
    返回:
    调整后的图形对象
    """
    if _type == "LD":
        ###x, y 轴范围
        fig.update_layout(
            autosize=True,
            plot_bgcolor='white', ###set the background color
            xaxis_tickvals=[0, 2000,4000,6000,8000],
            yaxis_tickvals=[0, 2000,4000,6000,8000],
            yaxis_title=None, xaxis_title=None,
            xaxis=dict(
                    range=[0, 9399],  # 根据实际数据范围调整
                    constrain="domain"),  # 限制轴范围
            yaxis=dict(
                scaleanchor="x",  # y轴比例锚定到x轴
                scaleratio=1,      # 比例1:1
                constrain="domain",  # 限制轴范围
                range=[0, 9399]))
                
    elif _type == "GM":
        ymax = max(clusterdf['ycoord'])
        ymin = min(clusterdf['ycoord'])
        xmax = max(clusterdf['xcoord'])
        xmin = min(clusterdf['xcoord'])        
        # 计算最大范围并加10作为缓冲
        max_range = max(ymax-ymin, xmax-xmin) + 10
        center_x = (xmax + xmin) / 2
        center_y = (ymax + ymin) / 2
        # 设置轴范围以保持居中
        xaxis_range = [center_x - max_range/2 -100, center_x + max_range/2 + 100]
        yaxis_range = [center_y - max_range/2 -100, center_y + max_range/2 + 100]
        # 生成2k间隔的刻度值
        start_x = round(xaxis_range[0] / 2000) * 2000
        end_x = round(xaxis_range[1] / 2000) * 2000
        x_tickvals = list(range(int(start_x), int(end_x)+2000, 2000))
        start_y = round(yaxis_range[0] / 2000) * 2000
        end_y = round(yaxis_range[1] / 2000) * 2000
        y_tickvals = list(range(int(start_y), int(end_y)+2000, 2000))
        
        fig.update_layout(
            autosize=True,
            plot_bgcolor='white',
            xaxis_tickvals=x_tickvals,
            yaxis_tickvals=y_tickvals,
            yaxis_title=None, 
            xaxis_title=None,
            xaxis=dict(
                range=xaxis_range,
                constrain="domain"),
            yaxis=dict(
                scaleanchor="x",
                scaleratio=1,
                constrain="domain",
                range=yaxis_range))
    elif _type == "UMAP":
        ymax = max(clusterdf['UMAP2'])
        ymin = min(clusterdf['UMAP2'])
        xmax = max(clusterdf['UMAP1'])
        xmin = min(clusterdf['UMAP1'])        
        # 计算最大范围并加10作为缓冲
        max_range = max(ymax-ymin, xmax-xmin) + 2
        center_x = (xmax + xmin) / 2
        center_y = (ymax + ymin) / 2
        # 设置轴范围以保持居中
        xaxis_range = [center_x - max_range/2 - 2, center_x + max_range/2 + 2]
        yaxis_range = [center_y - max_range/2 - 2, center_y + max_range/2 + 2]
        # 生成5间隔的刻度值
        start_x = round(xaxis_range[0] / 5) * 5
        end_x = round(xaxis_range[1] / 5) * 5
        x_tickvals = list(range(int(start_x), int(end_x)+5, 5))
        start_y = round(yaxis_range[0] / 5) * 5
        end_y = round(yaxis_range[1] / 5) * 5
        y_tickvals = list(range(int(start_y), int(end_y)+5, 5))
        
        fig.update_layout(
            autosize=True,
            plot_bgcolor='white',
            xaxis_tickvals=x_tickvals,
            yaxis_tickvals=y_tickvals,
            yaxis_title=None, 
            xaxis_title=None,
            xaxis=dict(
                range=xaxis_range,
                constrain="domain"),
            yaxis=dict(
                scaleanchor="x",
                scaleratio=1,
                constrain="domain",
                range=yaxis_range))        

    else:
        print("Not an available oligochip! Only LD or GM can be chosen")  

    return fig

def _add_slider(fig, default_size=3.0):
    # 添加滑动条调整点的大小
    steps = [
        {
            "args": [{"marker.size": [round(size, 1)]}],
            "label": "" if i % 5 != 0 else f"{size:.1f}",  # 每5步显示一个标签
            "method": "restyle"
        }
        for i, size in enumerate(np.arange(1.0, 6.1, 0.1))  # 1.0到6.0，步长0.1
    ]
    
    fig.update_layout(
        autosize=True,
        margin={"l": 50, "r": 50, "t": 80, "b": 50},
        sliders=[
            {
                "active": int((default_size - 1.0) * 10),  # 计算默认位置
                "currentvalue": {
                    "prefix": "Point Size: ",
                    "font": {"size": 12},
                    "xanchor": "left",
                    "offset": 20
                },
                "steps": steps,
                "x": 0,
                "len": 0.5,
                "xanchor": "left",
                "yanchor": "top",
                "y": 0,
                "pad": {"t": 50, "b": 20},
                "transition": {"duration": 0},
            }
        ]
    )
    return fig

def _add_combined_sliders(fig, clusterdf, default_size=3.0, x_col="xcoord", y_col="ycoord"):
    import numpy as np

    # --- 获取原始 trace 信息 ---
    trace = fig.data[0]

    x_orig = np.array(clusterdf[x_col])
    y_orig = np.array(clusterdf[y_col])
    center_x = np.mean(x_orig)
    center_y = np.mean(y_orig)

    color = trace.marker.color if hasattr(trace.marker, "color") else None
    size = trace.marker.size if hasattr(trace.marker, "size") else default_size
    hovertemplate = trace.hovertemplate if hasattr(trace, "hovertemplate") else None
    customdata = trace.customdata if hasattr(trace, "customdata") else None
    mode = trace.mode if hasattr(trace, "mode") else "markers"

    # --- 点大小滑动条 ---
    size_steps = [
        {
            "args": [{"marker.size": [round(s, 1)]}],
            "label": "" if i % 5 != 0 else f"{s:.1f}",
            "method": "restyle"
        }
        for i, s in enumerate(np.arange(1.0, 6.1, 0.1))
    ]

    # --- 旋转滑动条 ---
    rotation_steps = []
    for angle in range(0, 360, 15):
        theta = np.radians(angle)
        x_rot = (x_orig - center_x) * np.cos(theta) - (y_orig - center_y) * np.sin(theta) + center_x
        y_rot = (x_orig - center_x) * np.sin(theta) + (y_orig - center_y) * np.cos(theta) + center_y

        update_args = {
            "x": [x_rot.tolist()],
            "y": [y_rot.tolist()],
            "mode": [mode],
            "marker": {
                "color": color,
                "size": size,
            },
        }

        if hovertemplate is not None:
            update_args["hovertemplate"] = [hovertemplate]
        if customdata is not None:
            update_args["customdata"] = [customdata]

        rotation_steps.append({
            "args": [update_args],
            "label": f"{angle}°",
            "method": "update"
        })

    # --- 3. 合并布局 ---
    fig.update_layout(
        autosize=True,
        margin={"l": 50, "r": 50, "t": 80, "b": 150},  # 底部留更多空间
        sliders=[
            # 点大小滑动条
            {
                "active": int((default_size - 1.0) * 10),
                "currentvalue": {"prefix": "Point Size: ", "xanchor": "left", "offset": 20},
                "steps": size_steps,
                "x": 0,
                "len": 0.5,
                "xanchor": "left",
                "yanchor": "top",
                "y": 0,
                "pad": {"t": 50, "b": 20},
                "transition": {"duration": 0},
            },
            # 旋转滑动条
            {
                "active": 0,
                "currentvalue": {"prefix": "Rotation: ", "suffix": "°", "xanchor": "left", "offset": 20},
                "steps": rotation_steps,
                "x": 0.55,
                "len": 0.5,
                "xanchor": "left",
                "yanchor": "top",
                "y": 0,
                "pad": {"t": 50, "b": 10},
                "transition": {"duration": 50}
            }
        ]
    )
    return fig

# def spatial_scatter(_infile, _type = "UMI", oligochip="LD"):

#     clusterdf = pd.read_csv(_infile, header=0, sep="\t")
#     clusterdf["log_nUMI"] = np.log(clusterdf["nCount_Spatial"]+1)
#     clusterdf = clusterdf.sort_values(by="Cluster")
#     clusterdf["Pct"] = clusterdf["Cluster"].map(clusterdf["Cluster"].value_counts(normalize=True))
#     clusterdf["Cluster"] = clusterdf["Cluster"].astype("category")
#     ##sort the clusterdf by the Cluster column
#     clusterdf = clusterdf.sort_values(by="Cluster")
#     color_map = {}
#     for i in range(len(clusterdf["Cluster"].unique())):
#         color_map[clusterdf["Cluster"].unique()[i]] = my36colors[i]

#     if _type == "Cluster":
#         fig1 = _umap_theme(px.scatter(clusterdf, 
#                                       x="xcoord", 
#                                       y="ycoord", 
#                                       color="Cluster",
#                                       hover_data=["Pct"],
#                                       color_discrete_sequence=color_map,
#                                       ))    
#         fig1.update_layout(
#             title=dict(text="Cells Colored by Cluster", font=dict(size=15), x=0.5, y=0.95),
#             legend=dict(
#                 title="",
#                 borderwidth=0,
#             ),
#         )
#         fig1.update_traces(marker_size=3)
#         fig1 = check_layout(clusterdf, fig1, oligochip)
#         # fig1 = _add_combined_sliders(fig1, clusterdf, color_col="Cluster")  # 传递颜色列名
#         fig1 = _add_combined_sliders(fig1, clusterdf)
        
#         fig2 = _umap_theme(
#             px.scatter(
#                 clusterdf,
#                 x="UMAP1",
#                 y="UMAP2",
#                 color="Cluster",
#                 hover_data=["Pct"],
#                 color_discrete_sequence=color_map,
#             )
#         )
#         fig2.update_layout(
#             title=dict(text="UMAP Projection of Cells Colored by Cluster", font=dict(size=15), x=0.5, y=0.95),
#             legend=dict(
#                 title="",
#                 borderwidth=0,
#             ),
#         )
#         fig2.update_traces(marker_size=3)
#         fig2 = _add_slider(fig2, default_size=3)  # 添加滑动条
#         fig2 = check_layout(clusterdf, fig2, "UMAP")
#     elif _type == "UMI":
#         fig1 = _umap_theme(px.scatter(clusterdf, 
#                                       x="xcoord", 
#                                       y="ycoord", 
#                                       color="log_nUMI",
#                                       ))    
#         fig1.update_layout(
#             title=dict(text="Cells Colored by log(UMI counts)", font=dict(size=15), x=0.5, y=0.95),
#             legend=dict(
#                 title="",
#                 borderwidth=0,
#             ),
#         )
#         fig1.update_traces(marker_size=3)
#         fig1 = check_layout(clusterdf, fig1, oligochip)
#         # fig1 = _add_combined_sliders(fig1, clusterdf, default_size=3)
#         fig1 = _add_combined_sliders(fig1, clusterdf)
        
#         fig2 = _umap_theme(
#             px.scatter(
#                 clusterdf,
#                 x="UMAP1",
#                 y="UMAP2",
#                 color="log_nUMI",
#             )
#         )

#         fig2.update_layout(
#             title=dict(text="UMAP Projection of Cells Colored by log(UMI counts)", font=dict(size=15), x=0.5, y=0.95),
#             legend=dict(
#                 title="",
#                 borderwidth=0,
#             ),
#         )
#         fig2.update_traces(marker_size=3)
#         fig2 = _add_slider(fig2, default_size=3)  # 添加滑动条
#         fig2 = check_layout(clusterdf, fig2, "UMAP")
#     else:
#         print("Please input the correct type: UMI or Cluster")

#     return fig1, fig2

def plot_sb_cb_umi_knee(_sb_umi_file, _cb_umi_file):

    sb_umi_df = pd.read_csv(_sb_umi_file, sep = "\t", header = 0)
    # 创建 Log-Log 散点图
    plot1 = px.line(sb_umi_df, x='rank', y='umi_count_sum', log_x=True, log_y=True, markers=True,
                    title='Log-Log Plot of UMI Count Sum by SB',
                    labels={'rank': 'Rank (log scale)', 'umi_count_sum': 'Sum of UMI Count (log scale)'})    
    plot1.update_layout(
        width=600, height=500,
        plot_bgcolor='white', ###set the background color
        yaxis_title=None, xaxis_title=None,
        margin=dict(l=20, r=20, t=30, b=20))
    ##update the axis style
    plot1.update_xaxes(mirror=True,ticks='outside',showline=True,linecolor='black', gridcolor='whitesmoke')
    plot1.update_yaxes(mirror=True,ticks='outside',showline=True,linecolor='black', gridcolor='whitesmoke')

    cb_umi_df = pd.read_csv(_cb_umi_file, sep = "\t", header = 0)
    # 创建 Log-Log 散点图
    plot2 = px.line(cb_umi_df, x='rank', y='umi_count_mean', color='cluster', log_x=True, log_y=True, markers=True,
                    # color_discrete_map = {"blue": "single-cluster", "yellow": "multi-cluster", "red": "no-cluster"},
                    color_discrete_sequence = px.colors.qualitative.Pastel2,
                    title='Log-Log Plot of UMI Count Mean by CB',
                    labels={'rank': 'Rank (log scale)', 'umi_count_mean': 'Mean of UMI Count (log scale)'})    
    plot2.update_layout(
        width=600, height=500,
        plot_bgcolor='white', ###set the background color
        yaxis_title=None, xaxis_title=None,
        margin=dict(l=20, r=20, t=30, b=20))
    ##update the axis style
    plot2.update_xaxes(mirror=True,ticks='outside',showline=True,linecolor='black', gridcolor='whitesmoke')
    plot2.update_yaxes(mirror=True,ticks='outside',showline=True,linecolor='black', gridcolor='whitesmoke')
    
    return plot1, plot2


from bokeh.plotting import figure
from bokeh.layouts import row, column  # 确保导入row和column
from bokeh.models import (ColumnDataSource, Slider, HoverTool, 
                         CustomJS, LinearColorMapper, ColorBar,Legend)
from bokeh.palettes import Viridis256
import numpy as np
import pandas as pd
from bokeh.embed import file_html, components
from bokeh.resources import CDN
class BokehFigure:
    def __init__(self):
        self.plot = None
        self.layout = None
        
    def bokeh_to_html(self):
        """返回可嵌入的HTML字符串"""
        return file_html(self.layout, CDN, "Spatial Scatter")
    
    def to_components(self):
        """返回(script, div)元组"""
        return components(self.layout)



def prepare_legend(clusterdf, color_map, source, fig, xaxis, yaxis):
    # 创建分类图例项
    legend_items = []
    for cluster in sorted(clusterdf["Cluster"].unique()):
        cluster = str(cluster)  # 确保cluster是字符串
        subset = ColumnDataSource({
            xaxis: source.data[xaxis][np.array(source.data['cluster']) == cluster],
            yaxis: source.data[yaxis][np.array(source.data['cluster']) == cluster]
        })
        r = fig.scatter(
            x=xaxis, y=yaxis, 
            size=3,
            fill_color=color_map[cluster],
            line_color=None,
            fill_alpha=0.7,
            source=subset,
            name=str(cluster)  # 确保名称是字符串
        )
        legend_items.append((str(cluster), [r]))  # 确保图例标签是字符串
    
    # 添加外部图例
    legend = Legend(
        items=legend_items,
        location="center_right",
        title="Clusters",
        label_text_font_size="10px",
        glyph_height=15,
        glyph_width=15,
        spacing=5,
        label_standoff=8,
        margin=10
    )
    return legend, legend_items

def prepare_axis_range(clusterdf, _type="LD"):

    if _type == "LD":
        xaxis_range = (0, 9399)
        yaxis_range = (0,9399)
    elif _type == "GM":
        ymax = max(clusterdf['ycoord'])
        ymin = min(clusterdf['ycoord'])
        xmax = max(clusterdf['xcoord'])
        xmin = min(clusterdf['xcoord'])        
        # 计算最大范围并加10作为缓冲
        max_range = max(ymax-ymin, xmax-xmin) + 10
        center_x = (xmax + xmin) / 2
        center_y = (ymax + ymin) / 2
        # 设置轴范围以保持居中
        xaxis_range = (center_x - max_range/2 -100, center_x + max_range/2 + 100)
        yaxis_range = [center_y - max_range/2 -100, center_y + max_range/2 + 100]
    elif _type == "UMAP":
        ymax = max(clusterdf['UMAP2'])
        ymin = min(clusterdf['UMAP2'])
        xmax = max(clusterdf['UMAP1'])
        xmin = min(clusterdf['UMAP1'])        
        # 计算最大范围并加10作为缓冲
        max_range = max(ymax-ymin, xmax-xmin) + 2
        center_x = (xmax + xmin) / 2
        center_y = (ymax + ymin) / 2
        # 设置轴范围以保持居中
        xaxis_range = [center_x - max_range/2 - 2, center_x + max_range/2 + 2]
        yaxis_range = [center_y - max_range/2 - 2, center_y + max_range/2 + 2]
    else:
        print("Not an available oligochip! Only LD or GM can be chosen")  

    return xaxis_range, yaxis_range

# 更新 CLUSTER_ROTATE_CODE
CLUSTER_ROTATE_CODE = """
        const data = source.data;
        const ox = data['original_x'];
        const oy = data['original_y'];
        const center_x = ox.reduce((a,b) => a+b, 0) / ox.length;
        const center_y = oy.reduce((a,b) => a+b, 0) / oy.length;
        
        const angle = rotate_slider.value * Math.PI / 180;
        const cosA = Math.cos(angle);
        const sinA = Math.sin(angle);
        
        // 获取翻转参数
        const flip_horizontal = flip_h_slider.value;
        const flip_vertical = flip_v_slider.value;
        
        // 更新所有渲染器
        for (const r of renderers) {
            const source = r.data_source;
            const cluster = r.name;
            const x = source.data['x'];
            const y = source.data['y'];
            const orig_x = data['original_x'].filter((_, idx) => data['cluster'][idx] === cluster);
            const orig_y = data['original_y'].filter((_, idx) => data['cluster'][idx] === cluster);
            
            for (let i = 0; i < x.length; i++) {
                let dx = orig_x[i] - center_x;
                let dy = orig_y[i] - center_y;
                
                // 应用翻转
                if (flip_horizontal) dx = -dx;
                if (flip_vertical) dy = -dy;
                
                // 应用旋转
                x[i] = dx * cosA - dy * sinA + center_x;
                y[i] = dx * sinA + dy * cosA + center_y;
            }
            
            r.glyph.size = size_slider.value;
        }
        source.change.emit();
        """

# 更新 UMI_RORATE_CODE
UMI_RORATE_CODE = """
    const data = source.data;
    const ox = data['original_x'];
    const oy = data['original_y'];
    const center_x = ox.reduce((a,b) => a+b, 0) / ox.length;
    const center_y = oy.reduce((a,b) => a+b, 0) / oy.length;

    const angle = rotate_slider.value * Math.PI / 180;
    const cosA = Math.cos(angle);
    const sinA = Math.sin(angle);
    
    // 获取翻转参数
    const flip_horizontal = flip_h_slider.value;
    const flip_vertical = flip_v_slider.value;

    for (let i = 0; i < data['x'].length; i++) {
        let dx = ox[i] - center_x;
        let dy = oy[i] - center_y;
        
        // 应用翻转
        if (flip_horizontal) dx = -dx;
        if (flip_vertical) dy = -dy;
        
        // 应用旋转
        data['x'][i] = dx * cosA - dy * sinA + center_x;
        data['y'][i] = dx * sinA + dy * cosA + center_y;
    }

    renderer.glyph.size = size_slider.value;
    source.change.emit();
    """

def spatial_scatter(_infile, _type = "UMI", oligochip="LD"):
    # 准备数据源
    clusterdf = pd.read_csv(_infile, header=0, sep="\t")
    clusterdf["log_nUMI"] = np.log(clusterdf["nCount_Spatial"]+1)
    clusterdf = clusterdf.sort_values(by="Cluster")
    clusterdf["Pct"] = clusterdf["Cluster"].map(clusterdf["Cluster"].value_counts(normalize=True))
    clusterdf["Cluster"] = clusterdf["Cluster"].astype("category")
    ##sort the clusterdf by the Cluster column
    clusterdf = clusterdf.sort_values(by="Cluster")
    color_map = {str(cluster): my36colors[i % len(my36colors)] 
                 for i, cluster in enumerate(clusterdf["Cluster"].unique())}

    spatial_source = ColumnDataSource(data=dict(
        x=clusterdf["xcoord"],
        y=clusterdf["ycoord"],
        cluster=clusterdf["Cluster"].astype(str),  # 确保cluster是字符串
        pct=clusterdf["Pct"],
        umap1=clusterdf["UMAP1"],
        umap2=clusterdf["UMAP2"],
        log_numi=clusterdf["log_nUMI"],
        original_x=clusterdf["xcoord"].copy(),
        original_y=clusterdf["ycoord"].copy(),
        color=[color_map[str(c)] for c in clusterdf["Cluster"]]
    ))
    
    l_xaxis_range, l_yaxis_range = prepare_axis_range(clusterdf, _type=oligochip)
    r_xaxis_range, r_yaxis_range = prepare_axis_range(clusterdf, _type="UMAP")

    if _type == "Cluster":
    
        # ========== 左图，Cluster空间分布，点调节大小+图案旋转 ==========
        spatial_fig = figure(
            width=450, height=400,
            x_range=l_xaxis_range, 
            y_range=l_yaxis_range,
            x_axis_label="X",
            y_axis_label="Y",
            title="Spatial Coordinates (Colored by Cluster)",
            tools="pan,wheel_zoom,box_zoom,reset",
            toolbar_location="above"
        )
        
        legend, legend_items = prepare_legend(clusterdf, color_map, spatial_source, spatial_fig, 'x', 'y')
        spatial_fig.add_layout(legend, 'right')
        
        # 左图控制滑块
        size_slider = Slider(width = 150, height=25, start=1, end=7, value=3, step=0.5, title="Point Size")
        rotate_slider = Slider(width = 150, height=25, start=0, end=360, value=0, step=1, title="Rotation Angle")
        # 新增翻转滑块
        flip_h_slider = Slider(width = 150, height=25, start=0, end=1, value=0, step=1, title="Flip Horizontal")
        flip_v_slider = Slider(width = 150, height=25, start=0, end=1, value=0, step=1, title="Flip Vertical")

        # 左图回调
        spatial_callback = CustomJS(args={
            'source': spatial_source,
            'size_slider': size_slider,
            'rotate_slider': rotate_slider,
            'flip_h_slider': flip_h_slider,
            'flip_v_slider': flip_v_slider,
            'renderers': [r for (_, [r]) in legend_items]  # 获取所有渲染器
        }, code=CLUSTER_ROTATE_CODE)

        size_slider.js_on_change('value', spatial_callback)
        rotate_slider.js_on_change('value', spatial_callback)
        flip_h_slider.js_on_change('value', spatial_callback)
        flip_v_slider.js_on_change('value', spatial_callback)

        # ========== 右图：UMAP坐标 + cluster颜色+点大小调节 ==========
        umap_fig = figure(
            width=450, height=400,
            x_range=r_xaxis_range, 
            y_range=r_yaxis_range,
            x_axis_label="UMAP1",
            y_axis_label="UMAP2",
            title="UMAP Projection (Colored by Cluster)",
            tools="pan,wheel_zoom,box_zoom,reset"
        )
        legend, legend_items = prepare_legend(clusterdf, color_map, spatial_source, umap_fig, 'umap1', 'umap2')
        umap_fig.add_layout(legend, 'right')
        
        # 右图控制滑块
        umap_size_slider = Slider(width = 150, height=25, start=1, end=7, value=3, step=0.5, title="Point Size")
        # 右图回调
        umap_callback = CustomJS(args={
            'legend_items': legend_items,  # 直接传递图例项
            'slider': umap_size_slider
        }, code="""
        for (const item of legend_items) {
            const [_, renderers] = item;
            for (const r of renderers) {
                r.glyph.size = slider.value;
            }
        }
        """)
        umap_size_slider.js_on_change('value', umap_callback)
        
        # ========== 组合图形 ==========
        spatial_controls = column(row(size_slider, rotate_slider, spacing = 25), 
                                  row(flip_h_slider, flip_v_slider, spacing = 25))
        umap_controls = column(umap_size_slider)
        
        fig = BokehFigure()
        fig.layout = row(
            column(spatial_fig, spatial_controls),
            column(umap_fig, umap_controls),
            sizing_mode='fixed',
            spacing=50,  # 两图之间添加50像素间距
        )
    elif _type == "UMI":
        # ========== 右图：空间分布坐标 + log_nUMI颜色 ==========
        umi_mapper = LinearColorMapper(
            palette=Viridis256,
            low=clusterdf["log_nUMI"].min(),
            high=clusterdf["log_nUMI"].max()
        )
        spatial_fig = figure(
            width=450, height=400,
            x_range=l_xaxis_range, 
            y_range=l_yaxis_range,
            title="Spatial Coordinates (Colored by log UMI)",
            tools="pan,wheel_zoom,box_zoom,reset",
            toolbar_location="above"
        )
        # 绘制所有点（统一渲染器）
        spatial_renderer = spatial_fig.scatter(
            x='x', y='y', 
            size=3,
            fill_color={'field': 'log_numi', 'transform': umi_mapper},
            line_color=None,
            fill_alpha=0.7,
            source=spatial_source
        )
        # 左图控制滑块
        size_slider = Slider(width = 150, height=25, start=1, end=7, value=3, step=0.5, title="Point Size")
        rotate_slider = Slider(width = 150, height=25, start=0, end=360, value=0, step=1, title="Rotation Angle")
        # 新增翻转滑块
        flip_h_slider = Slider(width = 150, height=25, start=0, end=1, value=0, step=1, title="Flip Horizontal")
        flip_v_slider = Slider(width = 150, height=25, start=0, end=1, value=0, step=1, title="Flip Vertical")
        
        # 左图回调（旋转+大小）
        spatial_callback = CustomJS(args={
            'source': spatial_source,
            'size_slider': size_slider,
            'rotate_slider': rotate_slider,
            'flip_h_slider': flip_h_slider,
            'flip_v_slider': flip_v_slider,
            'renderer': spatial_renderer
        }, code=UMI_RORATE_CODE)

        size_slider.js_on_change('value', spatial_callback)
        rotate_slider.js_on_change('value', spatial_callback)
        flip_h_slider.js_on_change('value', spatial_callback)
        flip_v_slider.js_on_change('value', spatial_callback)
        # 添加颜色条
        color_bar = ColorBar(
            color_mapper=umi_mapper,
            label_standoff=12,
            location=(0,0),
            title="log(UMI)"
        )
        spatial_fig.add_layout(color_bar, 'right')        
        
        # ========== 右图：UMAP坐标 + log_nUMI颜色 ==========      
        umap_fig = figure(
            width=450, height=400,
            x_range=r_xaxis_range, 
            y_range=r_yaxis_range,
            title="UMAP Projection (Colored by log_nUMI)",
            tools="pan,wheel_zoom,box_zoom,reset",
            toolbar_location="above"
        )
        
        umap_renderer = umap_fig.scatter(
            x='umap1', y='umap2', size=3,
            fill_color={'field': 'log_numi', 'transform': umi_mapper},
            line_color=None,
            fill_alpha=0.7,
            source=spatial_source
        )
        
        umap_fig.add_layout(color_bar, 'right')
        
        # 右图控制滑块
        umap_size_slider = Slider(width = 150, height=25, start=3, end=20, value=8, step=0.5, title="Point Size")
        umap_callback = CustomJS(args={'renderer': umap_renderer, 'slider': umap_size_slider}, 
                            code="renderer.glyph.size = slider.value;")
        umap_size_slider.js_on_change('value', umap_callback)

        # ========== 组合图形 ==========
        spatial_controls = column(row(size_slider, rotate_slider, spacing = 25), 
                                  row(flip_h_slider, flip_v_slider, spacing = 25))
        umap_controls = column(umap_size_slider)
        
        fig = BokehFigure()
        fig.layout = row(
            column(spatial_fig, spatial_controls),
            column(umap_fig, umap_controls),
            sizing_mode='fixed',
            spacing=50,  # 两图之间添加50像素间距
        )
    else:
        print("Please input the correct type: UMI or Cluster")
    
    return fig
