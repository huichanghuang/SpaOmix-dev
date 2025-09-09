import math
import pandas as pd
import numpy as np
from scipy.spatial.distance import cdist
import matplotlib.pyplot as plt
import argparse

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--infile", type=str, required=True)
    parser.add_argument("--outfile", type=str, required=True)
    return parser.parse_args()

def plot(data, outfile):
    nearest_distances = data['nearest_distance']
    bin_width = np.std(nearest_distances)/2  # 改用标准差控制分箱
    bin_count = int((np.max(nearest_distances) - np.min(nearest_distances)) / bin_width)
    bin_count = max(10, min(bin_count, 30))  # 限制在10-30个bin之间

    plt.figure(figsize=(15, 6))
    n, bins, patches = plt.hist(nearest_distances, 
                            bins=bin_count,
                            color='skyblue',
                            edgecolor='black',
                            alpha=0.7)

    # 优化标签显示
    for i in range(len(patches)):
        if n[i] > 0:  # 只显示有数据的柱子
            plt.text(patches[i].get_x() + patches[i].get_width()/2, 
                    patches[i].get_height()+0.3,
                    f"{int(n[i])}",
                    ha='center', 
                    fontsize=8,
                    bbox=dict(facecolor='white', alpha=0.7, edgecolor='none'))

    # 添加统计信息
    mean_val = np.mean(nearest_distances)
    median_val = np.median(nearest_distances)
    plt.axvline(mean_val, color='red', linestyle='--', label=f'Mean: {mean_val:.1f} µm')
    plt.axvline(median_val, color='green', linestyle=':', label=f'Median: {median_val:.1f} µm')
    # 添加图例（现在会显示统计线的标签）
    plt.legend()

    # 添加分布曲线
    from scipy.stats import gaussian_kde
    kde = gaussian_kde(nearest_distances)
    x = np.linspace(np.min(nearest_distances), np.max(nearest_distances), 1000)
    plt.plot(x, kde(x)*len(nearest_distances)*bin_width, 'b-', linewidth=1.5)

    # 设置x轴刻度（显示取整后的刻度，但保持实际bin边界不变）
    binwidth = bins[1]- bins[0]
    rounded_up = math.ceil(binwidth / 10) * 10
    rounded_bins = np.arange(0, bin_count * rounded_up, rounded_up) 
    plt.xticks(rounded_bins, rotation=45, ha='right')  # 显示取整后的刻度
    plt.title('Nearest Neighbor Distance Distribution', fontsize=14)
    plt.xlabel('Distance(µm)', fontsize=12)
    plt.ylabel('Number of CBs', fontsize=12)
    plt.grid(axis='y', linestyle='--', alpha=0.5)
    plt.tight_layout()

    plt.savefig(outfile)
    return


def main():
    args = parse_args()
    infile = args.infile
    outfile = args.outfile

    data = pd.read_csv(infile, sep="\t", header=0)
    # 提取坐标
    points = data[['xcoord', 'ycoord']].values
    # 计算距离矩阵
    distance_matrix = cdist(points, points, 'euclidean')
    # 将对角线设为无穷大（排除自身）
    np.fill_diagonal(distance_matrix, np.inf)
    # 找到每个点的最近邻索引和距离
    nearest_indices = np.argmin(distance_matrix, axis=1)
    nearest_distances = distance_matrix[np.arange(len(nearest_indices)), nearest_indices]
    # 添加结果到 DataFrame
    data['nearest_cb'] = data['cb'].iloc[nearest_indices].values
    data['nearest_distance'] = nearest_distances.round(2)

    plot(data, outfile)

    return

if __name__ == "__main__":
    main()

