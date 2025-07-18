# 脚本目录说明

## 标准分析

标准分析使用脚本 [standard_analysis.R](standard_analysis.R)。

```bash
# 查看帮助文档
Rscript standard_analysis.R -h
```

- 参数 `--input` 需要样本配置文件，文件内容示例如下：

```bash
sample,path,group,age
pbmc3k,/path/to/filtered,a,18
pbmc3k_b,/path/to/filtered,b,18
```

其中 `sample`,`path` 列是必须的，path 对应的路径支持各种单细胞输入格式，可以是 cellranger 的输出、STARsolo的输出、dnbc4tools的输出、h5ad、rds、txt、csv、tsv等。
其余的列是任意的，例如你可以随意添加性别、组织、年龄等列。

运行：

```bash
Rscript standard_analysis.R --input 配置文件
```

- 参数 `--group` 为可选的，是运行差异分析的分组文件，文件示例如下：

```bash
name,colname,group1,group2
a-vs-b,group,a,b
1-vs-2,seurat_clusters,1,2
```

如果你想要例如类别1，2，3与4，5，6比较，可以以冒号（;）进行分割，示例如下：

```bash
name,colname,group1,group2
123-vs-456,seurat_clusters,1;2;3,4;5;6
```

文件四列都是必须的，其中第一列为组别名称，第二列为元数据的列，三、四列为比较组。

## 生成标准分析报告

在 `standard_analysis.R` 的输出目录下创建项目文件 `project.yaml`，项目文件的示例内容如下：

```yaml
项目类型: 单细胞 3’转录组测序
合同编号: 新的合同编号
项目名称: 新的项目名称
客户单位: 利德健康客户
报告时间: null
样本信息:
  样本名: [sample1,sample2,sample3,sample4]
  分组: [group1,group1,group2,group2]
  物种: [hsa,hsa,hsa,hsa]
  组织: [brain,brain,brain,brain]
```

运行下面的命令，其中 `$outdir` 指代 `standard_analysis.R` 脚本的输出目录

### html 报告

```bash
cp -rf report $outdir && cd $outdir/report
quarto render --to aria-html
```

### pdf 报告

#### 模板 1
此模板生成的格式为 html，但可以通过浏览器打印转为 pdf
```bash
quarto render --to paged-html -P paged:true
```

html 转 pdf:

```bash
python pprint -i report-paged.html -o report.pdf
```
> 备注：虽然可以命令行转，但是目前来说centos 7 不好配置环境，所以暂时推荐自己使用谷歌浏览器的打印功能转。

### word 报告

```bash
quarto render --to aria-docx
python script/add_docx.py -i report.docx -o report.docx
```


### 删除报告源码

交付的话，最好把报告源码删除

```bash
cp report.html ../ && cd ../ && rm -rf report
```
