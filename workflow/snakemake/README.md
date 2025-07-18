## snakemake 流程使用文档

### 使用

1. 更新配置文件；

默认的参数配置可以在 `config/config.yaml` 文件中找到，可以通过 snakemake `--configfile` 或者`--config` 参数对默认参数进行修改。

流程运行时必须在配置文件中包含 `samples` 键值以及 `outdir`，如 [examples/input.yaml](examples/input.yaml) 所示。

2. 运行

- 本地运行

```bash
snakemake --snakefile Snakefile \
    --configfile input.yaml \
    --cores 80
```

- 投递任务运行

运行前需安装好 snakemake 插件：

```bash
pip install snakemake-executor-plugin-cluster-generic
```


```bash
snakemake --snakefile Snakefile \
    --configfile input.yaml \
    --executor cluster-generic --cluster-generic-submit-cmd "qsub -V -l walltime=1000:00:00" \
    --jobs 5
```


### 配置文件参数说明

1. 最简单的情况，你可以参考下面的配置定义好每个样本的输入文件以及输出目录。

```yaml
samples:
  sampleA: 
    R1: 
      - "/path/to/test1_L001_R1.fastq.gz"
      - "/path/to/test1_L002_R1.fastq.gz"
    R2: 
      - "/path/to/test1_L001_R2.fastq.gz"
      - "/path/to/test1_L002_R2.fastq.gz
  sampleB:
    R1: "/path/to/test2_L001_R1.fastq.gz"
    R2: "/path/to/test2_L001_R1.fastq.gz"

outdir: "/path/to/outdir"
```

2. 修改物种参数，全局默认值为 human，可以修改为 mouse

```yaml
species: mouse
```

同时你可以为每个样本设置物种参数，如下所示：

```yaml
samples:
  sampleA:
    species: mouse
  sampleB:
    species: human
```

3. dev 参数

有的时候，需要区分显示内外部指标信息，你可以通过设置 dev 参数，dev为 fasle 时，则不显示全部指标。

```yaml
dev: false
```

同时也可以为每个样本分别设置 dev 参数：

```yaml
samples:
  sampleA:
    dev: false
  sampleB:
    dev: false
```


4. recall 重新识别细胞

有的时候，可能个别样本自动识别的细胞不符合预期，可以增加参数重新识别细胞：

> 注意：运行此参数，你需要将 `snakemake` 流程输出目录生成的标记文件删除：
> `rm 02.count/*/.count.done`

```yaml
samples:
  sampleA:
    recall: true
    min_umi: 1000
  sampleB:
    recall: true
    force_cell: 5000
```

