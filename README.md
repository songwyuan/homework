# GEO-RNAseq-Workflow
GEO转录组数据分析流程的实操作业

# 📚 RNA-seq 数据分析学习任务

## 任务目标

1. 熟悉RNA-seq数据从下载到分析的全流程，掌握Linux、R、组学分析基础技能。

2. 具备独立处理转录组数据的能力，并能对结果进行生物学解释。

3. 学会用Markdown和GitHub记录和展示分析流程，提升项目管理和展示技能。

## 任务背景

本项目基于指定的GEO编号`GSE255647`，通过RNA-seq数据的下载、预处理、比对、定量、差异表达分析和功能富集分析，完成转录组数据的全流程分析，并以Markdown形式记录过程，最终以GitHub仓库展示。

**任务要求**：分析SARS-CoV-2以1 MOI感染Calu-3/2B4细胞系12h后，转录组水平的变化。

## 📁 任务模块

### 模块1：任务准备

**目标**：基于`micromamba`配置环境，熟悉工具和资源。（包括新建R4.4.1的环境，并安装以下工具）

**任务**：

安装并熟悉以下工具：

1. fastqc、multiqc、fastp：质控与清洗

2. hisat2、stringtie：比对与定量

3. samtools：BAM文件操作

4. R与DESeq2、clusterProfiler、ggplot2

5. 注册GitHub账号，掌握基本操作

**思考题**：

1. 如何在没有管理员权限下安装fastqc？

2. GitHub的README文件有什么作用？

### 模块2：数据下载与预处理

**目标**：从GEO下载原始数据，进行质控和清洗。

**任务**：

1. 使用wget、curl或aspera下载指定GEO编号的FASTQ数据
   #下载脚本
#!/bin/bash
set -euo pipefail
export PATH=/mnt/alamo01/yuansongwei7@mgt01:/mnt/alamo01/users/yuansongwei7:$PATH

OUTDIR="/mnt/alamo01/users/yuansongwei7/download_data/GSE255647"
mkdir -p "$OUTDIR"

# SRR▒~V▒~O▒▒~H~W表▒~H▒~O▒▒| ▒▒~M▒▒~\~@▒~A▒~G▒▒~L修▒~T▒▒~I
SRR_LIST=(SRR27961772 SRR27961773 SRR27961774 SRR27961775 SRR27961776
SRR27961777 SRR27961778 SRR27961779 SRR27961780 SRR27961781
SRR27961782 SRR27961783 SRR27961784 SRR27961785 SRR27961786
SRR27961787 SRR27961788 SRR27961789 SRR27961790 SRR27961791
SRR27961792 SRR27961793 SRR27961794 SRR27961795 SRR27961796
SRR27961797 SRR27961798)

# ▒~K载 + 解▒~N~K + ▒~N~K缩▒~G▒▒~U▒
download_and_convert() {
    local srr="$1"
    local url="https://sra-pub-run-odp.s3.amazonaws.com/sra/${srr}/${srr}"
    local sra_file="${OUTDIR}/${srr}.sra"
    local done_flag="${OUTDIR}/${srr}.done"

    # 跳▒~G已▒~L▒~H~P▒~Z~D
    if [[ -f "$done_flag" ]]; then
        echo "▒~\~E $srr already processed. Skipping."
        return
    fi

    echo "▒~_~S▒ Downloading $srr..."
    wget -O "$sra_file" "$url" || {
        echo "▒~]~L Failed to download $srr"
        return 1
    }

    echo "▒~_~T~D Converting $srr to FASTQ..."
    fastq-dump --split-files --outdir "$OUTDIR" "$sra_file" || {
        echo "▒~]~L fastq-dump failed for $srr"
        rm -f "$sra_file"
        return 1
    }

     echo "▒~V~R~_~S▒~V~R Compressing FASTQ files for $srr..."
    if [[ -f "${OUTDIR}/${srr}_1.fastq" ]]; then
        gzip "${OUTDIR}/${srr}_1.fastq"
    fi
    if [[ -f "${OUTDIR}/${srr}_2.fastq" ]]; then
        gzip "${OUTDIR}/${srr}_2.fastq"
    fi
    if [[ -f "${OUTDIR}/${srr}.fastq" ]]; then
        gzip "${OUTDIR}/${srr}.fastq"
    fi

    echo "▒~V~R~\~E $srr done."
    touch "$done_flag"
    rm -f "$sra_file"
}

# ▒~V~R~I▒~V~R▒~V~R~G~O▒~V~R~D▒~V~R~P~F
for srr in "${SRR_LIST[@]}"; do
    download_and_convert "$srr"
done

echo "▒~V~R~_~N~I All downloads and conversions completed."

3. 质控：fastqc和multiqc

4. 清洗：fastp（去低质、去接头）

5. 输出：质控报告、清洗后数据

**思考题**：

如果质控报告中碱基质量偏低，如何调整fastp的参数？

### 模块3：比对与定量

**目标**：基于参考基因组构建索引并进行比对与定量。

**任务**：

1. 下载物种参考基因组和基因注释（GTF文件）

2. 用hisat2-build构建索引，hisat2比对

3. 使用stringtie进行转录本组装和基因定量

4. 输出：比对结果（BAM）、定量矩阵

**思考题**：

为什么需要构建参考基因组索引？

### 模块4：数据分析与可视化

**目标**：对定量数据进行质控并与GEO原作者的表达矩阵比较。

**任务**：

1. 表达矩阵的log2转换和标准化

2. 热图和PCA分析

3. 输出：可视化结果（热图、PCA图）

**思考题**：

如果PCA中样本未能按组分离，原因可能是什么？

### 模块5：差异表达分析与功能富集

**目标**：差异表达分析与功能富集，理解生物学意义。

**任务**：

1. DESeq2筛选差异表达基因

2. clusterProfiler进行GO/KEGG富集分析

3. 输出：差异表达结果、富集结果

**思考题**：

为什么要进行多重检验校正？

### 模块6：结果总结与项目展示

**目标**：总结分析结果，展示项目成果。

**任务**：

1. 撰写Markdown报告，包含背景、方法、结果、讨论

2. 在GitHub创建仓库，上传项目

GitHub仓库示例：

```
project_rnaseq/
│── data/
│── scripts/
│── results/
│── README.md
│── report.md
```

**思考题**：

如何优化这次RNA-seq分析项目？








# fastq转换格式
# fastp质控
# 下载物种参考基因组和基因注释（GTF文件）
# 用hisat2-build构建索引，hisat2比对
# 使用stringtie进行转录本组装和基因定量
# 输出：比对结果（BAM）、定量矩阵
# DESeq2进行差异表达分析。
#载入必要包
library(DESeq2)

#读取计数矩阵，跳过前几行注释
counts <- read.table("gene_counts.txt", header=TRUE, row.names=1, comment.char="#")

#删除注释列，featureCounts默认多了几列非样本列（比如 Chr, Start, End）
counts <- counts[, -(1:5)]  # 保留纯计数数据列，具体列数根据实际文件调整

#构建样本信息表 (这里请你根据实际样本分组填充)
sample_info <- data.frame(
  row.names = colnames(counts),
  condition = c(rep("mock", 9), rep("SARS-CoV2", 9), rep("SARS-CoV1", 9)) # 举例
)

#构建DESeq2数据集对象
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = sample_info,
                              design = ~ condition)

#运行差异分析
dds <- DESeq(dds)

#获取差异表达结果，比较组名请替换为你感兴趣的组
res <- results(dds, contrast=c("condition", "SARS-CoV2", "mock"))

#简单查看结果
head(res)

#保存结果
write.csv(as.data.frame(res), file="DEG_SARSCoV2_vs_Mock.csv")

# RNA-seq 分析结果富集与网络可视化

#示例：根据富集结果构建 igraph 网络图
library(clusterProfiler)
library(igraph)
library(tidyverse)

#假设你已经完成了富集分析，拿到了 enrichGO 结果
go_result <- read.csv("go_enrich_result.csv")  # 可改为 enrichGO() 结果对象

#示例构建图（根据 term 与 gene 的关系）
edges <- go_result %>%
  separate_rows(geneID, sep = "/") %>%
  select(Term = Description, Gene = geneID)

#构建 igraph 对象
g <- graph_from_data_frame(edges, directed = FALSE)

#简单可视化
plot(g,
     vertex.label.cex = 0.7,
     vertex.size = 5,
     vertex.label.color = "black",
     edge.color = "grey")


