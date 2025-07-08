# 下载脚本
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


