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

 # fastqc
#!/bin/bash
# fastqc_all.sh：质控所有 .fastq.gz 文件

set -euo pipefail

FASTQ_DIR="/mnt/alamo01/users/yuansongwei7/project_rnaseq/virus_projects/SARS-CoV-2/GSE255647/data/fastq"
OUT_DIR="/mnt/alamo01/users/yuansongwei7/results/fastqc/raw"

mkdir -p "$OUT_DIR"
cd "$FASTQ_DIR"

for fq in *.fastq.gz; do
    echo "Running fastqc for $fq ..."
    fastqc -o "$OUT_DIR" "$fq"
done

# multiqc
#!/bin/bash
set -euo pipefail

#multiqc 汇总 fastp 结果
FASTP_RESULT_DIR="/mnt/alamo01/users/yuansongwei7/project_rnaseq/virus_projects/SARS-CoV-2/GSE255647/results/fastp/clean"
OUT_DIR="/mnt/alamo01/users/yuansongwei7/project_rnaseq/virus_projects/SARS-CoV-2/GSE255647/results/multiqc"

mkdir -p "$OUT_DIR"

cd "$FASTP_RESULT_DIR"
multiqc . -o "$OUT_DIR"

echo "MultiQC report generated at $OUT_DIR"



6. 清洗：fastp（去低质、去接头）

#!/bin/bash
# 在download_data/GSE255647目录下运行

OUTDIR="/mnt/alamo01/users/yuansongwei7/download_data/fastq"
SRADIR="/mnt/alamo01/users/yuansongwei7/download_data/GSE255647"
mkdir -p "$OUTDIR" || { echo "wrong outdir"; exit 1; }
THREADS=4

for sra in "$SRADIR"/*.sra; do
    SRR=$(basename "$sra" .sra)
    echo "===== 转换 $SRR ====="

    # 检查SRA文件是否存在
    if [ ! -f "$sra" ]; then
        echo "错误：文件 $sra 不存在，跳过"
        continue
    fi

    fastq-dump "$sra" \
      --split-files \
      --gzip \
      --outdir "$OUTDIR"

    # 检查压缩后的文件是否存在
    if [[ -f "$OUTDIR/${SRR}_1.fastq.gz" && -f "$OUTDIR/${SRR}_2.fastq.gz" ]]; then
        echo "✓ $SRR 转换成功"
    else
        echo "❌ $SRR 转换失败！请检查"
    fi
done


6. 输出：质控报告、清洗后数据

**思考题**：

如果质控报告中碱基质量偏低，如何调整fastp的参数？

### 模块3：比对与定量

**目标**：基于参考基因组构建索引并进行比对与定量。
project_rnaseq/
├── mapping/
│   └── hisat2/
│       ├── SRR27961772.sorted.bam
│       └── ...
├── quantification/
│   └── stringtie/
│       ├── SRR27961772.gtf
│       └── ...
├── genome_index/
│   └── hg38/
│       └── gencode.v43.annotation.gtf
└── scripts/
    ├── stringtie_all.sh
    ├── SRR27961772_stringtie.sh
    └── logs/
        ├── SRR27961772_stringtie.o
        └── SRR27961772_stringtie.e
**任务**：

1. 下载物种参考基因组和基因注释（GTF文件）
   #!/bin/bash
set -euo pipefail

#下载与保存路径
GENOME_DIR="/mnt/alamo01/users/yuansongwei7/genome_index/hg38"
mkdir -p "$GENOME_DIR"

echo ">>> 开始下载人类参考基因组（FASTA）..."
FA_URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/GRCh38.primary_assembly.genome.fa.gz"
FA_GZ="$GENOME_DIR/GRCh38.primary_assembly.genome.fa.gz"
FA_FILE="$GENOME_DIR/GRCh38.primary_assembly.genome.fa"

#下载FASTA
wget -c -O "$FA_GZ" "$FA_URL"

#解压
echo ">>> 解压参考基因组..."
gunzip -c "$FA_GZ" > "$FA_FILE"

echo ">>> 下载人类基因注释文件（GTF）..."
GTF_URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/gencode.v43.annotation.gtf.gz"
GTF_GZ="$GENOME_DIR/gencode.v43.annotation.gtf.gz"
GTF_FILE="$GENOME_DIR/gencode.v43.annotation.gtf"

#下载GTF
wget -c -O "$GTF_GZ" "$GTF_URL"

#解压
echo ">>> 解压GTF文件..."
gunzip -c "$GTF_GZ" > "$GTF_FILE"

#检查
echo ">>> 检查文件是否存在..."
[[ -f "$FA_FILE" ]] && echo "✔ FASTA 准备完成：$FA_FILE" || echo "✘ FASTA 文件缺失！"
[[ -f "$GTF_FILE" ]] && echo "✔ GTF 准备完成：$GTF_FILE" || echo "✘ GTF 文件缺失！"


3. 用hisat2-build构建索引，hisat2比对

   #!/bin/bash
# hisat2索引构建脚本（适用于人类参考基因组 hg38）

set -euo pipefail

echo ">>> 开始构建 hisat2 索引..."

#路径配置
GENOME_FA="/mnt/alamo01/users/yuansongwei7/genome_index/hg38/GRCh38.primary_assembly.genome.fa"
OUTDIR="/mnt/alamo01/users/yuansongwei7/genome_index/hg38"
BASENAME="hg38"

#创建输出目录（如果不存在）
mkdir -p "$OUTDIR"

#构建索引
hisat2-build -p 16 "$GENOME_FA" "$OUTDIR/$BASENAME"

echo ">>> hisat2 索引构建完成！"


# hisatall。sh

#!/bin/bash
# hisat2_all.sh：自动批量提交比对作业
set -euo pipefail

# 输入输出路径
FASTQ_DIR="/mnt/alamo01/users/yuansongwei7/download_data/GSE255647_analysis/cleaned_data"
OUTDIR="/mnt/alamo01/users/yuansongwei7/project_rnaseq/mapping/hisat2"
INDEX="/mnt/alamo01/users/yuansongwei7/genome_index/hg38/hg38"
LOGDIR="./logs"
mkdir -p "$OUTDIR" "$LOGDIR"

# 启用更严格匹配
shopt -s nullglob

# 遍历所有 *_1.fastq.gz 文件
for f1 in "$FASTQ_DIR"/*_1.fastq.gz; do
    # 获取SRR编号
    base=$(basename "$f1" _1.fastq.gz)
    f2="$FASTQ_DIR/${base}_2.fastq.gz"

    # 检查配对文件存在
    if [[ ! -f "$f2" ]]; then
        echo "❌ 配对文件缺失: $f2"
        continue
    fi

    # 输出文件名
    BAM="$OUTDIR/${base}.sorted.bam"

    # 比对作业脚本内容
    cat <<EOF > "${base}_hisat2.sh"
#!/bin/bash
#PBS -N hisat2_$base
#PBS -q fast
#PBS -l cpu=64
#PBS -l mem=64G
#PBS -o $LOGDIR/hisat2_${base}.o
#PBS -e $LOGDIR/hisat2_${base}.e
#PBS -V
#PBS -cwd

source ~/.bashrc
micromamba activate rnaseq

hisat2 -p 64 -x $INDEX -1 $f1 -2 $f2 \
  | samtools sort -@ 32 -o "$BAM"

samtools index "$BAM"
EOF

    # 提交任务
    qsub "${base}_hisat2.sh"
    echo "✓ 已提交: $base"
done


5. 使用stringtie进行转录本组装和基因定量

#!/bin/bash
set -euo pipefail
#-------------------------------------------------------------------
#Script Name: stringtie_all.sh
#Description: 批量处理 HISAT2 比对结果，进行转录本定量分析。
#Author: ChatGPT for yuansongwei7
#Date: 2025-07-11
#-------------------------------------------------------------------

#设置路径变量
BAM_DIR="/mnt/alamo01/users/yuansongwei7/project_rnaseq/mapping/hisat2"
GTF="/mnt/alamo01/users/yuansongwei7/genome_index/hg38/gencode.v43.annotation.gtf"
OUT_ROOT="/mnt/alamo01/users/yuansongwei7/project_rnaseq/quantification/stringtie"
LOG_DIR="/mnt/alamo01/users/yuansongwei7/scripts/logs"
SCRIPT_DIR="/mnt/alamo01/users/yuansongwei7/scripts"

#创建日志目录
mkdir -p "$OUT_ROOT" "$LOG_DIR"

#遍历所有 BAM 文件
for bam_file in ${BAM_DIR}/*.sorted.bam; do
  base=$(basename "$bam_file" .sorted.bam)
  outdir="${OUT_ROOT}/${base}"
  SCRIPT="${SCRIPT_DIR}/${base}_stringtie.sh"

  mkdir -p "$outdir"

  cat <<EOF > "$SCRIPT"
#!/bin/bash
source ~/.bashrc
micromamba activate rnaseq

mkdir -p "$outdir"

stringtie "$bam_file" \\
  -G "$GTF" \\
  -o "$outdir/out.gtf" \\
  -p 32 -e -B
EOF

  chmod +x "$SCRIPT"

  qsub -cwd -V -l cpu=64:mem=64G -q fast \
       -o "$LOG_DIR/${base}_stringtie.o" \
       -e "$LOG_DIR/${base}_stringtie.e" \
       -N "stringtie_${base}" \
       "$SCRIPT"

  echo "✓ 已提交：$base"
done


6. 输出：比对结果（BAM）、定量矩阵

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


