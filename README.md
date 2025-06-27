# homework 下载脚本
#!/bin/bash
set -euo pipefail
export PATH=/mnt/alamo01/yuansongwei7@mgt01:/mnt/alamo01/users/yuansongwei7:$PATH
# 
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
