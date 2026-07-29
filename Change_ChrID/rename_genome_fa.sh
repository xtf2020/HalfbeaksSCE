ori_fasta=$1
rn_file=$2
out_fasta=$3


if [ $# -lt 3 ] ;then
    echo "$0 [ori_fasta] [rn_file] [out_fasta]"
    exit 0
fi

if [ ! -f "${ori_fasta}.fai" ]; then
    samtools faidx $ori_fasta
fi

# 清空输出文件（如果存在）
> $out_fasta

while read a b c;do
    if [ "$c" = "rev" ]; then
        # 对于 rev，提取序列、改名、反向互补
        seqkit grep -n -r -p "$a" $ori_fasta | \
        seqkit replace -p ".+" -r "$b" | \
        seqkit seq -r -p >> $out_fasta
    else
        # 对于 for，只提取序列和改名
        seqkit grep -n -r -p "$a" $ori_fasta | \
        seqkit replace -p ".+" -r "$b" >> $out_fasta
    fi
done < $rn_file