orth_dir=$1
all_cds=$2

sc_id=$orth_dir/Orthogroups/Orthogroups_SingleCopyOrthologues.txt
sc_gene_lst=$orth_dir/Orthogroups/Orthogroups.txt

mkdir unalign_cds_fa

for n in `cat $sc_id`
do
    seqkit faidx --quiet -j40 -r -l <(grep "$n" $sc_gene_lst |cut -d' ' -f 2- |sed 's/ /\n/g') \
	$all_cds >./unalign_cds_fa/$n.fa
done

mkdir align_result
for n in `cat $sc_id`
do 
    echo "macse --java_mem 2000m --in_seq_file ./unalign_cds_fa/$n.fa --genetic_code_number 1 --no_prefiltering --replace_FS_by_gaps --out_file_prefix $n --out_dir ./align_result/$n" >>cds_aln.cmd
done

ParaFly -c cds_aln.cmd


