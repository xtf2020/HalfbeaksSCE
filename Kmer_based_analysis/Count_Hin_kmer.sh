
########## the file "ind" contains IDs for all H.intermedius individual used for analysis, each row one individual
for n in `cat ind`
do
    mkdir $n && cd $n
    echo "/mnt/PRJ001/Zhenyu/jx_sex/fq_rn/"$n"_1.fq.gz" >"$n"_fq_lst
    echo "/mnt/PRJ001/Zhenyu/jx_sex/fq_rn/"$n"_2.fq.gz" >>"$n"_fq_lst

    /opt/bio/kmersGWAS/external_programs/kmc_v3 -t35 -k31 -ci3 @"$n"_fq_lst output_kmc_canon ./ 1> kmc_canon.1 2> kmc_canon.2
    /opt/bio/kmersGWAS/external_programs/kmc_v3 -t35 -k31 -ci0 -b @"$n"_fq_lst output_kmc_all ./ 1> kmc_all.1 2> kmc_all.2

    /opt/bio/kmersGWAS//bin/kmers_add_strand_information -c output_kmc_canon -n output_kmc_all -k 31 -o kmers_with_strand
    rm *.kmc*

    cd ../
done
