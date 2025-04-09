

/opt/bio/kmersGWAS/bin/filter_kmers -t ../../all_kmers_table -k male_spe_kmer_lst -o male_kmer_table

awk '{               
    sum1=0; sum2=0;
    for(i=2;i<=13;i++) sum1+=$i;
    for(i=14;i<=25;i++) sum2+=$i;
    print $0, sum1, sum2
}' male_kmer_table > male_kmer_table_cal

awk '$26>=11 && $27<=2 {print $0}' male_kmer_table_cal >male_specific.kmer_table

cat male_specific.kmer_table |awk '{print ">"NR"\t"$1}'|sed 's/\t/\n/' > ./male_spe/male_specific.kmer_fa

