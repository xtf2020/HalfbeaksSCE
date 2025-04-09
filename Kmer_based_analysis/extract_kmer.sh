for n in `cat /mnt/PRJ001/Zhenyu/jx_sex/kmer/ind |sed -n '1,12p'`
do
 /mnt/PRJ001/Zhenyu/jx_sex/kmer/run_gemma/bbmap/bbduk.sh \
  in1=/mnt/PRJ001/Zhenyu/jx_sex/fq_rn/"$n"_1.fq.gz \
  in2=/mnt/PRJ001/Zhenyu/jx_sex/fq_rn/"$n"_2.fq.gz \
  outm1="$n"_extract_1.fq outm2="$n"_extract_2.fq \
  k=31 mm=f \
  ref=male_specific.kmer_fa \
  stats="$n"_male_kmers_stats.txt
  
done
