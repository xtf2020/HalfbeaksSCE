####### Calculate gene-r2 using VCFtools

grep -w -i 'male' popmap |cut -f1 >male.ind
cat popmap |cut -f1 >all.ind

######### calculate r2 in the homozygous sex
vcftools --gzvcf vcf \
	--thin 10000 --maf 0.15 \
	--keep male.ind --max-missing 1 \
	--chr chr05 --geno-r2 -c|\
	bgzip -@2 >male.ld.gz
    
######### calculate r2 for all inviduals include males and females
vcftools --gzvcf vcf \
	--thin 10000 --maf 0.15 \
	--keep all.ind --max-missing 1 \
	--chr chr05 --geno-r2 -c|\
	bgzip -@2 >all.ld.gz

######## transfer site r2 to window r2

Site_r2_windowLD.r male.ld
mv Win_LD male_win_ld ####### rename the results file

sed -n '1d' Win_LD |\
	awk '{print $2"\t"$1"\t"$3}'>all.ld
	
cat male.ld all.ld >combine.ld

