########## /usr/bin/bash

######## identify the single copy genes among
#H. intermedius, Z/W of H. sajori, 
# Rhynchorhamphus georgii, 
# Hemiramphus far, 
# Hemiramphus archipelagicus, 
# H. gernaerti, H. limbatus, H. paucirastris

orthofinder -f /mnt/gene_tree/all_pep_dir/ -t 20

######## align cds sequence for sing_copy_gene
bash ./aln_scog.sh/ single_copy_gene.id all_cds.far

###### construct phylogenic tree for each single copy group
for n in `cat sing_copy_gene.id`
do 
	cp /align_result/"$n"_final_mask_align_NT.aln ./ind_tree/
	iqtree2 -s ./ind_tree/"$n"_final_mask_align_NT.aln --quiet
done

##### concat trees into one file

cat ./ind_tree/"$n"_final_mask_align_NT.aln.treefile >>sex_all_tree

########## simple the node infomation of tree file
sed 's/_[0-9]*-T[0-9]//g' all_tree > all_tree_rn

python Count_tree_topology.py sex_all_tree_rn