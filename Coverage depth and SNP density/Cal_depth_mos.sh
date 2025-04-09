#!/bin/bash

bam_dir=$1
pop=$2
win=$3

[ $# -lt 3 ] && echo $0 [bamfile_dir] [popfile] [windowsize] && exit 1

##check absulute path
if [ "${bam_dir:0:1}" = "/" ] && [ "${pop:0:1}" = "/" ]
then
	echo "path ok"
else
	echo "please use absulute path" && exit 1
fi

## produce window depth file using mosdepth for each individual
mkdir mos_res_temp && cd mos_res_temp
for i in `cut -f1 $pop`;do
    echo "calculating individual $i "
    mosdepth -t 12 -x -n -b $win\
        $i $bam_dir/$i.bam
    zcat $i.regions.bed.gz |\
        sed "1i chr\tstart\tend\t$i" > $i.depth
done

## concat male
paste `grep -w 'male' $pop |cut -f1|\
    sed ':label;N;s/\n/.depth /;b label'|\
    sed 's/$/.depth/g'` >male.out.temp

m=`grep -w 'male' $pop|wc -l`

cut -f1,2,3,`seq 1 $m |\
    awk '{print $1*4}'|\
    sed ':label;N;s/\n/,/;b label'` male.out.temp \
    >male.out

## concat female
paste `grep -w 'female' $pop |cut -f1|\
    sed ':label;N;s/\n/.depth /;b label'|\
    sed 's/$/.depth/g'` >female.out.temp

f=`grep -w 'female' $pop|wc -l`

cut -f`seq 1 $f|\
    awk '{print $1*4}'|\
    sed ':label;N;s/\n/,/;b label'` female.out.temp \
    >female.out

paste male.out female.out >mos_depth.out
mv mos_depth.out ../ && cd ../ && rm -r mos_res_temp
