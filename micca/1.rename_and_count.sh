#!/bin/sh

#################################################
## Script to rename samples and count input reads
#################################################

## setting the environmnent
currpath=$(pwd)
project_home="$HOME/useful_microbiome_chlorine"
sing_container="$HOME/software/qiime_docker:fischuu-qiime-1.9.1.sif"
datapath="$HOME/temp"
outdir="${project_home}/Analysis/micca/1.renamed"
core=8

## Rename samples
cd $project_home

if [ ! -d $outdir ]; then
	mkdir -p $outdir
fi

cd $datapath
este=""

for i in *.fastq.gz
do
  sample=$(echo "$i" | cut -d "_" -f1 | cut -d "-" -f1,2,3)
  read=$(echo "$i" | cut -d "_" -f2) ## this will change depending on the structure of sample names !!
  echo $sample
  cp $i ${outdir}/$sample"_"$read$este
  echo -e "$i\t-->\t$sample"_"$read$este" >> ${outdir}/log_renamer.txt
done

## Count reads

cd ${outdir}

for i in *.fastq.gz
do
        echo -n $i >> seq_count_16S_raw.txt
        echo -n " " >> seq_count_16S_raw.txt
        echo $(zcat $i | wc -l) / 4 | bc >> seq_count_16S_raw.txt
done

echo "DONE!"

