#!/bin/sh

## script that removes primers and adapaters
## follows the ampliseq structure described in the link below: 
## https://bioinformateachers.github.io/bioinformatics/microbiome/2022/06/11/quick-ampliseq.html

set -x

## setting the enviornmnent
currpath=$(pwd)
project_home="$HOME/bontempo_pigs_rectum"
outdir="${project_home}/Analysis/micca"
core=8
## primer sequences (primers follow the adapters)
fwd_primer="CCTACGGGNGGCWGCAG" #(adapter: TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG)
rev_primer="GACTACHVGGGTATCTAATCC" #(adapter: GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG)

export PATH=/home/biscarinif/.local/bin:$PATH

## make folders
if [ ! -d "${outdir}/cutadapt" ]; then
	mkdir -p ${outdir}/cutadapt
fi

if [ ! -d "${outdir}/trimmed" ]; then
	mkdir -p ${outdir}/trimmed
fi

if [ ! -d "${outdir}/join" ]; then
	mkdir -p ${outdir}/join
fi

if [ ! -d "${outdir}/quality_control" ]; then
	mkdir -p ${outdir}/quality_control
fi


## Create a file of names that will be used for looping. Only file/sample name, remove extension and R1/R2

cd "${outdir}/1.renamed"

## remove file with sample/file names if already present (old files, we don't want to append new lines)
#if [ -f "names.txt" ];
#	rm names.txt
#fi

for i in *.fastq.gz
do
	echo "$i" | cut -d "_" -f1 >> names.txt
done

sed 'n; d' names.txt > names_single.txt

cp names_single.txt ${outdir}/trimmed
cp names_single.txt ${outdir}/cutadapt
cp names_single.txt ${outdir}/join

# remove primers with cutadapt
# Primers (Sequences from Pindo and FMACH)
# forward: CCTACGGGNGGCWGCAG
# reverse: GACTACNVGGGTWTCTAATCC

cd "${outdir}/1.renamed"
echo $(pwd)

while read file
do
	echo "Running cutadapt on file '${file}'"
	## !!! watch out for the read indicator (e.g. R1/R2, or just 1/2) !!!
	cutadapt -g Forward=${fwd_primer} -G Reverse=${rev_primer} --discard-untrimmed --pair-filter=any -o "${outdir}/cutadapt/${file}_R1_cutadapt.fastq.gz" -p "${outdir}/cutadapt/${file}_R2_cutadapt.fastq.gz" "${file}_R1.fastq.gz" "${file}_R2.fastq.gz" >> "${outdir}/quality_control/cutadapt_report.txt"  

done < names_single.txt

# --discard-untrimmed, --trimmed-only
#                        Discard reads that do not contain an adapter.

#--pair-filter=(any|both|first)
#                        Which of the reads in a paired-end read have to match
#                        the filtering criterion in order for the pair to be
#                        filtered. Default: any


echo "DONE!!"

