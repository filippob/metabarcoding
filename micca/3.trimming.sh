#!/bin/sh

## setting the enviornmnent
currpath=$(pwd)
project_home="$HOME/useful_microbiome_chlorine"
outdir="${project_home}/Analysis/micca"
core=8
sickle_exe="$HOME/software/sickle/sickle"

if [ ! -d "${outdir}/trimmed" ]; then
	mkdir -p ${outdir}/trimmed
fi

## TRIMMING
# trim low quality part. 
# Q = 20 inspect quality, eventually for 16S Q can be set to 25

cd ${outdir}/cutadapt

while read file
do
	echo "Running sickle on file "${file}""
	echo "Running sickle on file "${file}"" >> ${outdir}/quality_control/stats_trim.txt
	## !!! watch out with the read indicator (e.g. R1/R2 or 1/2) !!!
	$sickle_exe pe -f "${file}_1_cutadapt.fastq.gz" -r "${file}_2_cutadapt.fastq.gz" -o ${outdir}/trimmed/"${file}_trimmed_1.fastq.gz" -p ${outdir}/trimmed/"${file}_trimmed_2.fastq.gz" -s ${outdir}/trimmed/"${file}_singles.gz" -t sanger -q 20 -g 1>> ${outdir}/quality_control/stats_trim.txt
done < names_single.txt

echo "DONE!!"

