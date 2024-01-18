#!/bin/sh

## setting the enviornmnent
currpath=$(pwd)
project_home="$HOME/bontempo_pigs_caecum_gastroherb"
outdir="${project_home}/Analysis/micca"
core=8
sickle_exe="$HOME/software/sickle/sickle"

r1="R1" ## delimiter R1
r2="R2" ## delimiter R2

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
	$sickle_exe pe -f "${file}_${r1}_cutadapt.fastq.gz" -r "${file}_${r2}_cutadapt.fastq.gz" -o ${outdir}/trimmed/"${file}_trimmed_${r1}.fastq.gz" -p ${outdir}/trimmed/"${file}_trimmed_${r2}.fastq.gz" -s ${outdir}/trimmed/"${file}_singles.gz" -t sanger -q 20 -g 1>> ${outdir}/quality_control/stats_trim.txt
	## pe: paired-end reads
	## -t: type of quality score (Illumina, Solexa, Sanger): Sanger quality uses a Phred score in [0, 93] ## more here: https://en.wikipedia.org/wiki/FASTQ_format#Quality
	## -q: quality threshold (Phred score) [default: 20]
	## -g: gzip output

done < names_single.txt

echo "DONE!!"

