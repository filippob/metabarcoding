#!/bin/bash

## script to run FastQC/MultiQC on specified samples

###########################
## setting the environmnent
###########################

## folder paths
currpath=$(pwd)
datapath="/home/ngs/Macrogen/1.Rawdata/"
outdir="useful_microbiome_chlorine/Analysis/0.fastqc"
temp_folder="$HOME/temp"

## samples to select
sample_start= #first sample to use (in the sequence)
sample_end= #last sample to use (in the sequence)
prefix="RAL"

## computing parameters
core=8
condapath="/usr/local/Miniconda3-py38_4.8.3-Linux-x86_64/etc/profile.d/conda.sh" ## to export conda functions (e.g. activate/deactivate)
condaenv="mycobiota"

## Create analysis folders

echo $HOME
echo $currpath

## copying files of interest to process
echo " - creating temporary folder for sequence data ${temp_folder}"
if [ ! -d "${temp_folder}" ]; then
        mkdir -p ${temp_folder}
	chmod g+rwx ${temp_folder}
fi

echo " - copying relevant fastq files from data folder 1"
if [ ! -z ${sample_start} ]
then
	for i in $(seq ${sample_start} ${sample_end});
	do
		echo "file ${datapath}/${i}_"
		rsync -av ${datapath}/${i}_*.fastq.gz ${temp_folder}
	done
else
	rsync -av --exclude "${prefix}*" ${datapath}/*.fastq.gz ${temp_folder}
fi

## FastQC
echo " - running FastQC"
$HOME/software/FastQC/fastqc ${temp_folder}/*.fastq.gz -o ${temp_folder} -t 8

## MULTI-QC
cd $currpath
if [ ! -d $outdir ]; then
	mkdir -p $outdir
fi

echo " - activate the conda env and run multiqc"
source $condapath
conda activate $condaenv
cd ${temp_folder}
multiqc .

cp multiqc_report.html ../$outdir/
cp -r multiqc_data ../$outdir/

## clean the temporary folder (in case there were previous left overs)
#echo " - cleaning the temp/temp_fastq folder"
#cd ${temp_folder}
#rm -rf *
#cd $currpath

echo "DONE!"

