#!/bin/bash

## script to run FastQC/MultiQC on specified samples

###########################
## setting the environmnent
###########################

## folder paths
currpath=$(pwd)
datapath="/home/ngs/220713_M04028_0149_000000000-KFG5W/"
outdir="$HOME/bontempo_pigs_rectum/Analysis/0.fastqc"
temp_folder="$HOME/temp"

## samples to select
## if start/end is not defined, prefix is used to exclude files when copying
sample_start=15 #first sample to use (in the sequence)
sample_end=104 #last sample to use (in the sequence)
prefix=""

## software
fastqcsing="/home/core/nxf_singularity_cache/depot.galaxyproject.org-singularity-fastqc-0.11.9--0.img"
multiqcsing="/home/core/nxf_singularity_cache/quay.io-biocontainers-multiqc-1.9--py_1.img"

## computing parameters
singularity=1 ## 1 to use the singularity containers (slurm on minicluster); 0 to use local software
core=8
condapath="/usr/local/Miniconda3-py38_4.8.3-Linux-x86_64/etc/profile.d/conda.sh" ## to export conda functions (e.g. activate/deactivate)
condaenv="mycobiota"

## Create analysis folders

echo $HOME
echo $currpath

## create the temporary folder where fastq files are place and processed
echo " - creating temporary folder for sequence data ${temp_folder}"
if [ ! -d "${temp_folder}" ]; then
        mkdir -p ${temp_folder}
	chmod g+rwx ${temp_folder}
fi

## clean the temporary folder (in case there were previous left overs)
if [ -d "${temp_folder}" ]; then
        echo " - cleaning the temp/temp_fastq folder"
        cd ${temp_folder}
        rm *
fi

## copying files of interest to process
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
if [ $singularity == 1 ]; then
	echo "running FASTQC through singularity"
	$fastqcsing fastqc ${temp_folder}/*.fastq.gz -o ${temp_folder} -t 8
else
	echo "running FASTQC from locally installed software"
	$HOME/software/FastQC/fastqc ${temp_folder}/*.fastq.gz -o ${temp_folder} -t 8
fi

## MULTI-QC
cd $currpath
if [ ! -d $outdir ]; then
	mkdir -p $outdir
fi

cd ${temp_folder}

if [ $singularity == 1 ]; then
	echo "running MULTIQC through singularity"
	## DEPENDING ON THE VERSION OF THE CONTAINER YOU MAY USE ONE OF THE FOLLOWING SYNTAX LINES
        #singularity run $multiqcsing multiqc .
        #singularity run $multiqcsing .
	$multiqcsing multiqc .
else
	echo " - activate the conda env and run multiqc"
	source $condapath
	conda activate $condaenv

	echo "running MULTIQC through locally installed software"
	multiqc .
fi

cp multiqc_report.html $outdir/
cp -r multiqc_data $outdir/

echo "DONE!"

