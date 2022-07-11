#!/bin/sh

## script that removes barcodes from already demultiplexed data
## input data format is fastq: R1 and R2 files (forward and reverse)
## by default the first 8 bases are trimmed: otherwise, pass a manifest/mapping file with the barcode sequences

## setting up the environment
currpath=$(pwd)
project_home="$HOME/useful_microbiome_predipping"
#data_folder="data/220225_M04028_0144_000000000-K6CMG"
output_dir="Analysis/qiime1.9/1.extract_barcode"
#sing_container="${project_home}/Qiime1.9.sif"
sing_container="$HOME/software/qiime_docker:fischuu-qiime-1.9.1.sif"
temp_folder="$HOME/temp"

cd $project_home
echo "project folder is $project_home"

## make output folder if it does not exist
if [ ! -d "${output_dir}" ]; then
	mkdir -p ${output_dir}
	chmod g+rxw ${output_dir}
fi

## using the Singularity container
echo " - calling the singularity container"
singularity run ${sing_container} multiple_extract_barcodes.py --input_dir=${temp_folder} --output_dir=${output_dir} --read1_indicator _R1 --read2_indicator _R2

echo " - removing barcodes files"
cd $output_dir
find . -name \*barcodes.fastq -type f -delete
cd $currpath

echo " - removing temprorary files"
rm ${temp_folder}/*

echo "DONE!"

