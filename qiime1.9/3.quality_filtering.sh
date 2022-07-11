#!/bin/sh

## script that applies quality filters and assembles paired reads into a single Fasta file
## input data format is fastq files from 2.join_paired_ends.sh
## 

## setting up the environment
currpath=$(pwd)
project_home="$HOME/useful_microbiome_predipping"
input_folder="Analysis/qiime1.9/2.join_paired_ends"
output_dir="Analysis/qiime1.9/3.quality_filtering"
#sing_container="${project_home}/Qiime1.9.sif"
sing_container="$HOME/software/qiime_docker:fischuu-qiime-1.9.1.sif"
paramfile="Config/parameters/filter.parameters"

cd $project_home
echo "project folder is $project_home"

## make folder if it does not exist
if [ ! -d "${output_dir}" ]; then
	mkdir -p ${output_dir}
	chmod g+rxw ${output_dir}
fi

## generating qiime_parameters file
#echo "split_libraries_fastq:max_bad_run_length 3 split_libraries_fastq:min_per_read_length_fraction 0.75 split_libraries_fastq:sequence_max_n 0 split_libraries_fastq:phred_quality_threshold 19" >| "$PWD/qiime_parameters_quality"

## using the Singularity container
echo " - calling the singularity container for quality filtering"
singularity run ${sing_container} multiple_split_libraries_fastq.py -p "${project_home}/$paramfile" --input_dir=${project_home}/$input_folder \
	--read_indicator=_R1 --include_input_dir_path --output_dir=${project_home}/$output_dir
chmod g+rxw -R ${project_home}/$output_dir

## fix sample names in the output fasta file (<prefix>_ needs to be removed, otherwise prefix will be the sample name)
cd ${project_home}/$output_dir
#sed -i 's/seqprep_/s/g' seqs.fna

echo "DONE!

