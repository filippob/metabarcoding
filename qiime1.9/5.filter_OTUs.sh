#!/bin/bash

###Script that performs filtering on by total count across samples greater than 15 of the number of OTUs in at least 2 samples

## setting up the environment
currpath=$(pwd)
project_home="$HOME/useful_microbiome_predipping"
data_folder="Analysis/qiime1.9/4.OTU_picking"
output_dir="Analysis/qiime1.9/5.filter_OTUs"
sing_container="$HOME/software/qiime_docker:fischuu-qiime-1.9.1.sif"

cd $currpath
echo "project folder is $project_home"

## make folder if it does not exist
if [ ! -d "$project_home/${output_dir}" ]; then
        mkdir -p $project_home/${output_dir}
        chmod g+rxw ${project_home}/${output_dir}
fi

## using the Singularity container

echo " - calling the singularity container"

singularity run ${sing_container} filter_otus_from_otu_table.py \
        -i $project_home/${data_folder}/otu_table.biom \
        -n 15 \
        -s 2 \
        -o $project_home/${output_dir}/otu_table_filtered.biom

chmod g+xrw -R ${project_home}/$output_dir

echo " - converting biom file to tsv file"
cd ${project_home}/$output_dir
singularity run ${sing_container} biom convert -i otu_table_filtered.biom -o otu_table_filtered.txt --to-tsv --header-key taxonomy

echo "DONE!"

