#!/bin/sh

## setting the enviornmnent
currpath=$(pwd)
project_home="$HOME/bontempo_pigs_rectum"
analysis_dir="${project_home}/Analysis/micca"
inpdir="${analysis_dir}/clean"
outdir="${analysis_dir}/otu_table"
sing_container="$HOME/software/micca_latest.sif"
core=8

if [ ! -d ${outdir} ]; then
	mkdir -p ${outdir}
fi

## OTU picking
# pick otu
echo " - producing the OTU table"
singularity run $sing_container micca otu -m denovo_unoise -i ${inpdir}/assembled_16S_clean.fasta -o ${outdir}/ -t 8 --rmchim

echo "DONE!!"

