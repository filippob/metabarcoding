#!/bin/sh

## setting the enviornmnent
currpath=$(pwd)
project_home="$HOME/useful_microbiome_chlorine"
outdir="${project_home}/Analysis/micca"
sing_container="$HOME/software/micca.sif"
core=8

if [ ! -d "${outdir}/micca_16S" ]; then
	mkdir -p ${outdir}/micca_16S
fi

## OTU picking
# pick otu
echo " - producing the OTU table"
singularity run $sing_container micca otu -m denovo_unoise -i ${outdir}/micca_16S/WP1_assembled_16S.fasta -o ${outdir}/micca_16S/ -t 8 --rmchim

echo "DONE!!"

