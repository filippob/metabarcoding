#!/bin/sh

## script that generates the OTU/ASV table
## you can choose which (OTU or ASV) by
## modifying the code below 

set -x

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
## -m: denovo_greedy, denovo_unoise, denovo_swarm, closed_ref, open_ref
# * denovo_greedy clustering: useful for for the identification of 97% OTUs;
# * denovo_unoise: denoise Illumina sequences using the UNOISE3 protocol;
# * denovo_swarm): a robust and fast clustering method (deprecated, it will be removed in version 1.8.0);
# * closed_ref: closed-reference clustering, sequences are clustered against an external reference database and reads that could not be matched are discarded
# * open-reference clustering (open_ref): sequences are clustered against an external reference database and reads that could not be matched are clustered with the 'de novo greedy' protocol

echo "DONE!!"

