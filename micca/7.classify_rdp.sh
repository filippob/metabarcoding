#!/bin/sh

## setting the enviornmnent
currpath=$(pwd)
project_home="$HOME/useful_microbiome_chlorine"
outdir="${project_home}/Analysis/micca"
sing_container="$HOME/software/micca.sif"
dbpath="$HOME/databases/SILVA_132_QIIME_release"
core=8

if [ ! -d "${outdir}/micca_16S" ]; then
	mkdir $outdir/micca_16S
fi


# classify RDP
# export RDPPATH=/home/fvitali/Documenti/Personal_PATH_folder/rdp_classifier_2.13/
# singularity run $sing_container micca classify -m rdp -i $outdir/micca_16S/otus.fasta --rdp-gene 16srrna -o $outdir/micca_16S/taxa.txt

## CLASSIFY WITH VSEARCH AND SILVA!
# QIIME compatible SILVa DB should be downloaded
echo " - classify with SILVA 132"
$sing_container micca classify -m cons -i $outdir/micca_16S/otus.fasta  -o $outdir/micca_16S/taxa_SILVA.txt --ref $dbpath/rep_set/rep_set_16S_only/97/silva_132_97_16S.fna --ref-tax $dbpath/taxonomy/16S_only/97/taxonomy_7_levels.txt


echo "DONE!"

