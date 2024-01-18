#!/bin/sh

## setting the enviornmnent
currpath=$(pwd)
project_home="$HOME/bontempo_pigs_rectum"
analysis_dir="${project_home}/Analysis/micca"
inpdir="${analysis_dir}/otu_table"
outdir="${analysis_dir}/taxa_classification"
sing_container="$HOME/software/micca_latest.sif"
dbpath="$HOME/databases/SILVA_132_QIIME_release"
core=8

if [ ! -d ${outdir} ]; then
	mkdir $outdir
fi

# classify RDP
# export RDPPATH=/home/fvitali/Documenti/Personal_PATH_folder/rdp_classifier_2.13/
# singularity run $sing_container micca classify -m rdp -i $outdir/micca_16S/otus.fasta --rdp-gene 16srrna -o $outdir/micca_16S/taxa.txt

## CLASSIFY WITH VSEARCH AND SILVA!
# QIIME compatible SILVa DB should be downloaded
echo " - classify with SILVA 132"
$sing_container micca classify -m cons -i $inpdir/otus.fasta  -o $outdir/taxa_SILVA.txt --ref $dbpath/rep_set/rep_set_16S_only/97/silva_132_97_16S.fna --ref-tax $dbpath/taxonomy/16S_only/97/taxonomy_7_levels.txt

# -m: rdp, cons, otuid
# * cons: VSEARCH-based consensus classifier --> input sequences are searched in the reference database with VSEARCH
# (https://github.com/torognes/vsearch). For each query sequence the method retrives up to 'cons-maxhits' hits (i.e. identity >=
# 'cons-id'). Then, the most specific taxonomic label that is associated with at least 'cons-minfrac' of the hits is assigned. The method is similar to the UCLUST-based consensus
# taxonomy assigner presented in doi: 10.7287/peerj.preprints.934v2 and available in QIIME.


echo "DONE!"

