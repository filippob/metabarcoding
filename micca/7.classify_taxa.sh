#!/bin/sh

## setting the enviornmnent
currpath=$(pwd)
project_home="$HOME/delower"
analysis_dir="${project_home}/Analysis/micca"
inpdir="${analysis_dir}/otu_table"
outdir="${analysis_dir}/taxa_classification"
sing_container="$HOME/software/micca_latest.sif"

## SELECT REFERENCE DATABASE
dbpath="$HOME/databases/SILVA_132_QIIME_release"
#dbpath="16srrna" ## ribosomal database project: 16srrna, fungallsu, fungalits_warcup, fungalits_unite
#dbpath="$HOME/databases/gg_13_8_otus"
#dbpath="$HOME/databases/SILVA_138" ## no qiime 1.9 format available (find out whether the qiime2 format can be used)

core=8

if [ ! -d ${outdir} ]; then
	mkdir $outdir
fi

<<<<<<< HEAD
## check if SILVA
tmp=`echo $dbpath | grep 'SILVA'`
#echo "tmp is $tmp"

if [ "$tmp" != '' ]; then
	## CLASSIFY WITH VSEARCH AND SILVA!
	echo " - classify with SILVA"
	echo "dbpath is $dbpath"
	# QIIME compatible SILVa DB should be downloaded
	singularity run $sing_container micca classify -m cons -i $inpdir/otus.fasta  -o $outdir/taxa_SILVA.txt --ref $dbpath/rep_set/rep_set_16S_only/97/silva_132_97_16S.fna \
		--ref-tax $dbpath/taxonomy/16S_only/97/taxonomy_7_levels.txt
	# singularity run $sing_container micca classify -m cons -i $inpdir/otus.fasta  -o $outdir/taxa_SILVA.txt --ref $dbpath/SILVA_138_SSURef_NR99_tax_silva.fasta --ref-tax $dbpath/taxonomy
fi
=======
# classify RDP
# export RDPPATH=/home/fvitali/Documenti/Personal_PATH_folder/rdp_classifier_2.13/
# singularity run $sing_container micca classify -m rdp -i $outdir/micca_16S/otus.fasta --rdp-gene 16srrna -o $outdir/micca_16S/taxa.txt
>>>>>>> b2fac252f6f324db90d72f28bd20c54713275db0

## check if greengenes
tmp=`echo $dbpath | grep 'gg_'`
#echo "tmp is $tmp"

if [ "$tmp" != '' ]; then
	## GREENGENES
        echo " - classify with Greengenes"
	singularity run $sing_container micca classify -m cons -i $inpdir/otus.fasta -o $outdir/taxa_gg.txt --ref $dbpath/rep_set/97_otus.fasta --ref-tax $dbpath/taxonomy/97_otu_taxonomy.txt
fi

## else RDP
tmp=`echo $dbpath | grep 'gg_'`
tmp+=`echo $dbpath | grep 'SILVA'`
#echo "tmp is $tmp"

if [ "$tmp" == '' ]; then
	# classify RDP
	echo " - classify with RDP"
	echo "dbpath is $dbpath"
	#export RDPPATH= <rdppath already exported in the Micca singularity container> (/opt/rdp_classifier)
	singularity run $sing_container micca classify -m rdp -i $inpdir/otus.fasta --rdp-gene $dbpath -o $outdir/taxa_rdp.txt
fi

# -m: rdp, cons, otuid
# * cons: VSEARCH-based consensus classifier --> input sequences are searched in the reference database with VSEARCH
# (https://github.com/torognes/vsearch). For each query sequence the method retrives up to 'cons-maxhits' hits (i.e. identity >=
# 'cons-id'). Then, the most specific taxonomic label that is associated with at least 'cons-minfrac' of the hits is assigned. The method is similar to the UCLUST-based consensus
# taxonomy assigner presented in doi: 10.7287/peerj.preprints.934v2 and available in QIIME.


echo "DONE!"

