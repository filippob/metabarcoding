#!/bin/sh

## setting the enviornmnent
currpath=$(pwd)
project_home="$HOME/useful_microbiome_chlorine"
outdir="${project_home}/Analysis/micca/micca_16S"
taxaf="taxa_SILVA.txt"
otuf="otutable.txt"
qiime_container="$HOME/software/qiime_docker:fischuu-qiime-1.9.1.sif"

echo " - preparing the otu_table + taxonomy text file "
## add header to taxonomy from classification (Micca)
echo -e "#OTU ID\ttaxonomy" | cat - $outdir/$taxaf > $outdir/temp

## make columns: OTU ID and taxonomy
cut -f1 $outdir/temp > $outdir/col_1
cut -f2 $outdir/temp > $outdir/col_n

## remove the OTU column from the Micca otu table (named OTU instead of #OTU ID)
## then add the new OTU column
cut -f2- $outdir/$otuf > $outdir/temp2
paste --delimiters='\t'  $outdir/col_1 $outdir/temp2 > $outdir/temp

## now 
paste --delimiters='\t' $outdir/temp $outdir/col_n > $outdir/otu_table.csv

## house cleaning
rm $outdir/temp*

## convert - filter - convert
echo " - converting and filtering the OTU table "
singularity run $qiime_container biom convert -i $outdir/otu_table.csv -o $outdir/otu_table.biom --table-type="OTU table" --to-hdf5 --process-obs-metadata taxonomy
singularity run $qiime_container filter_otus_from_otu_table.py -i $outdir/otu_table.biom -n 15 -s 2 -o $outdir/otu_table_filtered.biom
singularity run $qiime_container biom convert -i $outdir/otu_table_filtered.biom -o $outdir/otu_table_filtered.txt --to-tsv --header-key taxonomy

echo "DONE!"
