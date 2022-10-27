#!/bin/sh

## script to prepare a text-file OTU (ASV) table from the output of Micca
## adding taxonomies (results from classification) and header (#OTU ID\ttaxonomy)
## the generated OTU/ASV table is then converted to the biom format
## then it is filtered (min n. of counts/samples)
## and converted back to text format
## the filtered biom file thus generated is ready as input for the R package phyloseq
## (to continue with down stream calculations of alpha and beta diversity and 
## and for the comparison of specific taxa)

## setting the enviornmnent
currpath=$(pwd)
project_home="$HOME/bontempo_pigs_rectum"
analysis_dir="${project_home}/Analysis/micca"
inpdir="${analysis_dir}/taxa_classification"
outdir="${analysis_dir}/filtered_otu"
taxaf="taxa_SILVA.txt"
otuf="${analysis_dir}/otu_table/otutable.txt"

qiime_container="$HOME/software/qiime-1.9.1_latest.sif"
ncounts=20 ## min n. of counts overall
nsamples=4 ## min. n. of samples with non-zero counts

if [ ! -d ${outdir} ]; then
        mkdir $outdir
fi

echo " - preparing the otu_table + taxonomy text file "
## add header to taxonomy from classification (Micca)
echo -e "#OTU ID\ttaxonomy" | cat - $inpdir/$taxaf > $outdir/temp

## make columns: OTU ID and taxonomy
cut -f1 $outdir/temp > $outdir/col_1
cut -f2 $outdir/temp > $outdir/col_n

## remove the OTU column from the Micca otu table (named OTU instead of #OTU ID)
## then add the new OTU column
cut -f2- $otuf > $outdir/temp2
paste --delimiters='\t'  $outdir/col_1 $outdir/temp2 > $outdir/temp

## now 
paste --delimiters='\t' $outdir/temp $outdir/col_n > $outdir/otu_table.csv

## convert - filter - convert
echo " - converting and filtering the OTU table "
singularity run $qiime_container biom convert -i $outdir/otu_table.csv -o $outdir/otu_table.biom --table-type="OTU table" --to-hdf5 --process-obs-metadata taxonomy
singularity run $qiime_container filter_otus_from_otu_table.py -i $outdir/otu_table.biom -n $ncounts -s $nsamples -o $outdir/otu_table_filtered.biom
singularity run $qiime_container biom convert -i $outdir/otu_table_filtered.biom -o $outdir/otu_table_filtered.txt --to-tsv --header-key taxonomy

## house cleaning
cd $outdir
rm temp*
rm col_*

echo "DONE!"
