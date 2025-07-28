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
project_home="$HOME/beet-its-2025"
analysis_dir="${project_home}"
inpdir="${analysis_dir}/data"
outdir="${analysis_dir}/filtered_features"
featurefile="feature-table.tsv"
taxaf="taxonomy.tsv"

qiime_container="$HOME/software/qiime-1.9.1_latest.sif"
ncounts=20 ## min n. of counts overall
nsamples=3 ## min. n. of samples with non-zero counts

if [ ! -d ${outdir} ]; then
        mkdir $outdir
fi

echo " - preparing the otu_table + taxonomy text file "
## remove the original header from the taxonomy file
tail -n +2 "$inpdir/$taxaf" > $outdir/temp
## add the correct header to taxonomy from classification
echo -e "#OTU ID\ttaxonomy" | cat - $outdir/temp > $outdir/temp2

## make columns: OTU ID and taxonomy
cut -f1 $outdir/temp2 > $outdir/col_1
cut -f2 $outdir/temp2 > $outdir/col_n

## now 
paste --delimiters='\t' $inpdir/$featurefile $outdir/col_n > $outdir/feature-table.csv

## convert - filter - convert
echo " - converting and filtering the OTU table "
singularity run $qiime_container biom convert -i $outdir/feature-table.csv -o $outdir/feature-table.biom --table-type="OTU table" --to-hdf5 --process-obs-metadata taxonomy

## !! THE CODE BELOW NEEDS TO BE CONVERTED FROM QIIME 1.9 TO QIIME 2 !! ## 
#singularity run $qiime_container filter_otus_from_otu_table.py -i $outdir/feature-table.biom -n $ncounts -s $nsamples -o $outdir/feature-table_filtered.biom
#singularity run $qiime_container biom convert -i $outdir/feature-table_filtered.biom -o $outdir/feature-table_filtered.txt --to-tsv --header-key taxonomy

## house cleaning
cd $outdir
rm temp*
rm col_*

echo "DONE!"
