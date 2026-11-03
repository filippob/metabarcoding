#!/bin/sh

## setting the environment
currpath=$(pwd)
project_home="$HOME/tunisian_rats/colitis"
analysis_dir="${project_home}/Analysis/micca"
inpdir="${analysis_dir}/otu_table"
outdir="${analysis_dir}/taxa_classification"
sing_container="$HOME/software/micca_latest.sif"

## SELECT REFERENCE DATABASE
dbpath="$HOME/databases/SILVA_132_QIIME_release"
#dbpath="$HOME/databases/Silva_138_2"
#dbpath="16srrna" ## ribosomal database project: 16srrna, fungallsu, fungalits_warcup, fungalits_unite
#dbpath="$HOME/databases/gg_13_8_otus"
#dbpath="$HOME/databases/SILVA_138" ## no qiime 1.9 format available (find out whether the qiime2 format can be used)

core=8

mkdir -p "${outdir}"

case "$dbpath" in

    *SILVA_132*)
        echo " - classify with SILVA 132"

        singularity run "$sing_container" \
            micca classify \
            -m cons \
            -i "$inpdir/otus.fasta" \
            -o "$outdir/taxa-SILVA-132.txt" \
            --ref "$dbpath/rep_set/rep_set_16S_only/97/silva_132_97_16S.fna" \
            --ref-tax "$dbpath/taxonomy/16S_only/97/taxonomy_7_levels.txt"
        ;;

    *Silva_138*)
        echo " - classify with SILVA 138"

        singularity run "$sing_container" \
            micca classify \
            -m cons \
            -i "$inpdir/otus.fasta" \
            -o "$outdir/taxa-SILVA-138.txt" \
            --ref "$dbpath/rep_set/SILVA_138.2_SSURef_NR99_tax_silva.fasta" \
            --ref-tax "$dbpath/taxonomy/taxonomy_7_levels.txt"
        ;;

    *gg_*)
        echo " - classify with Greengenes"

        singularity run "$sing_container" \
            micca classify \
            -m cons \
            -i "$inpdir/otus.fasta" \
            -o "$outdir/taxa_gg.txt" \
            --ref "$dbpath/rep_set/97_otus.fasta" \
            --ref-tax "$dbpath/taxonomy/97_otu_taxonomy.txt"
        ;;

    *)
        echo " - classify with RDP"

        singularity run "$sing_container" \
            micca classify \
            -m rdp \
            -i "$inpdir/otus.fasta" \
            --rdp-gene "$dbpath" \
            -o "$outdir/taxa_rdp.txt"
        ;;
esac

