#!/bin/sh

## setting the enviornmnent
currpath=$(pwd)
project_home="$HOME/bontempo_pigs_caecum_gastroherb"
analysis_dir="${project_home}/Analysis/micca"
outdir="${analysis_dir}/join"
inputdir="${analysis_dir}/trimmed"
sing_container="$HOME/software/micca_latest.sif"
core=8

r1="R1"
r2="R2"

if [ ! -d "${outdir}" ]; then
	mkdir -p ${outdir}
fi


# Count trimmed data
cd $inputdir

echo " - counting sequences "
for i in *.fastq.gz
do
        echo -n $i >> seq_count_16S_QC.txt
        echo -n " " >> seq_count_16S_QC.txt
        echo $(zcat $i | wc -l) / 4 | bc >> seq_count_16S_QC.txt
done


## Join reads (MICCA)

cd ${inputdir}

# remove singles reads from sickle

if ls *singles.gz 1> /dev/null 2>&1; then
	rm *singles.gz
fi

echo " - uncompressing trimmed fastq files "
gunzip *.fastq.gz

## SINGULARITY CONTAINER ##
echo " - joining reads"
## -l: minimum overlap length (bps, default = 32)
## -d max n. of allowed mismatches in the overlap region (default = 8)
## -t: n. of threads to use (1 to 256)
singularity run $sing_container micca mergepairs -i ${inputdir}/*_${r1}.fastq -o ${outdir}/assembled_16S.fastq -p _${r1} -e _${r2} -l 32 -d 8 -t 7

# -l : minimum overlap between reads
# -d : maximum mismatch in overlap region

# Counting reads in assembled file
echo " - counting reads after joining "
grep -c '^@M' ${outdir}/assembled_16S.fastq

# 11534421. Total sum of trimmed reads = 26670364; Teoric 100% assembly = 26670364/2 = 13335182
# Read loss from QC reads to assembled = 13335182 - 11534421 = 1800761;
# Read loss from QC reads to assembled % = 1800761/13335182*100 = 13%

# zipping back trimmed files
echo " - recompressing trimmed fastq files "
cd ${inputdir}

gzip *.fastq

echo "DONE!!"

