#!/bin/sh

## script that filters assembled reads based on quality
## you can tune parameters in the code below
## e.g. discard reads with more than the specified number of allowed Ns (uncalled bases) 
## or discard sequences with more than the specified expected error rate % (e.g. <=1% is less or equal than one error per 100 bases)

set -x

## setting the environment
currpath=$(pwd)
project_home="$HOME/suini_insetti"
analysis_dir="${project_home}/Analysis/micca"
inpdir="${analysis_dir}/join"
outdir="${analysis_dir}/clean"
sing_container="$HOME/software/micca_latest.sif"
core=8

if [ ! -d ${outdir} ]; then
	mkdir -p ${outdir}
fi

# Remove N from assembly
echo " - filtering reads"
## max N's; max error rate (default: 1 error / 100 bps)
singularity run $sing_container micca filter -i ${inpdir}/assembled_16S.fastq -o ${outdir}/assembled_16S_clean.fasta --maxns 0 --maxeerate 1
## --maxns: discard reads with more than the specified number of allowed Ns (uncalled bases) 
## --maxeerate: discard sequences with more than the specified expected error rate % (e.g. <=1% is less or equal than one error per 100 bases: same as Phred score).

# count
echo " - countig reads after filtering"
grep -c '>' ${outdir}/assembled_16S_clean.fasta

echo "DONE!!"

