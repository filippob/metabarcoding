
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

# count
echo " - countig reads after filtering"
grep -c '>' ${outdir}/assembled_16S_clean.fasta

echo "DONE!!"

