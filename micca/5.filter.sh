
## setting the environment
currpath=$(pwd)
project_home="$HOME/bontempo_pigs_rectum"
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
singularity run $sing_container micca filter -i ${inpdir}/assembled_16S.fastq -o ${outdir}/assembled_16S_clean.fasta --maxns 0

# count
echo " - countig reads after filtering"
grep -c '>' ${outdir}/assembled_16S_clean.fasta

echo "DONE!!"

