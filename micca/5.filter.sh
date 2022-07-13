
## setting the environment
currpath=$(pwd)
project_home="$HOME/useful_microbiome_chlorine"
outdir="${project_home}/Analysis/micca"
sing_container="$HOME/software/micca.sif"
core=8

if [ ! -d "${outdir}/micca_16S" ]; then
	mkdir -p ${outdir}/micca_16S
fi

# Remove N from assembly
echo " - filtering reads"
singularity run $sing_container micca filter -i ${outdir}/micca_16S/WP1_assembled_16S.fastq -o ${outdir}/micca_16S/WP1_assembled_16S.fasta --maxns 0

# count
echo " - countig reads after filtering"
grep -c '>' ${outdir}/micca_16S/WP1_assembled_16S.fasta

echo "DONE!!"

