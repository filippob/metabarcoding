#!/bin/bash
#SBATCH --job-name=fastqc    		# Job name
#SBATCH --get-user-env
#SBATCH --partition=long
#SBATCH --nodes=1
#SBATCH --ntasks=1                    # Run on a single task
#SBATCH --cpus-per-task=1             # Declare 1 CPUs per task
#SBATCH --mem=16gb                     # Job memory request
#SBATCH --output=/home/biscarinif/bontempo_pigs_rectum/log/fastqc_%j.log   # Standard output and error log

repofolder="$HOME/metabarcoding_pipeline"
start=`date +%s`
bash $repofolder/0.fastqc.sh
end=`date +%s`
runtime=$((end-start))

echo "execution time was $runtime"
echo "DONE!"

