#!/bin/bash -l

#SBATCH --partition=angsd_class 
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=feature
#SBATCH --time=03:00:00 # HH/MM/SS
#SBATCH --mem=30G # memory requested, units available: K,M,G,T

# Usage: run_feature_counts.sh 

log_name="run_feature_counts.log"
echo "Starting at:" `date` >> $log_name

cwid="xil4009" # specify cornell user id 

SMG6_project_dir="/athena/angsd/scratch/${cwid}/SMG6_RNA_SEQ" 
alignment_dir="${SMG6_project_dir}/alignment_star"
annotation="/athena/angsd/scratch/${cwid}/genomes/Mus_musculus.GRCm39.109.gtf"

feature_counts_file="${SMG6_project_dir}/feature_counts/featureCounts.txt"


cd /athena/angsd/scratch/${cwid}/ # store any tmp files
mamba activate angsd 


if [ -r $feature_counts_file ]; then
    echo "deleting old feature counts folder ${feature_counts_file%/*}"
    rm -r "${feature_counts_file%/*}"
fi

if [ ! -r $feature_counts_file ]; then
    echo "making new feature counts file $feature_counts_file"
    mkdir -p "${feature_counts_file%/*}" && touch "$feature_counts_file"
fi


# run featureCounts

bam_files=($(find "$alignment_dir" -name "*.bam" -type f))

echo $feature_counts_file

featureCounts \
-a ${annotation} \
-o ${feature_counts_file} \
"${bam_files[@]}"


echo "--------------- Script Termination ----------------"

echo "This is job #:" $SLURM_JOB_ID >> $log_name

exit


