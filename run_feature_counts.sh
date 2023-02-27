#!/bin/bash -l

#SBATCH --partition=angsd_class 
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=feature
#SBATCH --time=03:00:00 # HH/MM/SS
#SBATCH --mem=30G # memory requested, units available: K,M,G,T

# Usage: run_feature_counts.sh 


cwid="xil4009" # specify cornell user id 

SMG6_project_dir="/athena/angsd/scratch/${cwid}/SMG6_RNA_SEQ" 
alignment_dir="${SMG6_project_dir}/alignment_star"
annotation="/athena/angsd/scratch/${cwid}/genomes/Mus_musculus.GRCm39.109.gtf"

feature_counts_path="${SMG6_project_dir}/feature_counts/featureCounts.txt"


cd /athena/angsd/scratch/${cwid}/ # store any tmp files
mamba activate angsd 


# run feature count 
if [ ! -r "${feature_counts_path}" ]; then

    bam_files=($(find "$alignment_dir" -name "*.bam" -type f))

    featureCounts \
    -a ${annotation} \
    -o ${feature_counts_path} \
    "${bam_files[@]}"

fi

echo "--------------- Script Termination ----------------"
exit


