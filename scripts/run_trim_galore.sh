#!/bin/bash -l

#SBATCH --partition=angsd_class 
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=trim_galore
#SBATCH --time=05:00:00 # HH/MM/SS
#SBATCH --mem=50G # memory requested, units available: K,M,G,T

# Usage: run_trim_galore.sh <fastq_dir > 

# read input 
fastq_dir=$1

cwid="xil4009" # specify cornell user id 

trim_dir="${fastq_dir}/../fastq_trimmed"

if [[ ! -d $trim_dir ]]; then
    echo "making fastqc dir $trim_dir" 
    mkdir -p $trim_dir
fi


cd /athena/angsd/scratch/${cwid}/ # store any tmp files

mamba activate angsd  # activate environment

for file in "$fastq_dir"/*.fastq.gz; do 

    ID=$(echo "$file" | egrep -o 'SRR[0-9]{3,}')

    # Run trim-galore, if not already present
    if [ ! -r ${trim_dir}/${ID}_trimmed.fq.gz ]; then
        echo "---------------TrimGlaore ${ID}----------------"
        mamba activate trim-galore
        trim_galore -o $trim_dir --gzip --illumina $file 
    fi

    echo " "
done

echo "--------------- Script Termination ----------------"
exit