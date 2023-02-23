#!/bin/bash -l

#SBATCH --partition=angsd_class 
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=fastqc
#SBATCH --time=03:00:00 # HH/MM/SS
#SBATCH --mem=30G # memory requested, units available: K,M,G,T

# Usage: fastq2bam.sh <fastq_dir > 

# read input 
fastq_dir=$1

cwid="xil4009" # specify cornell user id 

fastqc_dir="${fastq_dir}/fastqc"

echo "the fastqc dir will be created within fastq dir $fastq_dir"

if [[ ! -d $fastqc_dir ]]; then
    echo "making fastqc dir $fastqc_dir" 
    mkdir -p $fastqc_dir
fi


cd /athena/angsd/scratch/${cwid}/ # store any tmp files
mamba activate angsd  # activate environment to access certain tools

for file in "$fastq_dir"/*q.gz; do 

    ID=$(echo "$file" | egrep -o 'SRR[0-9]{3,}')

    # Run FastQC, if not already present
    if [ ! -r ${fastqc_dir}/${ID}_fastqc ]; then
        echo "---------------FastQC ${ID}----------------"
        fastqc ${file} --extract --outdir $fastqc_dir
    fi

    echo " "
done

echo "--------------- Script Termination ----------------"
exit