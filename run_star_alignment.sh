#!/bin/bash -l

#SBATCH --partition=angsd_class 
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=star_alignment
#SBATCH --time=08:00:00 # HH/MM/SS
#SBATCH --mem=50G # memory requested, units available: K,M,G,T

# Usage: fastq2bam.sh <fastq_dir >

# read input 
fastq_dir=$1 # can align on trimmed or raw fastq

cwid="xil4009" # specify cornell user id 

SMG6_project_dir="/athena/angsd/scratch/${cwid}/SMG6_RNA_SEQ" 
GENOMEDIR="/athena/angsd/scratch/${cwid}/genomes"
starIndex_dir="${GENOMEDIR}/m39_STARindex"

# alignment output directory 
alignment_dir="${SMG6_project_dir}/alignment_star"


if [[ ! -d $alignment_dir ]]; then
    echo "making new alignment dir"
    mkdir -p $alignment_dir
fi


cd /athena/angsd/scratch/${cwid}/ # store any tmp files

mamba activate angsd  # activate environment to access certain tools

for file in "$fastq_dir"/*q.gz; do 

    ID=$(echo "$file" | egrep -o 'SRR[0-9]{3,}')

    # run alignment 
    if [ ! -r "${alignment_dir}/${ID}.Aligned.sortedByCoord.out.bam" ]; then

        echo "---------------STAR alining ${ID}----------------"
        # echo "${alignment_dir}/${ID}"
        STAR --runMode alignReads \
            --runThreadN 6 \
            --genomeDir "${starIndex_dir}" \
            --readFilesIn "${file}" \
            --readFilesCommand zcat \
            --outFileNamePrefix "${alignment_dir}/${ID}". \
            --outSAMtype BAM SortedByCoordinate # output in bam format 
        
        echo "---------------Indexing Bam ${ID}----------------"
        # Index the BAM file
        samtools index ${alignment_dir}/${ID}.Aligned.sortedByCoord.out.bam
    fi

    echo " "
done

echo "--------------- Script Termination ----------------"
exit