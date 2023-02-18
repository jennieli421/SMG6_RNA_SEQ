#!/bin/bash -l

#SBATCH --partition=angsd_class 
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=hi_slurm
#SBATCH --time=01:30:00 # HH/MM/SS
#SBATCH --mem=10G # memory requested, units available: K,M,G,T

# Usage: fastq2bam.sh <fastq_dir > <fastqc_dir > <alignment_dir >



home_dir="~xil4009"
repo_root=$(dirname $home_dir)
echo $repo_root


raw_data_dir="/athena/angsd/scratch/xil4009/SMG6_raw_data" 
genome_dir="/home/xil4009/angsd_exercises/exercise4/referenceGenomes/sacCer3_STARindex"


# Check that we have our command line argument(s)
arg_count=$#
if [ $arg_count -lt 3 ]; then
    echo "Not enough command line arguments. Exiting ..."
    echo "Usage: fastq2bam.bash <fastq_dir > <fastqc_dir > <alignment_dir >" exit
fi

# Read arguments from command line
fastq_dir=$1
fastqc_dir=$2
alignment_dir=$3

if [[ -d $alignment_dir ]]; then
    echo "removing old alignment dir!!!"
    rm -r $alignment_dir
fi

if [[ ! -d $alignment_dir ]]; then
    echo "making new alignment dir!!!"
    mkdir -p $alignment_dir
fi



mamba activate angsd


for file in "$fastq_dir"/*.fastq.gz; do 

    ID=$(echo "$file" | egrep -o 'ERR[0-9]{3,}')
    echo "$file"

    echo "---------------FastQC ${ID}----------------"

    # Run FastQC, if not already present
    if [ ! -d ${fastqc_dir}/${ID}_fastqc ]; then
        fastqc ${fastq_dir}/${ID}.fastq.gz --extract --outdir $fastqc_dir
    fi

    

    echo "---------------STAR alining ${ID}----------------"

    # run alignment 
    if [ ! -r "${alignment_dir}/${ID}.Aligned.sortedByCoord.out.bam" ]; then
        STAR --runMode alignReads \
            --runThreadN 1 \
            --genomeDir "${genome_dir}" \
            --readFilesIn "${file}" \
            --readFilesCommand zcat \
            --outFileNamePrefix "${alignment_dir}/${ID}". \
            --outSAMtype BAM SortedByCoordinate 
    fi

    echo "---------------Indexing Bam ${ID}----------------"
    # Index the BAM file
    samtools index ${alignment_dir}/${ID}.Aligned.sortedByCoord.out.bam

    
    echo "---------------Converting Sam to Bam ${ID}----------------"
    # convert BWA SAM to a BAM file 
    samtools view -b ${ID}.bwa.sam -o ${ID}.bwa.bam
    # remove the BWA SAM file 
    rm ${ID}.bwa.sam
    # sort the BAM file 
    samtools sort ${ID}.bwa.bam -o ${ID}.bwa.sorted. bam

    echo " "

done

echo "--------------- Script Termination ----------------"

exit