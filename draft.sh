#!/bin/bash

# start by getting the absolute path to the directory this script is in, which will be the top level of the repo
# this way script will work even if the repo is downloaded to a new location, rather than relying on hard coded paths to where I put the repo. 
full_path=$(realpath $0)
repo_root=$(dirname $full_path)

# export the path to the repo for scripts called by this script to also use - will unset at end
export repo_root

# gather user settings, first asking which study the code should run on - this is only setting currently for the viz side
# NOTE: take command line argument as subject name variable 
study=$1
export study



data_dir="WT_2" # the dir to save downloaded fastq files

if [[ ! -d $data_dir ]]; then
	mkdir -p $data_dir # create output folder if there isn't already
fi


# Remove old metadata file if exists. Avoid concatonation error. 
if [[ -f $TM_folder/$subject_level_text_loc/FRESH_17_subject_level_text.csv ]]; then
	echo "Removing old subject-level text summary"
	rm $TM_folder/$subject_level_text_loc/FRESH_17_subject_level_text.csv
fi



for file in "$fastq_dir"/*.fastq.gz; do 

    # if file already exist, skip download and continue 
    if [[ -f "$data_dir/$ID.fastq.gz" ]]; then
        echo "$ID fastq already exist, skipping"
        continue
    fi 


    echo "$ID"
    link_add="ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR458/${ID}/${ID}.fastq.gz"
    wget -P $data_dir $link_add

    echo " "

    # run alignment 
    if [ ! -r "${alignment_dir}/${sample}_${ID}.Aligned.sortedByCoord.out.bam" ]; then
        STAR --runMode alignReads \
            --runThreadN 1 \
            --genomeDir "${genome_dir}" \
            --readFilesIn "${file}" \
            --readFilesCommand zcat \
            --outFileNamePrefix "${alignment_dir}/${sample}_${ID}". \
            --outSAMtype BAM SortedByCoordinate 
    fi
    
    # if STAR exist but there is no log, means it didn't finish properly

done