#!/bin/bash -l

#SBATCH --partition=angsd_class 
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=fastq
#SBATCH --time=03:00:00 # HH/MM/SS
#SBATCH --mem=10G # memory requested, units available: K,M,G,T

# Usage:
cwid="xil4009" # specify cornell user id 

home_dir="/home/${cwid}"
repo_root="${home_dir}/project_scripts"

raw_data_dir="/athena/angsd/scratch/${cwid}/SMG6_RNA_SEQ/fastq_raw" # the dir to save downloaded fastq files

if [[ ! -d $raw_data_dir ]]; then
	mkdir -p $raw_data_dir # create output folder if there isn't already
fi


mamba activate angsd  # activate environment to access certain tools
cd /athena/angsd/scratch/${cwid}/ # store any tmp files

echo "--------------- Downloading Fastq ----------------"

runid_ftp="${repo_root}/metadata/filereport_read_run_PRJNA776901_tsv.txt" # cols: run_accession	fastq_ftp

tail -n +2 ${runid_ftp} | while read line; do 

    vars=($line)

    ID=${vars[0]}
    ftp="ftp://${vars[1]}"

    # if file already exist, skip download and continue 
    if [[ -f "$raw_data_dir/$ID.fastq" ]] ||  [[ -f "$raw_data_dir/$ID.fastq.gz" ]]; then
        echo "$ID fastq already exist, skipping"
        continue
    fi 

    echo "downloading $ID"
    echo $ftp

    if [ $ID == "SRR16684379" ]; then
        fasterq-dump -O $raw_data_dir ${ID}
        gzip $raw_data_dir/${ID}.fastq
        echo "Downloaded $ID fastq using dump"
    else 
        wget -P $raw_data_dir $ftp
        echo "Downloaded ${ID} fastq"
    fi 
    
    echo " "
done

echo "--------------- Script termination ----------------"
exit