#!/bin/bash -l

#SBATCH --partition=angsd_class 
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=alignqc
#SBATCH --time=08:00:00 # HH/MM/SS
#SBATCH --mem=40G # memory requested, units available: K,M,G,T

#SBATCH --mail-type=ALL  # BEGIN, END, FAIL, or ALL
#SBATCH --mail-user=<xil4009@med.cornell.edu>

# Usage: run_alignqc.sh

cwid="xil4009" # specify cornell user id 

SMG6_project_dir="/athena/angsd/scratch/${cwid}/SMG6_RNA_SEQ" 
alignment_dir="${SMG6_project_dir}/alignment_star"
annotation="/athena/angsd/scratch/${cwid}/genomes/Mus_musculus.GRCm39.109.gtf"

alignment_qc_dir="${SMG6_project_dir}/alignment_qc"

if [[ -d $alignment_qc_dir ]]; then
    echo "deleting old alignment_qc_dir"
    rm -r $alignment_qc_dir
fi

if [[ ! -d $alignment_qc_dir ]]; then
    echo "making new alignment_qc_dir"
    mkdir -p $alignment_qc_dir
fi


cd /athena/angsd/scratch/${cwid}/ # store any tmp files

mamba activate qorts 

for file in "$alignment_dir"/*.bam; do 
    ID=$(echo "$file" | egrep -o 'SRR[0-9]{3,}')
    echo "---------------AlignQC ${ID}----------------"
    # echo $file
    # echo $annotation
    # echo $alignment_qc_dir
    qorts -Xmx20g QC \
        --generatePlots \
        --singleEnded \
        --stranded \
        $file \
        ${annotation}  \
        ${alignment_qc_dir}/$ID

    echo " "
done


echo "--------------- Script Termination ----------------"
exit