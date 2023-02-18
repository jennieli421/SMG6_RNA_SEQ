#!/bin/bash -l

#SBATCH --partition=angsd_class 
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=starindex
#SBATCH --time=02:00:00 # HH/MM/SS
#SBATCH --mem=10G # memory requested, units available: K,M,G,T

# Usage:


home_dir="/home/xil4009"
repo_root="${home_dir}/project_scripts"
echo $repo_root

cwid="xil4009" # specify cornell user id 

SMG6_project_dir="/athena/angsd/scratch/${cwid}/SMG6_RNA_SEQ" 


mamba activate angsd  # activate environment to access certain tools

cd /athena/angsd/scratch/${cwid}/ # store any tmp files


# after genome has been generated
GENOMEDIR="${SMG6_project_dir}/genome/"
m39_fa="${GENOMEDIR}/Mus_musculus.GRCm39.dna.primary_assembly.fa"
m39_gtf="${GENOMEDIR}/Mus_musculus.GRCm39.109.gtf"

m39_STARindex="$GENOMEDIR/m39_STARindex"
mkdir -p $m39_STARindex

echo "--------------- Generating Index for Ref Genome ----------------"

# Genome index is created using the command: 
STAR --runMode genomeGenerate \
    --runThreadN 1 \
    --genomeDir "$m39_STARindex" \
    --genomeFastaFiles "$m39_fa" \
    --sjdbGTFfile "$m39_gtf" \
    --sjdbOverhang 99


echo "--------------- Script Termination ----------------"