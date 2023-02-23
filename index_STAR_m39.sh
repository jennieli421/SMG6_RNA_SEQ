#!/bin/bash -l

#SBATCH --partition=angsd_class 
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=starindex
#SBATCH --time=03:00:00 # HH/MM/SS
#SBATCH --mem=50G # memory requested, units available: K,M,G,T

# Usage:


cwid="xil4009" # specify cornell user id 

mamba activate angsd  # activate environment to access certain tools
cd /athena/angsd/scratch/${cwid}/ # store any tmp files

# after genome has been generated
GENOMEDIR="/athena/angsd/scratch/${cwid}/genomes"
m39_fa="${GENOMEDIR}/Mus_musculus.GRCm39.dna.primary_assembly.fa"
m39_gtf="${GENOMEDIR}/Mus_musculus.GRCm39.109.gtf"

m39_STARindex="${GENOMEDIR}/m39_STARindex"
echo "making directory $m39_STARindex"
mkdir -p $m39_STARindex

echo "--------------- Generating Index for Ref Genome ----------------"

# Genome index is created using the command: 
STAR --runMode genomeGenerate \
    --runThreadN 6 \
    --genomeDir "${m39_STARindex}" \
    --genomeFastaFiles "${m39_fa}" \
    --sjdbGTFfile "${m39_gtf}" \
    --sjdbOverhang 99

echo "--------------- Script Termination ----------------"