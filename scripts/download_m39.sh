#!/bin/bash -l

#SBATCH --partition=angsd_class 
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=download_genome
#SBATCH --time=01:00:00 # HH/MM/SS
#SBATCH --mem=50G # memory requested, units available: K,M,G,T


cwid="xil4009" # specify cornell user id 

mamba activate angsd  # activate environment to access certain tools

cd /athena/angsd/scratch/${cwid}/ # store any tmp files

# genome folder 
GENOMEDIR="/athena/angsd/scratch/${cwid}/genomes"

if [[ ! -d $GENOMEDIR ]]; then
    mkdir -p $GENOMEDIR
fi

# download the "Genome sequence, primary assembly (GRCm39)" fasta file
if [[ ! -f "${GENOMEDIR}/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz" ]]; then
    wget -P ${GENOMEDIR} "ftp://ftp.ensembl.org/pub/release-109/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz"
    echo "m39 genome downloaded"
    gunzip "${GENOMEDIR}/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz"
    echo "gzipped genome" 
fi 

# Note: wget -P means "directory-prefix"


# download the annotations that correspond to it 
if [[ ! -f "${GENOMEDIR}/Mus_musculus.GRCm39.109.gtf.gz" ]]; then
    wget -P ${GENOMEDIR} "ftp://ftp.ensembl.org/pub/release-109/gtf/mus_musculus/Mus_musculus.GRCm39.109.gtf.gz"
    echo "m39 annotation downloaded"
    gunzip "${GENOMEDIR}/Mus_musculus.GRCm39.109.gtf.gz"
    echo "gzipped annotation" 
fi

echo "--------------- Script termination ----------------"