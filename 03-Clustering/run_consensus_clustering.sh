#!/bin/sh
#PBS -l walltime=10:00:00
#PBS -l select=1:ncpus=1:mem=10gb
#PBS -N clustering
#PBS -q med-bio
#PBS -J 1-300

# PBS array job running 09-Clustering_metab_hpc.R, one job per candidate
# number of clusters. Submit from the project root:
#   qsub Scripts/03-Clustering/run_consensus_clustering.sh

cd "$PBS_O_WORKDIR"

module load anaconda3/personal
source activate r_41

Rscript Scripts/03-Clustering/09-Clustering_metab_hpc.R "$PBS_ARRAY_INDEX"
