#!/bin/sh
#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=128:mem=512gb
#PBS -N spls_calibration

# PBS job running one of the parallel sPLS calibration scripts. Submit from the
# project root and pass the script to run, e.g.:
#   qsub -v SCRIPT=Scripts/04-sPLS/12-sPLS_perf_centroid.R Scripts/04-sPLS/spls.sh
#   qsub -v SCRIPT=Scripts/04-sPLS/13-sPLS_cor_centroid.R Scripts/04-sPLS/spls.sh

cd "$PBS_O_WORKDIR"

module load anaconda3/personal
source activate r_41

Rscript "$SCRIPT"
