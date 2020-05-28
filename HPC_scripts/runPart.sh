#!/bin/bash

#PBS -l nodes=1:ppn=8,walltime=5:00:00
#PBS -q old
#PBS -N runPart

cd /home/Estancliffe/decoID/HPC_scripts/

export PATH=/export/Anaconda3-5.2.0/bin:$PATH
source activate py3


dataDir=$1   #/scratch/Estancliffe/Asp-Mal_1uM_3Da.mzML


python runPart.py $dataDir $2 8 $3