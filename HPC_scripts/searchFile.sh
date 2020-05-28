#!/bin/bash

#PBS -l nodes=1:ppn=8,walltime=2:00:00
#PBS -q old
#PBS -N searchFile

cd /home/Estancliffe/decoID/HPC_scripts/

export PATH=/export/Anaconda3-5.2.0/bin:$PATH
source activate py3

dataDir=$1   #/scratch/Estancliffe/Asp-Mal_1uM_3Da.mzML
peakfile=$2
db=$3

ms1=$4
dda=$5
iso=$6
unk=$7
lam=$8
ma=$9

python searchFile.py $db $dataDir 8 $peakfile $ms1 $dda $iso $lam $unk $ma