#!/bin/bash


#PBS -l nodes=1:ppn=8,walltime=0:30:00
#PBS -q old
#PBS -N testDeco

cd /home/Estancliffe/decoID/specificApplicationScripts/

export PATH=/export/Anaconda3-5.2.0/bin:$PATH
source activate py3


dataDir=$1

python searchDataFileCMD.py $dataDir 1