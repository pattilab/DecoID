#!/bin/bash

filename=$1
numFiles=$2
lam=$3

cd /home/Estancliffe/decoID/HPC_scripts/


tmpName="unknowns"

qsub -F "$1 $tmpName" runPart.sh

for (( tmpName=0; tmpName<$numFiles; tmpName++ ))
do
    qsub -F "$1 $tmpName $lam" runPart.sh
done