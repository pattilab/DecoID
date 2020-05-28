#!/bin/bash

cd /home/Estancliffe/decoID/HPC_scripts/

directory=$1
peakfile=$2
db=$3

ms1=$4
dda=$5
iso=$6
unk=$7
lam=$8
ma=$9

for filename in $directory*.mzML; do
    qsub -F "$filename $peakfile $db $ms1 $dda $iso $unk $lam $ma" ./searchFile.sh
done
