#!/bin/bash

directory=$1
peakfile=$2
db=$3

ms1=$4
dda=$5
iso=$6
unk=$7
lam=$8
ma=$9
lab=${10}
frag_cutoff=${11}
rtTol=${12}

for filename in $directory*.mzML; do
    qsub -F "$filename $peakfile $db $ms1 $dda $iso $unk $lam $ma $lab $frag_cutoff $rtTol" /home/Estancliffe/DecoID/searchFile.sh
done
