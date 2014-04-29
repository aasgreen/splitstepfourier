#!/bin/bash

mkdir Runsbmn
cd Runsbmn/
for N in {2..4}
do
    for distance in {1..2}
    do
        mkdir "N$N-D$distance"
        cd  "N$N-D$distance"
        cp ../../*.py ./
        cp ../../paratemplate ./
        python paramgen.py $distance 1 $N 0 0
        mv par-spec.py param.py
        python ssnl2.py
        cd ..
    done
done
