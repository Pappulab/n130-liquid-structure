#!/bin/bash

thisDir=$(pwd)
runsDir="XX_XX"
topoDir="XX_XX"

sysNames=(
"lk_8"
"L_10"
"L_16"
"L_20"
"no_A0"
"no_A1"
"no_A2"
)

for sysName in "${sysNames[@]}"
do
    for j in 1 2 3 4 5
    do
        mkdir -p "${thisDir}/${sysName}/Run_${j}/"
    done
done
    
