#!/bin/bash

thisDir=$(pwd)
runsDir="XX_XX"
topoDir="XX_XX"

trajDir="${runsDir}"

sysNames=("lk_8" "L_10" "L_16" "L_20" "no_A0" "no_A1" "no_A2")

for sysName in "${sysNames[@]}"
do
    for repNum in 1 2 3 4 5
    do
        sed -e "s@##SYSNAME##@${sysName}@g" \
            -e "s@##REPSNUM##@${repNum}@g"  \
            -e "s@##RUNSDIR##@${trajDir}@g" \
            -e "s@##TOPODIR##@${topoDir}@g" \
            network_template_charge.py > \
             "${sysName}/Run_${repNum}/network_charge.py"
             done
done
     
