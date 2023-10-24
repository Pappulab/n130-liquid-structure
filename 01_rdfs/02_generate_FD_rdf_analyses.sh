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
         repNum="${j}"
         

    sed -e "s@##SYSNAME##@${sysName}@g" \
        -e "s@##REPSNUM##@${repNum}@g"  \
        -e "s@##RUNSDIR##@${trajDir}@g" \
        -e "s@##TOPODIR##@${topoDir}@g" \
        rdf_template_FD.tk > \
         "${sysName}/Run_${repNum}/rdf_FD.tk"
         
    sed -e "s@##SYSNAME##@${sysName}@g" \
        -e "s@##REPSNUM##@${repNum}@g"  \
        -e "s@##RUNSDIR##@${trajDir}@g" \
        -e "s@##TOPODIR##@${topoDir}@g" \
        run_rdf_template_FD.sh > \
         "${sysName}/Run_${repNum}/run_gofr_FD.sh"
         
         
    done
done
    
