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



basicTractSels=(
"44 45 46 47 50 56 57 58" 
"46 47 48 49 60 61 62" 
"46 47 48 49 66 67 68" 
"46 47 48 49 70 71 72"
"44 45 46 47 50 56 57 58"
"44 45 46 47 50 56 57 58"
"44 45 46 47 50 56 57 58"
)

trajDir="${runsDir}"

for i in 0 1 2 3 4 5 6
do
    sysName="${sysNames[i]}"
    basicSel="${basicTractSels[i]}"

    for j in 1 2 3 4 5
    do
        repNum="${j}"
        
        sed -e "s@##SYSNAME##@${sysName}@g" \
            -e "s@##REPSNUM##@${repNum}@g"  \
            -e "s@##RUNSDIR##@${trajDir}@g" \
            -e "s@##TOPODIR##@${topoDir}@g" \
            -e "s@##BASIC##@${basicSel}@g"  \
            rdf_template_charged.tk > \
             "${sysName}/Run_${repNum}/rdf_charged.tk"
             
        sed -e "s@##SYSNAME##@${sysName}@g" \
            -e "s@##REPSNUM##@${repNum}@g"  \
            -e "s@##RUNSDIR##@${trajDir}@g" \
            -e "s@##TOPODIR##@${topoDir}@g" \
            run_rdf_template_charged.sh > \
             "${sysName}/Run_${repNum}/run_rdf_charged.sh"
         
    done
done
    
