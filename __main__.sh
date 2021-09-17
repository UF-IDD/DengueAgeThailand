#!/bin/bash

# Describe data
Rscript --vanilla Scripts/04-describeData.R


# Fit stan models

    # define directories
    datadir="04-dataForStan"
    fitdir="05-estimates"
    casever="DHF_add"
    
    # get the list of provinces
    ls ${datadir}/DHF_add | cut -d'.' -f1 > provinces.txt

    # function to run stan
    runStan_province()
    {
        prov=$(cat provinces.txt | sed -n ${1}p )
        for model in "m1h" "m1h.rt" "m1h.ra" "m3h" "m3h.rt" "m3h.ra" "m1h.rta" "m3h.rta" "m1h.rAT" "m3h.rAT"; do
        for qeMode in "" "QeT" ; do
            model=${model}${qeMode}
            if [[ ! -f "${fitdir}/${model}/${prov}.RDS" ]]; then
                Rscript --vanilla Scripts/04-fitStan.R ${model} ${prov} -casever ${casever}
            fi
        done
        done
    }
    for iprov in $(seq 1 72 ) ; do
        runStan_province ${iprov}
    done


# Check convergence of model fits and do model comparison
Rscript --vanilla Scripts/05-checkStanFit.R "DHF_add" -od "06-assessFits"


# Perform inferences, and plot figures
Rscript --vanilla 06-inferenceFromFits.R


# Simulate theoretical population (Figure S8)
Rscript --vanilla Scripts/10-betaDemogTransitionEffects.R
