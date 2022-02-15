#!/bin/bash

# Describe data
Rscript --vanilla Scripts/04-describeData.R


# Fit stan models

    # define directories
    datadir="04-dataForStan"
    fitdir="05-estimates"
    fitver=""
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
                Rscript --vanilla Scripts/04-fitStan${fitver}.R ${model} ${prov} -casever ${casever} \
                    -stanbase ${fitdir}${fitver}\
                    -stanfile Scripts/stan/piecewise${fitver}.stan
            fi
        done
        done
    }
    for iprov in $(seq 1 72 ) ; do
        runStan_province ${iprov}
    done


# Check convergence of model fits and do model comparison
checkFit(){
    Rscript --vanilla Scripts/05-checkStanFit.R "DHF_add" \
        -stanroot "${fitdir}${fitver}"
        -od "06-assessFits${fitver}"    
}
checkFit

# Plot priors vs posteriors for the best fit model
Rscript --vanilla Scripts/05-PriorVsPosterior.R DHF_add

# Assess sensitivity of parameter estimates

    # Fit best model using weaker priors on Q(i,t)
    fitver="_widePrior"
    # use same R file as before, but a different stan file which has wider priors
    cp Scripts/04-fitStan.R Scripts/04-fitStan${fitver}.R
    for iprov in $(seq 1 72 ) ; do
        runStan_province ${iprov}
    done
    checkFit
    
    # Refit all models with less degrees of freedom for tau(t)
    fitver="_16tau"
    for iprov in $(seq 1 72 ) ; do
        runStan_province ${iprov}
    done
    checkFit

    # Do the sensitivity assessment
    Rscript --vanilla Scripts/07-sensitivityAnalysis.R


# Perform inferences, and plot figures
Rscript --vanilla 06-inferenceFromFits.R

# Simulate theoretical population (Figure S15)
Rscript --vanilla Scripts/10-betaDemogTransitionEffects.R
