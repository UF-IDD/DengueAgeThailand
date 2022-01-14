#!/usr/bin/env Rscript --vanilla

library(argparse)
library(tidyverse)
library(rstan)

parser = ArgumentParser(description = "Fit stan model for a province.")

parser$add_argument('model.name')
parser$add_argument('place')
parser$add_argument('-datadir', default = "04-dataForStan")
parser$add_argument('-casever', required = T)
parser$add_argument('-stanbase', default = '05-estimates_16tau')
parser$add_argument('-stanfile', default = 'Scripts/stan/piecewise_16tau.stan')

parser$add_argument('-ncpu', type = "integer", default = 3)
parser$add_argument('-n_iter', type = "integer", default = 15000)
parser$add_argument('-n_warmup', type = "integer", default = 3000)


# parse arguments
inputArg = parser$parse_args()
attach(inputArg)


# prepare output directory
outdir = file.path(stanbase,casever,model.name)
dir.create(outdir, recursive=T)


# define some functions
fun = new.env()
fun$readDataCsv = function(x){
    as.matrix(read.csv(x, stringsAsFactors=F))
}
fun$getModelSetup = function(model.name){
    out = list()
    out$useQe   = grepl('^m[0-9]h',model.name)
    out$useQeT = grepl('QeT$',model.name)
    out$useQeA = grepl('QeA$',model.name)
    out$useK    = grepl('^m[2-3]',model.name)
    out$usePhiT = grepl('\\.rt',model.name)
    out$usePhiA = grepl('\\.r[t]*a',model.name)
    out$usePhiAT = grepl('\\.rAT',model.name)
    out
}



# import data and their configurations
case = fun$readDataCsv(paste0(datadir,"/",casever,"/",place,".csv"))
pop = fun$readDataCsv(paste0(datadir,"/pop/",place,".csv"))
config = readRDS( paste0(datadir,"/config.RDS") )




# get model setup
modelSetup = fun$getModelSetup(model.name)

# make stan input
stanInput = list()
stanInput$len_cohorts = length(config$cohorts)
stanInput$len_times = length(config$times)
stanInput$len_ages = length(config$ages)
stanInput$len_Ag = nrow(config$ageGroups)

# for PhiAT
nPhiAgeBins = 8
nPhiTimeBins = 4
# for all others
nAgeBins = nTimeBins = 16


# stanInput$len_binnedQeT = 1 + (4-1) * modelSetup$useQeT
stanInput$len_binnedQeT = 4 * modelSetup$useQeT # 4 epoches
stanInput$binwidth_QeT = 9

stanInput$len_binnedQeA = 4 * modelSetup$useQeA # 4 agebins
stanInput$binwidth_QeA = 15

stanInput$len_binnedK = nAgeBins * modelSetup$useK
stanInput$binwidth_K = 4

stanInput$len_binnedPhiT = nTimeBins * modelSetup$usePhiT
stanInput$binwidth_PhiT = 2

stanInput$len_binnedPhiA = nAgeBins * modelSetup$usePhiA
stanInput$binwidth_PhiA = 4

stanInput$dim_binnedPhiAT = c(nPhiAgeBins,nPhiTimeBins) * modelSetup$usePhiAT
stanInput$binwidth_PhiAT = c(8,9)



    #   Compile and Run
    #   ...............


file.copy(stanfile, outdir, overwrite=T)

compiled = stan(
    file = paste0(outdir,"/",basename(stanfile))
    ,model_name = model.name
    ,data = c( stanInput, list(
        naVal = config$naVal
        ,mapAge = config$ageGroups
        ,countObserved = case
        ,countPop = pop
    ))
    ,chains = 0
)



initVals = function(...){
    list(
        omega = runif(1, 1, 2)
        ,theta = runif(1, 1, 2)
        ,lambda = runif(16, 0, .6)
        ,binnedK = runif( stanInput$len_binnedK, 1, 1)
        ,binnedPhiT = runif( stanInput$len_binnedPhiT, .2, .8)
        ,binnedPhiA = runif( stanInput$len_binnedPhiA, .2, .8)
        ,binnedPhiAT = matrix( runif( prod(stanInput$dim_binnedPhiAT), .2, .8 )
            ,stanInput$dim_binnedPhiAT[1]
            ,stanInput$dim_binnedPhiAT[2]
        )
        ,binnedQe = matrix(
            runif( with(stanInput, max(len_binnedQeT,len_binnedQeA,1)) * 4
                ,c(0.05,.8,0.03,0.03)
                ,c(.1,.95,.05,.05)
            )
            , 4
            , with(stanInput, max(len_binnedQeT,len_binnedQeA,1))
        ) %>% t
    )
}



    

stanRDS = paste0(outdir,"/",place,".RDS")
if(!file.exists(stanRDS)){
    set.seed(ncpu*2)
    message(paste("Initializing", ncpu, "chains..."))
    fit = stan( fit=compiled
        ,data = c( stanInput, list(
            naVal = config$naVal
            ,mapAge = config$ageGroups
            ,countObserved = case
            ,countPop = pop
        ))
        ,iter = n_iter
        ,warmup = n_warmup
        ,chains = ncpu
        ,cores = ncpu
        ,seed = runif(1)*10^6
        ,refresh = 20
        ,thin = 5
        ,init = lapply(1:ncpu, initVals)
        ,include = FALSE
        ,pars = c(
            "lambda_times_expanded"
            ,"ageCountExpected"
            ,"countExpected"
            ,"nEntries"
            ,"Qe_t"
            ,"PhiA"
            ,"PhiT"
            ,"PhiAT"
            ,"K"
        )
    )
    saveRDS( fit
        ,file = stanRDS
    )
} else {
    message(paste("Fit exists:",stanRDS))
    message("Exiting...")
}