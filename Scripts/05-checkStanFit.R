#!/usr/bin/env Rscript --vanilla

library(argparse)
library(tidyverse)
library(rstan)
#rstan_options(auto_write = TRUE)

parser = ArgumentParser(description = "Check stan fits and compare.")


parser$add_argument('dataver')
parser$add_argument('-stanroot', default = "05-estimates")
parser$add_argument('-od', default = "06-assessFits")
parser$add_argument('-ncpu', type = "integer", default = 3)

# parse arguments
inputArg = parser$parse_args()


# function to set plot directory
outdir = with(inputArg, file.path(od, dataver))
setPlotdir = function(...){
    plotdir <<- file.path(outdir, ...)
    message(paste("Set plot directory:",plotdir))
    dir.create(plotdir, recursive=T)
}
# create root for outputs
setPlotdir()




    #   Plot fit file existence
    #   .......................
    
x = list.files( with(inputArg, file.path(stanroot, dataver))
        , pattern = "\\.RDS$"
        , recursive=T
    ) %>%
    strsplit("/") %>%
    do.call(what=rbind) %>%
    as.data.frame
    colnames(x) = c('model','province')
x$province = gsub('\\.RDS$','',x$province)
x = x %>% mutate(avail = 1)

g = x %>%
    ggplot(aes(
        x = model
        ,y = province
        ,fill = avail
    ))+
    geom_tile()+
    guides(fill = F)+
    theme(
        axis.text.x = element_text(angle=45, hjust=1)
        ,panel.grid = element_blank()
    )
ggsave( g
    ,filename = file.path(outdir,"fitAvail.pdf")
    ,width = 5
    ,height = 10
)




    #   Check convergence of each stan output
    #   .....................................

x$rds = list.files( with(inputArg, file.path(stanroot, dataver))
    , pattern = "\\.RDS$"
    , recursive=T
    , full.names = T
) 
# split by province
x = split(x, x$province)

library(parallel)
cl = makeCluster(inputArg$ncpu, outfile=file.path(plotdir,"clust.out"))


# connection to log corrupted files
out.corrupt = file.path(outdir,"corruptedFiles.txt")
file.remove(out.corrupt)
out.corrupt = file(out.corrupt, open = "a+b")
clusterExport(cl, c(
    "out.corrupt"
    , "outdir"
    , "setPlotdir"
))


# for each province
LooList = lapply(x, function(xProv){

    # compute psis-loo for each model
    provLooList = parApply(cl, X = xProv
        , 1, FUN = function(xProvModel){

            attach(as.list(xProvModel))
            message(rds)
            library(rstan)
            library(dplyr)
            library(ggplot2); theme_set(theme_bw())
            library(ggpubr)
            library(loo)

            fit = try(readRDS(rds))
            if(class(fit)=="try-error"){
                writeLines( rds, con = out.corrupt)
                return(NULL)
            }
            
            
            fitsum.raw = summary(fit)$summary
            if(is.null(fitsum.raw)){ return(NULL) }


            vars = grep('^log_lik|^lp__|^rmse', rownames(fitsum.raw), value=T, invert=T)
        
            getFitSummary = function(fitsum, vars, ess.thres=300){
                fitsum = fitsum[vars, c('mean','2.5%','97.5%','n_eff','Rhat')]
                colnames(fitsum)[2:4] = c('lower','upper','ESS')

                fitsum %>%
                    as.data.frame %>%
                    mutate(
                        var = gsub('\\[[0-9,]+\\]$','',vars)
                        , ivar = as.integer(gsub('[^0-9]+','',vars))
                    ) %>%
                    group_by(var) %>%
                    mutate( ivar = seq_along(var)) %>%
                    ungroup %>%
                    mutate(
                        converged = Rhat < 1.1 & ESS >= ess.thres
                        , ess.thres = ess.thres
                    ) %>%
                    filter(!is.na(converged))
            }
            performance = getFitSummary(fitsum.raw, c('rmse', 'log_lik_sum'))
            fitsum = getFitSummary(fitsum.raw, vars)
        
            # discard more as warmup
            if( !all(fitsum$converged) ){          
                sims = rstan::extract(fit, permuted=F, inc_warmup=F)
                nwarm = floor(dim(sims)[1]/2)
                    message(paste("Discarded as warmup:",nwarm))
                    fitsum.raw = monitor(
                        sims
                        , warmup = nwarm
                        , print = F
                    ) %>% as.matrix
                    performance = getFitSummary(fitsum.raw, c('rmse', 'log_lik_sum'))
                    fitsum = getFitSummary( fitsum.raw, vars)
                nwarm = nwarm * fit@stan_args[[1]]$thin
            } else {
                nwarm = 0 
            }
            fitsum = fitsum %>%
                mutate(nwarm = nwarm + fit@stan_args[[1]]$warmup)
            
            # if still cannot meet convergence criteria,
            # plot the trace for diagnostics
            if( !all(fitsum$converged) ){
                setPlotdir('trace_notConverged/', model)
                g = traceplot(fit, pars = vars[!fitsum$converged]) +
                    theme( axis.text.x = element_blank() )
                ggsave( g
                    ,filename = file.path(plotdir, paste0(province,'.jpg'))
                    ,width = 12
                    ,height = 8
                )
            } else {
            # otherwise, save summary for later use
                setPlotdir('summary_converged/', model)
                write.csv( fitsum
                    ,file = file.path(plotdir, paste0(province,'.csv'))
                    ,row.names = F
                )
            
                # save performance
                setPlotdir('performance_converged/', model)
                write.csv( performance
                    ,file = file.path(plotdir, paste0(province,'.csv'))
                    ,row.names = F
                )
            
                # compute loo object
                ll = extract_log_lik(fit, merge_chains=F)
                exclude = ll[1, 1, ] > 1
                ll = ll[ , ,!exclude]        
                fLoo = loo(
                    ll
                    ,r_eff = relative_eff(exp(ll))
                    ,cores = 1
                )
                # store problematic data points if exists
                setPlotdir('high_paretoK_converged/', model)
                i = sapply(pareto_k_ids(fLoo, threshold = 0.7), function(i){
                    min(which(cumsum(!exclude)==i))
                })
                if(length(i)>0){
                    write.csv(
                        data.frame(
                            t = ceiling(i/17)
                            ,ag = ifelse((i %% 17) == 0, 17, i %% 17)
                        )
                        ,file = file.path(plotdir, paste0(province,'.csv'))
                        ,row.names = F
                    )
                }

                # and return to compare ELPD across models later
                return(fLoo)
            }
        return(NULL)
        }
    )
    names(provLooList) = xProv$rds

    setPlotdir('elpd_converged','loo')
    saveRDS(provLooList, file=paste0(plotdir,"/",xProv$province[1],".RDS"))
    
    # compare LOO between models
    library(loo)
    setPlotdir('elpd_converged')
    provLooList = provLooList[!sapply(provLooList, is.null)]
    if(length(provLooList)==0){
        return(NULL)
    } else if(length(provLooList)<2){
        x = data.frame(
            elpd_diff = 0
            ,se_diff = 0
        )
        rownames(x) = names(provLooList)
    } else {
        x = loo_compare(x = provLooList) %>%
            as.data.frame %>%
            select(elpd_diff, se_diff)
    }
    meta = strsplit(rownames(x),"/") %>%
        do.call(what=rbind)
    x$rank = with(x, (se_diff * 2) <= -elpd_diff ) %>% cumsum
    write.csv( x %>%
        mutate(
            dataver = meta[ ,2]
            ,model = meta[ ,3]
            ,province = xProv$province[1]
        )

        ,file = paste0(plotdir,'/',xProv$province[1],'.csv')
        ,row.names = F
    )
    return(NULL)
})



stopCluster(cl)
close(out.corrupt)
