#!/usr/bin/env Rscript --vanilla
rm(list=ls())



fitdir = "06-assessFits"
dataver = "DHF_add"
indir = file.path(fitdir, dataver)

library(tidyverse); theme_set(theme_bw())
library(ggpubr)


outdir = file.path("07-inferenceFromFits", dataver)
setPlotdir = function(...){
    plotdir <<- file.path(outdir, ...)
    message(paste("Set plot directory:",plotdir))
    dir.create(plotdir, recursive=T)
}
setPlotdir()




    #   Import data
    #   ...........

# spatial data
sp = new.env()
source("Scripts/01-extractSpatialData.R", local=sp)



# case, pop data
dat = new.env()    
dat$datadir = "04-dataForStan"
dat$readDataCsv = function(x){
    as.matrix(read.csv(x, stringsAsFactors=F))
}

places = sp$provinceCentroids$province72 %>% unique
names(places) = places
dat$pop = lapply( places, function(place){
    dat$readDataCsv(paste0(dat$datadir,"/pop/",place,".csv"))
})
dat$case = lapply( places, function(place){
    dat$readDataCsv(paste0(dat$datadir,"/",dataver,"/",place,".csv"))
})
config = readRDS( paste0(dat$datadir,"/config.RDS") )
config$ageGroupLabels = apply(config$ageGroups, 1, function(x){
    x = range(config$ages[x==1])
    if(x[1]==x[2]){
        return(x[1])
    } else if(x[2]==99){
        return(paste0(x[1],'+'))
    } else {
        paste0(x[1],'-',x[2])
    }
})



    #   Define some functions
    #   .....................
    fun = new.env()
    
fun$models = c(
    "m1h"   # null model
    ,"m3h", "m1h.ra", "m1h.rt"  # single component models
    ,"m3h.ra", "m3h.rt", "m1h.rta", "m1h.rAT" # two components
    ,"m3h.rta", "m3h.rAT" # three components
)
fun$models = c(
    fun$models
    , paste0(fun$models,"QeT")
)
fun$getModelSetup = function(model.name){
    out = list()
    out$useQe   = grepl('^m[0-9]h',model.name) & ! grepl('QeT$',model.name)
    out$useQeT = grepl('QeT$',model.name)
    out$useK    = grepl('^m[2-3]',model.name)
    out$usePhiT = grepl('\\.rt',model.name)
    out$usePhiA = grepl('\\.r[t]*a',model.name)
    out$usePhiAT = grepl('\\.rAT',model.name)
    out$useLambdaT = !grepl('m0',model.name)
    out
}
fun$modelComponents = c(
    "QeT" = expression(Q(i,t))
    ,"Qe" = expression(Q(i))
    ,"PhiAT" = expression(phi(a,t))
    ,"PhiA" = expression(phi(a))
    ,"PhiT" = expression(phi(t))
    ,"K" = expression(kappa(a))
    ,"LambdaT" = expression(bar(tau)(t))
)
fun$getModelComponent = function(model.name){
    setup = fun$getModelSetup(model.name)
    data.frame(
        model = model.name
        ,parname = names(fun$modelComponents)
        ,usePar = unlist(setup[ paste0("use",names(fun$modelComponents)) ])
    )
}



    #   plot performance comparisons
    #   ............................

# import: RMSE, LL
x = lapply(
    list.files(
        file.path(indir,"performance_converged")
        ,pattern = "\\.csv$"
        ,recursive = T
        ,full.names = T
    )
    ,function(f){
        meta = strsplit(f,"/")[[1]]
        read.csv(f, stringsAsFactors=F) %>%
            mutate(
                model = meta[4]
                ,province = gsub('\\.csv$','',meta[5])
            )
    }
) %>%
do.call(what=rbind)

# double check that all did converge
summary(x)
x %>% filter(ESS < 300, Rhat > 1.01) # filter for non-converging entries, if exists

x = x %>%
    select(province, model, var, mean, lower, upper) %>%
    rename(val = mean)

x = split(x, x$var)

# RMSE of the null model
x$rmse %>%
    filter(model == "m1h") %>%
    summarize(
        Mean = mean(val)
        , Median = median(val)
        , lower = quantile(val, .025)
        , upper = quantile(val, .975)
    ) %>%
    round(2)



# function to determine country-level best model
# based on a given metric
getCountryBestModel = function(x, value_col = val, smaller_better = 1, maxRank = 1){
    x %>%
        mutate(val = {{value_col}} ) %>%
        group_by(province) %>%
        mutate(ranking = rank(smaller_better * val, ties.method = "min")) %>%
        filter(ranking <= maxRank) %>%

        # fractional votes if there are ties
        group_by(province, ranking) %>%
        mutate(vote = 1/n()) %>%       
        
        group_by(model, ranking) %>%
        summarize( count = n() ) %>%
        group_by(ranking) %>%
        filter(count == max(count)) %>%
        ungroup %>%
        arrange(ranking, factor(model, levels = fun$models)) %>%
        select(model) %>%
        head(maxRank) %>%
        unlist
}

# import: ELPD
elpd = lapply(
        list.files(
            file.path(indir,"elpd_converged")
            ,pattern = "\\.csv"
            ,full.names = T
        )
        ,read.csv
        ,stringsAsFactors = F
    ) %>%
    do.call(what=rbind)
    

# Top 3: elpd
getCountryBestModel(elpd
    , value_col=elpd_diff
    , smaller_better = -1
    , maxRank = 3
)

# make things relative to country-level best model
elpd = elpd %>%
    left_join(
        elpd %>%
            filter(model == getCountryBestModel(elpd
                , value_col=elpd_diff
                , smaller_better = -1
            )) %>%
            select(province, elpd_diff)
        ,by = "province"
        ,suffix = c("",".best")
    ) %>%
    mutate(
        # make it the lower the better
        delta = -(elpd_diff - elpd_diff.best)
        , lowerDelta = delta - 2*se_diff
        , upperDelta = delta + 2*se_diff
    ) %>%
    mutate(var = "ELPD") %>%
    group_by(province) %>%

    mutate(
        rankVal = rank(rank, ties.method="min")
        # make the lower the better
        , val = -elpd_diff
        , lower = -elpd_diff - 2*se_diff
        , upper = -elpd_diff + 2*se_diff
    ) %>%
    select(-rank, -elpd_diff, -elpd_diff.best, -se_diff) %>%
    ungroup

# fix "Nan" province issue
elpd$province = with(elpd, ifelse(province=="NaN","Nan",province))

# import estimates
elpd$est = lapply(
    with(elpd, file.path(
        indir,'summary_converged', model
        , paste0(province, '.csv')
    ))
    ,read.csv
    ,stringsAsFactors = F
)
elpd$province = factor(elpd$province, levels=places)
elpd$npar = sapply(elpd$est, nrow)


# compute BIC
fun$calcBic =  function(prov, ll, np){
    # count number of non-missing data points
    ndata = sum(dat$case[[prov]] != -1)
    log(ndata)*np - 2*ll
}
x[['bic']] = x[['log_lik_sum']] %>%
    rename(lowerLL = lower, upperLL = upper) %>%
    left_join( elpd %>% select(model, province, npar)
        ,by = c('model','province')
    ) %>%
    mutate(
        var = 'bic'
        ,val = mapply( fun$calcBic
            ,prov = province
            ,ll = val
            ,np = npar
        )
        ,lower = mapply( fun$calcBic
            ,prov = province
            ,ll = upperLL # upperbound for LL --> lowerbound for BIC
            ,np = npar
        )
        ,upper = mapply( fun$calcBic
            ,prov = province
            ,ll = lowerLL
            ,np = npar
        )
    ) %>%
    select(-lowerLL, -upperLL)

# Top 3: BIC
getCountryBestModel( x[['bic']], maxRank = 3 )

# make relative to country-level best model
x[['bic']] = x[['bic']] %>%
    left_join(
        x[['bic']] %>%
            filter(model == getCountryBestModel( x[['bic']])) %>%
            select(province, val)
        ,by = "province"
        ,suffix = c("",".best")
    ) %>%
    mutate(
        delta = val - val.best
        , lowerDelta = lower - val.best
        , upperDelta = upper - val.best
    ) %>%
    select(-val.best, -npar)


# Top 3: loglikelihood
getCountryBestModel( x[['log_lik_sum']]
    , smaller_better = -1
    , maxRank = 3
)

# revalue such that the smaller the better
x[['log_lik_sum']] = x[['log_lik_sum']] %>%
    left_join(
        x[['log_lik_sum']] %>%
            filter(model == getCountryBestModel( x[['log_lik_sum']]
                , smaller_better = -1
            )) %>%
            select(province, val)
        ,by = "province"
        ,suffix = c("",".best")
    ) %>%
    mutate(
        # make it the lower the better
        delta = -(val - val.best)
        , lowerDelta = -(upper - val.best)
        , upperDelta = -(lower - val.best)
        , upper.tmp = -lower
        , lower = -upper
        , upper = upper.tmp
    ) %>%
    select(-val.best, -upper.tmp)



# Top 3: RMSE
getCountryBestModel( x[['rmse']], maxRank = 3 )


signedLog = function(x, FUN = log){
    ifelse(x==0, 0, sign(x) * FUN(abs(x)))
}
# compute delta RMSE
x[['rmse']] = x[['rmse']] %>%
    left_join(
        x[['rmse']] %>%
            filter(model == getCountryBestModel( x[['rmse']] )) %>%
            select(province, val)
        ,by = "province"
        ,suffix = c("",".best")
    ) %>%
    mutate(
        delta = val - val.best
        , lowerDelta = lower - val.best
        , upperDelta = upper - val.best
    ) %>%
    select(-val.best) %>%
    # convert RMSE to log10(RMSE) so that plots look more meaningful
    mutate_at(vars(val:upperDelta), signedLog, FUN = log10) %>%
    mutate(var = 'log10(RMSE)')


x = x %>%
    do.call(what=rbind)

metrics = c(
    'rmse' = 'RMSE'
    ,'log10(RMSE)' = 'log10(RMSE)'
    ,'log_lik_sum' = '-LL'
    ,'bic' = 'BIC'
    ,'ELPD' = '-ELPD'
)


x = elpd %>%
    select(-dataver, -est, -npar, -rankVal) %>%
    rbind( ungroup(x)[ ,setdiff(colnames(elpd),c('dataver','est','npar','rankVal'))] ) %>%
    mutate(
        province = factor(province, levels=places)
        ,var = factor(metrics[var], levels=metrics)
    )


#   Figure S2: Model performance rankings
modelBreaks = c(1,4,8,11,14,18) # splits between models with different
                                # number of components
qetBreak = 10
modelSelected = 'm3h.rtaQeT'
modelLetter = LETTERS[seq_along(fun$models)]
names(modelLetter) = fun$models

gCountry = x %>%
    mutate(
        model = factor(model, levels=fun$models)
        ,varEx = factor(var, levels=metrics, labels=paste0('Delta~plain((',metrics,'))'))
    ) %>% 
    ggplot(aes( x = factor(model, levels=fun$models)))+
    geom_vline(xintercept = modelBreaks + .5, size=.1)+
    geom_vline(xintercept = qetBreak + .5, size=.3)+
    geom_jitter(aes( y = delta, color = delta < 0 )
        , alpha = .4
        , shape = 1
        , size = 2
        , width = .2
        , height = 0
    )+
    facet_wrap( ~varEx
        ,scale="free_y",ncol=1, strip.position="right"
        ,labeller = label_parsed
    )+
    theme(
        axis.text.x = element_blank()
        ,axis.ticks.x = element_blank()
        ,axis.title.x = element_blank()
        ,panel.grid.major = element_blank()
        ,panel.grid.minor = element_blank()
        ,plot.margin = unit(c(0.5,0,0,0),"lines")
    )+
    scale_x_discrete("Model", expand=c(.05,.05))+
    scale_y_continuous("Difference from\ncountry-level best model"
        , expand=c(.2,.2)
    )+
    scale_color_manual(values = c("black","red"), guide = "none")

gComponent = lapply(fun$models, fun$getModelComponent) %>%
    do.call(what=rbind) %>%
    mutate(
        var = "Component"
        ,parname = factor(parname, levels=names(fun$modelComponents))
        ,model = factor(model, levels = fun$models)
    ) %>%
    ggplot(aes( x = model, y = parname, alpha = usePar))+
    geom_point(aes(color = model == modelSelected), size=3, shape=15)+
    scale_alpha_manual(values = c(0,1), guide = "none")+
    scale_color_manual(values = c('TRUE' = 'blue', 'FALSE' = 'black'), guide = "none")+
    scale_y_discrete(labels = fun$modelComponents)+
    scale_x_discrete("Candidate Model"
        , expand=c(.05,.05)
        , labels = modelLetter
    )+
    facet_wrap(~var, strip.position="right")+
    geom_vline(xintercept = modelBreaks + .5, size=.1)+
    geom_vline(xintercept = qetBreak + .5, size=.3)+
    theme(
        axis.ticks.x = element_blank()
        ,axis.title.y = element_blank()
        ,panel.grid.major = element_blank()
        ,panel.grid.minor = element_blank()
        ,panel.spacing = unit(0,"lines")
        ,plot.margin = unit(c(0,0,0,0),"lines")
    )


g = ggarrange( gCountry, gComponent
        ,nrow = 2
        ,ncol = 1
        ,align = "v"
        ,heights = c(1.9,.6)
)
ggsave( g
    ,filename = paste0(plotdir,"/performance.pdf")
    ,width = 6
    ,height = 6
)

# variable to store the four models which were best performing across all metrics 
bestModels = grep('m3h\\.rta|m3h\\.rAT', fun$models, value = T)

g = x %>%
    filter(
        model %in% bestModels
        , var == '-ELPD'
    ) %>%
    mutate(model = factor(model, levels = bestModels)) %>%
    ggplot(aes(x = model, y = val, color = model == modelSelected))+
    geom_point()+
    geom_linerange(aes(ymin = lower, ymax = upper))+
    facet_wrap( ~ province, ncol = 8 )+
    scale_color_manual(values = c('TRUE' = 'blue', 'FALSE' = 'black'), guide = "none")+
    scale_x_discrete("Candidate model"
        , labels = modelLetter
    )+
    ylab("Difference in -ELPD from\ncountry-level best model")+
    theme(
        panel.grid = element_blank()
        , strip.text = element_text(size = rel(0.5))
    )
ggsave( g
    ,filename = paste0(plotdir,"/performance_bestModels.pdf")
    ,width = 9
    ,height = 12
)
  

# In how many provinces did the most complex models win?
x %>%
    filter(var == '-ELPD') %>%
    group_by(province) %>%
    filter(val == min(val)) %>%
    ungroup %>%
    with(table(model %in% bestModels))

# In how many provinces did each model converge?
# And, was it the best fitting model when it did converge?
x %>%
    filter(var == '-ELPD') %>%
    filter(model %in% bestModels) %>%
    group_by(province) %>%
    mutate(isBest = val == min(val)) %>%
    with(table(isBest, model)) %>%
    addmargins(1) %>%
    addmargins(2)


    #   Supplement: Model evaluation and selection
    #   ..........................................

# null model
x %>%
    filter(model == 'm1h') %>%
    group_by(var) %>%
    summarize(
        Median = median(val) %>% round(2)
        , lower = quantile(val, 0.025) %>% round(2)
        , upper = quantile(val, 0.975) %>% round(2)
    ) %>%
    as.data.frame


# single component model
fun$checkPerformanceIntervals = function(m1,m2){
    by(x, x$var, function(xx){
        xx %>% filter(model==m1) %>%
            left_join( xx %>% filter(model==m2)
                ,by = c('province')
                ,suffix = c('.m1','.m2')
            ) %>%
            mutate(
                m1Wins = upper.m1 < lower.m2
                , m2Wins = upper.m2 < lower.m1
                , Ties = !(m1Wins|m2Wins)
            ) %>%
            summarize(
                m1 = m1
                , m2 = m2
                , m1Wins = sum(m1Wins, na.rm = T)
                , m2Wins = sum(m2Wins, na.rm = T)
                , Ties = sum(Ties, na.rm = T)
                , nprov_eval = m1Wins + m2Wins + Ties
            )
    }) 
}
# suggests age more important than time
fun$checkPerformanceIntervals('m1h.ra','m1h.rt')

# but not for Phi(a) versus K(t)
fun$checkPerformanceIntervals('m1h.ra','m3h')


# adding second component: age + age
fun$checkPerformanceIntervals('m1h.ra','m3h.ra')

# adding Phi(t) to
fun$checkPerformanceIntervals('m1h.ra','m1h.rta') # Phi(a)
fun$checkPerformanceIntervals('m3h','m3h.rt')     # K(a)

# # combine all
fun$checkPerformanceIntervals('m1h.rta','m3h.rta') # compared to omitting K(a)
fun$checkPerformanceIntervals('m3h.rt','m3h.rta')  # compared to omitting Phi(a)


# replace Phi(a)*Phi(t) with Phi(a,t)
fun$checkPerformanceIntervals('m3h.rta','m3h.rAT')


# selected model compared to all others
fun$checkPerformanceIntervals('m3h.rta','m3h.rtaQeT')
fun$checkPerformanceIntervals('m3h.rAT','m3h.rtaQeT')
fun$checkPerformanceIntervals('m3h.rATQeT','m3h.rtaQeT')
 
 
 
        #   Congruence of estimates
        #   .......................
        setPlotdir('estimates')


    # subset to only best model(s)
    elpd.full = elpd
    elpd = elpd %>% filter(model %in% c('m3h.rta','m3h.rtaQeT'))

    # factor to scale the plotting of lambda and K
    lambda.scaler = 1/5

# function to expand binned variables
fun$expandBinned = function(binMax, binwidth, iMax, leftHang = 0){
    i = c(
        rep(1,leftHang) %>% cumsum
        ,rep((1+leftHang):binMax, each=binwidth)
    )
    if( length(i) >= iMax ){ return( head(i, iMax) ) }
    c(i, rep(binMax, iMax - length(i)))
}

# function to extract parameter values
fun$extractParam = function(x, parname, varname=NULL, ...){
    x = x %>%
        filter(var == parname) %>%
        # --------------------- IMPORTANT ----------------------
        # Use medians instead of means because variables may or may not
        # be normally distributed. Just to be on the safe side.
        mutate(mean = mid) %>%
        # ------------------------------------------------------
        select(var, ivar, mean, lower, upper)

    ctrl = list(...)
    if(length(ctrl) == 0 ){ return(x) }
        
    i = do.call(fun$expandBinned, ctrl)
    x[i, ] %>%
        mutate(
            var = varname
            ,ivar = seq_along(ivar)
        )
}


paramCols = c(
    'Lambda' = 'black'
    , 'K' = "#009e05"
    , 'Qe' = '#7300ff'
    , 'PhiA' = "#c77b00"
    , 'PhiT' = '#c4600e'
)
# function to plot the series of piecewise constant parameters
fun$plotParamSeries = function(x, ylabel
    , legendPosition=c(2,2)
    , Color.qi = '#8f8f8f'
    , Color.qit = paramCols[['Qe']]
){
    elpd$x = x
    el = elpd %>% 
        select( -intersect(colnames(x[[1]]), colnames(elpd)) ) %>%
        mutate(model = factor(model, levels=fun$models)) %>%
        unnest(cols=x)

    el %>%
        ggplot(aes( x = xvar))+
            geom_ribbon(aes(
                    ymin = lower 
                    ,ymax = upper
                    , group = paste(province,model)
                    , fill = model
                )
                , alpha=.005
            )+
            geom_errorbar(data = el %>%
                group_by(model, xvar, ivar) %>%
                summarize(
                     lower = quantile(mean, .025)
                    ,upper = quantile(mean, .975)
                )
                ,aes(
                    ymin=lower
                    ,ymax=upper
                    , col = model 
                ) )+   
            geom_point(data = el %>%
                group_by(model, xvar, ivar) %>%
                summarize( midVal = median(mean) )
                ,aes(
                    y = midVal
                    , col = model 
                ), shape=16, size=2 )+
        ylab(ylabel)+
        scale_color_manual(
            "Model components"
            ,values = c(Color.qi, Color.qit)
            ,label = c(
                expression(paste(phi(t),", ",phi(a),", ",kappa(a),", ",Q(i)))
                ,expression(paste(phi(t),", ",phi(a),", ",kappa(a),", ",Q(i,t)))
            )
        )+
        scale_fill_manual(
            values = c(Color.qi, Color.qit)
            ,guide = "none"
        )+
        theme(
            panel.grid = element_line(color="#f7f5f0")
            ,strip.text = element_text(margin = margin(0,0,0,0, "lines"))
            ,legend.position = legendPosition
            ,plot.margin = margin(2,.5,0,0, "lines")
            ,axis.title.y = element_text(color = Color.qit, size = rel(1.2))
            ,axis.text.y = element_text(color = Color.qit)
        )
}


        
# time varying parameters
gTime = ggarrange(
    # lambda
    fun$plotParamSeries( 
        elpd$expandedLambda <- lapply(elpd$est, function(est){
            fun$extractParam(est, "lambda") %>% 
            mutate(xvar = 1980 - 1 + ivar) %>%
            # scaling
            mutate(
                mean = mean * lambda.scaler
                ,lower = lower * lambda.scaler
                ,upper = upper * lambda.scaler
            )
        })
        #, expression(bar(tau)(t))
        , "Time-varying per-serotype\ninfection intensity"
        , legendPosition=c(.8,.8)
        , Color.qit = paramCols[['Lambda']]
    )+
    scale_x_continuous("Year", limits = c(1980,2017) + c(-.5,.5))

    # PhiT
    ,fun$plotParamSeries(
        elpd$expandedPhiT <- lapply(elpd$est, function(est){
            fun$extractParam(est, "binnedPhiT", varname="PhiT", binMax=16, binwidth=2, iMax=37) %>% mutate(xvar = 1980 + ivar)
        })
        #, expression(phi(t))
        , "Time-varying\nreporting rate"
        , Color.qit = paramCols[['PhiT']]
    )+
    scale_x_continuous("Year", limits = c(1980,2017) + c(-.5,.5))

    ,ncol = 1
    ,nrow = 2
    ,heights = c(1.4,2)
    ,align = "v"
)


# age varying parameters: K, PhiA
gAge = ggarrange(
    # K
    fun$plotParamSeries(
        elpd$expandedK <- lapply(elpd$est, function(est){
            fun$extractParam(est, "binnedK", varname="K", binMax=16, binwidth=4, iMax=65) %>% 
            mutate(xvar = ivar - 1) %>%
            # scaling
            mutate(
                mean = mean / lambda.scaler
                ,lower = lower / lambda.scaler
                ,upper = upper / lambda.scaler
            )
        })
        # , expression(kappa(a))
        , "Age-specific modifier for\ninfection intensity"
        , Color.qit = paramCols[['K']]
    )+
    coord_cartesian(ylim=c(0, 8))+
    scale_x_continuous("Age", limits = c(0,65) + c(-.5,.5))

    # PhiA
    ,fun$plotParamSeries(
        elpd$expandedPhiA <- lapply(elpd$est, function(est){
            fun$extractParam(est, "binnedPhiA", varname="PhiA", binMax=16, binwidth=4, iMax=65) %>% mutate(xvar = ivar - 1)
        })
        # , expression(phi(a))
        , "Age-specific\nreporting rate"
        , Color.qit = paramCols[['PhiA']]
    )+
    scale_x_continuous("Age", limits = c(0,65) + c(-.5,.5))

    ,ncol = 1
    ,nrow = 2
    ,heights = c(1.4,2)
    ,align = "v"
)




# Qe
fun$extractQe = function(varname = 'binnedQe'){
    lapply(elpd$est, function(est){
        fun$extractParam(est, varname)
    })
}
elpd$Qe = fun$extractQe()

qeDat = elpd %>%
    select(model, Qe) %>%
    mutate(model = factor(model, levels=fun$models)) %>%
    unnest(cols= Qe) %>%
    group_by(model,var,ivar) %>%
    summarize(
        midVal = median(mean)
        ,upper = quantile(mean, .975)
        ,lower = quantile(mean, .025)
    ) %>%
    mutate(
        epoch = ceiling(ivar/4)
        ,ivar = (ivar %% 4)
        ,ivar = ivar + (ivar==0)*4
    )

fun$plotQe = function(
    qeDatI
    ,xlabel = "i-th Infection"
    ,ylabel = "Clinical detectability (Q)"
    ,ylimit = c(0,1)
    ,fill.labs = c(
                "Constant detectability"
                ,paste0(
                    seq(1981,2012,by=9)
                    ,"-"
                    ,c(seq(1981+8,2012,by=9),2017)
                )
            )
    ,fill.colors = c(
                '#444444'
                , '#d5c3eb'
                , '#a988d1'
                , '#874bd1'
                , paramCols[['Qe']]
            )
){
    qeDatI %>%
        ggplot(aes( x = ivar, y = midVal, fill = paste(model,var,epoch) , col = paste(model,var,epoch)  ))+
        geom_bar(stat="identity", position = "dodge", size=0.2, width=.8
            , color = '#444444'
        )+
        geom_errorbar(aes(ymin = lower, ymax = upper), position = position_dodge(width=.8), col="black", width=.4, size=0.4)+
        scale_fill_manual(
            values = fill.colors
            ,labels = fill.labs
        )+
        scale_x_continuous( xlabel
            ,breaks = 1:4
        )+
        scale_y_continuous( ylabel
            ,expand = c(0,0)
        )+
        coord_cartesian(ylim = ylimit)
}
gQe = ggarrange(
    fun$plotQe(qeDat, ylimit = c(.91,1), xlabel="", ylabel = "") + 
        theme(
            legend.title = element_blank()
            ,panel.grid.major.x = element_blank()
            ,panel.grid.minor.x = element_blank()
            ,panel.border = element_blank()
            ,axis.text.x = element_blank()
            ,axis.ticks.x = element_blank()
            ,plot.margin = margin(1.5,.5,-0.2,0, "lines")
        ) +
        annotate("text"
            , x = 4.4, y = .975
            , label = "Clinical detectability"
            , hjust=1, vjust=1
            , color = "#7300ff"
            , size = 8
        )

    ,fun$plotQe(qeDat, ylimit = c(0,.31), ylabel = "") + 
        theme(
            legend.title = element_blank()
            ,panel.grid.major.x = element_blank()
            ,panel.grid.minor.x = element_blank()
            ,panel.border = element_blank()
            ,plot.margin = margin(0,.5,0,0, "lines")
        )
    ,nrow = 2
    ,ncol = 1
    ,heights = c(1,1.2)
    ,align = "v"
    ,common.legend = T
    ,legend = "right"
)

g = ggarrange(
        ggarrange( gTime, gAge
            ,ncol = 2
            ,nrow = 1
        )
        ,gQe
    ,nrow = 2
    ,ncol = 1
    ,heights = c(3,1)
)

ggsave( g
    ,filename = paste0(plotdir,"/piecewise_colored.pdf")
    ,width = 10
    ,height = 9.5
)


#   Congruence of estimates between
#   Q(i) and Q(i,t) variant of the best model
#   .........................................


elpd.model = split(elpd %>% 
    select(province, model, starts_with('expanded'))
    , elpd$model
)

g = mapply(function(param, transfunc){
        ep = paste0('expanded',param)
        lapply(elpd.model, function(el){
            el %>% 
                select_at(vars(c('model', ep))) %>% 
                unnest(cols = ep) %>% 
                select(-ivar)
        }) %>%
        Reduce(f = function(x,y){
            left_join(x, y, by = c('var','xvar'), suffix = c('.1','.2'))
        }) %>%
        ggplot(aes(x = mean.1, y = mean.2))+
        geom_point(shape = 1, alpha = 0.02)+
        geom_abline(slope = 1, intercept = 0, color = '#cccccc', linetype = 2)+
        scale_x_continuous('Model with constant Q(i)', trans = transfunc)+
        scale_y_continuous('Model with time-varying Q(i,t)', trans = transfunc)+
        theme(panel.grid = element_blank())+
        coord_fixed(ratio = 1)+
        ggtitle(fun$modelComponents[param])
    }
    , param = c('Lambda','K','PhiT','PhiA')
    , transfunc = c('log10','log10','identity','identity')
    , SIMPLIFY = F
)
g = ggarrange(plotlist = g
    , nrow = 2
    , ncol = 2
)
ggsave( g
    ,filename = paste0(plotdir,"/pairs.png")
    ,width = 10
    ,height = 9.5
    ,dpi = 300
)



    #   filter to only best model
    #   .........................

    elpd = elpd %>%
        select(-var, -rankVal) %>%
        filter(model == modelSelected)



    #   Patterns in the estimates
    #   (Figure 2)
    #   .........................

# tau(t)
lapply(1:37, function(time){
    elpd$expandedLambda %>%
        do.call(what = rbind) %>%
        filter(ivar == time) %>%
        with(quantile(mean, c(.025, .5, .9750)))
}) %>%
do.call(what = rbind)

# trough in Phi(a)
elpd$expandedPhiA %>%
    do.call(what = rbind) %>%
    filter(ivar == 13) %>%
    with(quantile(mean, c(.025, .5, .975))) %>%
    round(2)
# Phi(a=50) compared to Phi(a=20)
elpd %>%
    select(province, expandedPhiA) %>%
    unnest(cols = expandedPhiA) %>%
    filter(xvar %in% c(20,50)) %>%
    select(province, xvar, mean) %>%
    spread(xvar, mean) %>%
    mutate(mult = `50`/`20`) %>%
    with(c(mean(mult), quantile(mult, c(.025, .5, .975)))) %>%
    round(2)

# K(a=12) compared to K(a=11)
elpd %>%
    select(province, expandedK) %>%
    unnest(cols = expandedK) %>%
    filter(xvar %in% c(11,12)) %>%
    select(province, xvar, mean) %>%
    spread(xvar, mean) %>%
    mutate(mult = `12`/`11`) %>%
    with(c(median(mult), quantile(mult, c(.025, .5, .975)))) %>%
    round(2)
# K in older adults
elpd %>%
    select(province, expandedK) %>%
    unnest(cols = expandedK) %>%
    filter(xvar %in% 30:70) %>%
    group_by(xvar) %>%
    summarize(median(mean)) %>%
    as.data.frame


dCase = dat$case %>%
    lapply(function(x){
        colSums(x * (x!=-1))
    }) %>%
    do.call(what = rbind) %>%
    as.data.frame %>%
    mutate(province = names(dat$case)) %>%
    gather(Year, count, -province) %>%
    mutate(Year = gsub("[^0-9]","",Year) %>% as.integer) %>%
    mutate(province = factor(province, levels = names(dat$case)))
gCase = dCase %>%
    ggplot(aes(x = Year, y = count * 10 / 1000))+
    #geom_line(aes(group=province), color = '#bbbbbb', alpha = 0.4, size = 0.4)+
    geom_line(data = dCase %>% 
            group_by(Year) %>%
            summarize(count = sum(count))
        , aes(y = count / 1000)
        , size = 0.8
    )+
    scale_x_continuous(limits = c(1980,2017) + c(-.5,.5))+
    scale_y_continuous("1000 DHF cases\nreported")+
    theme(
        panel.grid = element_blank()
        , axis.text.x = element_blank()
        , axis.title.x = element_blank()
    )


gQe = qeDat %>%
    filter(model == modelSelected) %>%
    fun$plotQe(
        fill.labs = paste0(
                    seq(1981,2012,by=9)
                    ,"-"
                    ,c(seq(1981+8,2012,by=9),2017)
                )
        , fill.colors = c(
                '#d5c3eb'
                , '#a988d1'
                , '#874bd1'
                , paramCols[['Qe']]
            )
        # , ylabel = expression('Q(i,t)')
        , ylabel = "Clinical\ndetectability"
    )+
    guides(fill=guide_legend(keyheight=unit(0.4,"lines")))+
    theme(
        legend.title = element_blank()
        ,legend.position = c(1,1)
        ,legend.justification = c(1,1)
        ,panel.grid.major.x = element_blank()
        ,panel.grid.minor.x = element_blank()
        ,panel.border = element_blank()
        # ,plot.margin = margin(0,.5,0,0, "lines")
        ,axis.title.y = element_text(color = paramCols[['Qe']], size = rel(1.2))
        ,axis.text.y = element_text(color = paramCols[['Qe']])
    )


# time varying parameters
gTime = ggarrange(
    gCase
    # tau(t)
    , fun$plotParamSeries( 
        elpd$expandedLambda
        # , expression(bar(tau)(t))
        , "Time-varying per-serotype\ninfection intensity"
        , Color.qit = paramCols[['Lambda']]
        , Color.qi = character(0)
    )+
    scale_x_continuous("Year", limits = c(1980,2017) + c(-.5,.5))+
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

    # PhiT
    ,fun$plotParamSeries(
        elpd$expandedPhiT
        # , expression(phi(t))
        , "Time-varying\nreporting rate"
        , Color.qit = paramCols[['PhiT']]
        , Color.qi = character(0)
    )+
    scale_x_continuous("Year", limits = c(1980,2017) + c(-.5,.5))+
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

    ,ncol = 1
    ,nrow = 3
    ,heights = c(1.2, 2, 2)
    ,align = "v"
    ,labels = 'auto'
)


# age varying parameters: K, PhiA
# (plot as age bins)
ageBinLabel = data.frame(ageBin = fun$expandBinned(binMax = 16, binwidth = 4, iMax = 65)) %>%
    mutate(age = seq_along(ageBin) - 1) %>%
    group_by(ageBin) %>%
    nest %>%
    mutate(ageBinLabel = sapply(data, function(x){
        x = range(x$age)
        if(x[1]>=60){ return('60+')}
        paste0(x[1],'-',x[2])
    })) %>%
    with(ageBinLabel)
gAge = ggarrange(
    gQe
    # K
    ,fun$plotParamSeries(
        lapply(elpd$est, function(est){
            fun$extractParam(est, "binnedK") %>%
            mutate(xvar = factor(ageBinLabel[ivar], levels = ageBinLabel)) %>%
            # scaling
            mutate(
                mean = mean / lambda.scaler
                ,lower = lower / lambda.scaler
                ,upper = upper / lambda.scaler
            )
        })
        # , expression(kappa(a))
        , "Age-specific modifier for\ninfection intensity"
        , Color.qit = paramCols[['K']]
        , Color.qi = character(0)
    )+
    coord_cartesian(ylim=c(0, 8))+
    xlab('Age')+
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

    # PhiA
    ,fun$plotParamSeries(
        lapply(elpd$est, function(est){
            fun$extractParam(est, "binnedPhiA") %>% 
            mutate(xvar = factor(ageBinLabel[ivar], levels = ageBinLabel))
        })
        #, expression(phi(a))
        , "Age-specific\nreporting rate"
        , Color.qit = paramCols[['PhiA']]
        , Color.qi = character(0)
    )+
    xlab('Age')+
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

    ,ncol = 1
    ,nrow = 3
    ,heights = c(1.2,2,2)
    ,align = "v"
    ,labels = letters[4:6]
)

g = ggarrange(
      gTime
    , gAge
    , ncol = 2
    , nrow = 1
    , widths = c(2, 1.5)
)

ggsave( g
    ,filename = paste0(plotdir,"/piecewise_colored_best.pdf")
    ,width = 7
    ,height = 8
)




    #   Correlation in the estimates across provinces
    #   (Figure S6)
    #   .............................................

parnames = c(
    "expandedLambda" = expression( bar(tau)(t) )
    , "expandedK" = expression( kappa(a) )
    , "expandedPhiT" = expression( phi(t) )
    , "expandedPhiA" = expression( phi(a) )
)

plotCor = function(elpd){

    getCorr = function(
            parname
            , distRescalefun=function(x){x}
    ){
        # compute pair-wise correlation coefficients
        d = elpd[ ,c("province",parname)] %>%
            unnest(!!parname) %>%
            select(province,ivar,mean) %>%
            spread(province, mean) %>%
            ungroup %>%
            select(-ivar)
        corr = cor(d)
        data.frame(cor = corr[upper.tri(corr, diag=F)])
    }

    parCor = lapply(names(parnames), getCorr)
    gCor = mapply(function(corr, parname){
            g = ggplot(corr, aes(x = cor))+
                geom_histogram(aes(y = ..density..), binwidth = 0.05)+
                geom_vline(aes(xintercept = median(cor)), color = "blue", lty = "21") +
                scale_x_continuous("Correlation coef."
                    , limits = c(0,1)
                )+
                annotate("text", x = 0, y = 3
                    , label = paste('Median =', round(median(corr$cor),2))
                    , color = 'blue'
                    , vjust = 1
                    , hjust = 0
                )+
                ggtitle(parnames[[parname]])
            return(g)
        }
        , corr = parCor
        , parname = names(parnames)
        , SIMPLIFY = FALSE
    )
    ggarrange(
        plotlist = gCor
        ,ncol = 1
        ,nrow = length(gCor)
        ,align = 'hv'
    )
}

library(geosphere)
plotCorDist = function(
        parname
        , posVar=c('cenX','cenY')
        , xlabel="Euclidean distance (km)"
        , distfun=distm
        , distRescalefun=function(x){x}
){
    # compute pair-wise distances
    d = sp$provinceCentroids %>%
        mutate(province72 = factor(province72, levels=places)) %>%
        group_by(province72) %>%
        summarize(
            cenX = sum(cenX*area)/sum(area)
            ,cenY = sum(cenY*area)/sum(area)
        )
    rownames(d) = d$province72
    
    distance = distfun(d[ , posVar]) / 1000 # calculate distance and convert to km
    distance = as.matrix(distance)

    # compute pair-wise correlation coefficients
    d = elpd[ ,c("province",parname)] %>%
        unnest(parname) %>%
        select(province,ivar,mean) %>%
        spread(province, mean) %>%
        ungroup %>%
        select(-ivar)

    d = data.frame(
        cor = cor(d)[upper.tri(distance, diag=F)]
        ,distance = distance[upper.tri(distance, diag=F)]
    ) %>%
    mutate( distance = distRescalefun(distance) )

    g = ggplot(d, aes(
            x = distance
            ,y = cor
        ))+
        geom_point(alpha=.3, shape=1)+
        geom_smooth(method="gam")+
        labs(
            x = xlabel
            ,y = "Corr. coef."
        )+
        ggtitle(parnames[[parname]])
    return(g)
}

g = ggarrange(
    plotCor(elpd)
    , ggarrange(
        plotCorDist("expandedLambda")
        ,plotCorDist("expandedK")
        ,plotCorDist("expandedPhiT")
        ,plotCorDist("expandedPhiA")
        ,ncol = 1
        ,nrow = 4
    )
    , ncol = 2
    , nrow = 1
    , labels = 'auto'
)
ggsave( g
    ,filename = paste0(plotdir,'/est-correlation.pdf')
    ,width = 6
    ,height = 8
)



    #   State reconstruction
    #   ....................
    library(abind)
    setPlotdir('reconstruction')

infectionCols = c(
    "#edd312"
    ,"#ed5412"
    ,"#4a7500"
    ,"#6a12ed"
)

# state recon for one particular cohort
fun$calc.profile <- function(lambdas){
    prop = exp(-cumsum(lambdas))    # proportion remain susceptible to serotype i
    x = prop^4        # proportion susceptible to all serotypes
    z1 = 4 * prop^3 * (1 - prop)^1
    z2 = 6 * prop^2 * (1 - prop)^2
    z3 = 4 * prop^1 * (1 - prop)^3
    z4 = (1 - prop)^4
    rbind(x,z1,z2,z3,z4)
    # out: i x time
}
fun$calc.infection <- function(lambdas){
    profile = cbind(c(1,0,0,0,0),fun$calc.profile(lambdas))
    infection = matrix(apply(
        apply(profile,2,function(x){
            1-cumsum(x)
        })
    ,1, diff), nrow=length(lambdas), ncol=5)
    matrix(
        infection[,1:4]
        , ncol=4
        , dimnames=list(NULL,paste0('i',1:4))
    )
    # out: time x i
}

# for all cohorts born since the 1st given lambda
fun$getInfectedProportionCohort = function(lambdas, K=rep(1,length(config$ages))){
    # for each cohort h
    lapply( seq_along(lambdas), function(h, iMaxAge=length(config$ages)){
        iLambdas = h : length(lambdas)
        iLambdas = iLambdas[ -( (iMaxAge+1):(length(iLambdas)+1) ) ]
        out = fun$calc.infection(c(
            rep(0,h-1)  # before the cohort was born
            ,lambdas[ iLambdas ] * K[ seq_along(iLambdas) ] # during when they lived
        ))
        out %>%
        rbind( matrix(0, length(lambdas)-nrow(out), 4) ) # after they die
    }) %>%
    abind(along=3)
    # out: time x i x cohort
}

# convert (time x cohort) to (time x age)
fun$th2ta = function(x, iMaxAge = length(config$ages)){
    times = 1:nrow(x)
    iMaxCohort = ncol(x)
    # for each age
    lapply( 1:iMaxAge, function(ia){
        # for each time
        sapply( times, function(it){
            ih = it - ia + 1
            ifelse( (ih > iMaxCohort)|(ih < 1), 0, x[ it, ih])

        })
        # out: time x i
    }) %>%
    do.call(what=cbind)
    # out: time x age
}

# function to case counts over time, by age
# : for each i
fun$getCaseAge = function(lambdas, K, phiTA=1, Qit, pop){
    xx = fun$getInfectedProportionCohort(lambdas, K) %>% # time x i x cohort
        aperm( c(1,3,2) )   # time x cohort x i

    # population counts
    pop = cbind(
            matrix(0, nrow(pop), nrow(pop) - ncol(pop))
            ,pop    # cohort x time
        ) %>% t # time x cohort

    # a list of counts of each i-th infection: cohort x time
    xx = lapply(1:4, function(i){ 
        xx[ , , i] * pop
    })

    # a list of counts of each i-th infection: time x age
    xx = lapply(xx, fun$th2ta)

    # a list of counts of each i-th "reported" infection: time x age
    # : using default phiTA=1 returns
    # : a list of counts of each i-th "clinical" infection: time x age
    mapply(function(xi, qi){
            # trim years to only observed years
            xi = xi[ nrow(xi) - (length(config$times):1) + 1, ]
            (xi * qi) * phiTA
        }
        ,xi = xx    # time x age
        ,qi = Qit
        ,SIMPLIFY = FALSE
    )
}

# function to get case mean age from given parameter values
# : for each i as well as combined
fun$getCaseMeanAge = function(lambdas, K, phiTA=1, Qit, pop){
    xx = fun$getCaseAge(lambdas = lambdas, K = K, phiTA = phiTA, Qit = Qit, pop = pop)

    # compute mean age for reported cases of each i
    # : time x i
    meanAge = lapply(xx, function(xi){
        (xi %*% config$ages) / rowSums(xi)
    })
    
    # compute mean age for reported cases (all i combined)
    # : time
    meanAge[['combined']] =
        (
            lapply(xx, function(xi){
                (xi %*% config$ages)
            }) %>%
            do.call(what=cbind) %>%
            rowSums
        ) /
        (
            lapply(xx, rowSums) %>%
            do.call(what=cbind) %>%
            rowSums
        )
    
    return(meanAge)
}


# appropriately repeat expanded estimates to match full length
# of those parameters
fun$extendExpandedToFullLength = function(exLambda, exK, exPhiA, exQ){
    if((length(exQ) %% 4) != 0){
        stop('Q not in multiples of 4.')
    } else {
        exQ = split(exQ, seq_along(exQ) %% 4)[c('1','2','3','0')] %>%
            lapply(function(x){ x[fun$expandBinned(length(exQ)/4, 9, iMax = length(config$times))] })
    }

    list(
        exLambda = c(
            rep(exLambda[1], length(config$cohorts) - length(config$times))
            ,exLambda[-1]
        )
        , exK = c(
            exK
            ,rep(exK[length(exK)], length(config$ages) - length(exK))
        )
        , exPhiA = c(
            exPhiA
            ,rep(exPhiA[length(exPhiA)], length(config$ages) - length(exPhiA))
        )
        , exQ = exQ
    )
}

# get mean age from expanded estimates
# (i.e. wrapper function around fun$getCaseMeanAge)
# : assuming "m3h.rta" model
fun$getCaseMeanAgeFromExpanded = function(exLambda, exK, exPhiA, exPhiT, exQ, pop){
    fullLength = fun$extendExpandedToFullLength(exLambda, exK, exPhiA, exQ)
    with(fullLength, fun$getCaseMeanAge(
        exLambda
        ,exK
        ,as.matrix(exPhiT) %*% t(exPhiA)
        ,exQ
        ,pop
    ))
}

# get age from expanded estimates
# (i.e. wrapper function around fun$getCaseAge)
# : assuming "m3h.rta" model
fun$getCaseAgeFromExpanded = function(exLambda, exK, exPhiA, exPhiT, exQ, pop){
    fullLength = fun$extendExpandedToFullLength(exLambda, exK, exPhiA, exQ)
    with(fullLength, fun$getCaseAge(
        exLambda
        ,exK
        ,as.matrix(exPhiT) %*% t(exPhiA)
        ,exQ
        ,pop
    ))
}


# get the dataframe version of the mean ages, with year
fun$getMeanAgeDf = function(params, ver, prov){    
    do.call(fun$getCaseMeanAgeFromExpanded, params) %>% 
        do.call(what=cbind) %>%
        as.data.frame %>%
        mutate(year = seq_along(combined) + 1980) %>%
        mutate(ver = ver, province = prov)
}

# compute reconstruction
fun$getReconstruction = function(lambdas, K, prov){

    pop = dat$pop[[prov]]
    
    # get S going into time t: S(t-1/26)
    S = lapply(seq_along(config$cohorts), function(h){
        ht = h + config$ages
        ht = ht[ht <= length(lambdas)]

        lam = lambdas[ht]
        lam = lam * K[1:length(lam)]
        lam = c(rep(0, max(0,length(config$times) - length(lam))), lam)
        S = fun$calc.profile(cumsum(lam) - lam/26) %>% t
        # filter to only years in the data
        it = nrow(S) - (length(config$times):1) + 1
        S = S[ it, ]
        # convert to counts
        S = S * pop[h, ]
        
        list(
            S = S
            # store lambda * count
            # to compute average lambda later
            ,SxLambda = ((lam[it]/26) %*% t(4:0)) * S
        )
    })
    
    SxLambda = lapply(S, function(x){ x[['SxLambda']] }) %>%
        Reduce(f="+")
    
    S = lapply(S, function(x){ x[['S']] }) %>%
        Reduce(f="+")

    # get lambda_bar(t)    
    lambdas = rowSums(SxLambda) / rowSums(S[ ,-5])
    
    # get I(t) = lambda(t) * S(t-1/26)
    I = S * ((lambdas/26) %*% t(4:0))

    # solve for beta(t) = lambda(t)/I(t)
    Beta = (lambdas/26)/ (rowSums(I)/rowSums(S))

    data.frame(
        province = prov
        ,"S" = S
        ,"I" = I
        ,Beta = Beta
        ,lambda = lambdas
        ,Year = seq_along(Beta) + 1980
    )
}



# do the actual reconstruction computation
recon = mapply( fun$getReconstruction
    ,lambdas = lapply(elpd$expandedLambda, function(x){
        lambdas = x$mean
        c(
            rep(lambdas[1], length(config$cohorts) - length(config$times))
            ,lambdas[-1]
        )
    })
    ,K = lapply(elpd$expandedK, function(x){
        K = x$mean
        c(
            K
            ,rep(K[length(K)], length(config$ages) - length(K))
        )
    })
    ,prov = as.character(elpd$province)
    ,SIMPLIFY = FALSE
)


# take a look at how well the reconstructions fit the observed data
Error = with(elpd, {
    mapply(
        fun$getCaseAgeFromExpanded
        , exLambda = lapply(expandedLambda, function(x) x$mean )
        , exK = lapply(expandedK, function(x) x$mean )
        , exPhiA = lapply(expandedPhiA, function(x) x$mean )
        , exPhiT = lapply(expandedPhiT, function(x) x$mean )
        , exQ = lapply(Qe, function(x) x$mean )
        , pop = dat$pop[province]
        , SIMPLIFY = FALSE
    ) %>% 
    mapply(FUN = function(x, prov){
            case = dat$case[[prov]]
            # for each i-th infection, calculate counts by age group
            predicted = lapply(x, function(xi){
                    as.matrix(config$ageGroups) %*% t(xi)
                }) %>%
                # sum counts from all i's
                abind(along=3) %>%
                apply(c(1,2), sum)

            tibble(
                ageGroup = rep(1:nrow(config$ageGroups), times = ncol(case))[case != config$naVal]
                , year = rep(config$times, each = nrow(config$ageGroups))[case != config$naVal]
                , predicted = predicted[case != config$naVal]
                , observed = case[case != config$naVal]
            ) %>%
            group_by(year) %>%
            mutate(error = sum(abs(predicted - observed))) %>%
            ungroup %>%
            mutate(province = prov)
            
        }
        , prov = province
        , SIMPLIFY = FALSE
    )
})

plotGof = function(worst = T, dat = Error){
    g = lapply(dat, function(x, errorMin = !worst){
        x %>%
            arrange(error * ifelse(errorMin, 1, -1)) %>%
            filter(error == error[1])
    }) %>%
    do.call(what = rbind) %>%
    mutate(`Age Group` = factor(ageGroup
        , levels = seq_along(config$ageGroupLabels)
        , labels = config$ageGroupLabels
    )) %>%
    ggplot(aes(x = `Age Group`))+
    geom_col(aes(y = observed), fill = '#cccccc')+
    geom_point(aes(y = predicted), size = 0.5, shape = 1)+
    facet_wrap( ~ province, ncol = 8)+
    ylab('Number of DHF cases')+
    theme_classic()+
    theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = rel(0.7))
        , strip.text = element_text(size = rel(0.7))
    )
    return(g)
}

# worst fitting year of each province
g = plotGof(worst = T)
ggsave( g
    ,filename = paste0(plotdir,"/gof_bestModel_worst.pdf")
    ,width = 12
    ,height = 16
)
# best fitting year of each province
g = plotGof(worst = F)
ggsave( g
    ,filename = paste0(plotdir,"/gof_bestModel_best.pdf")
    ,width = 12
    ,height = 16
)

g = ggarrange(
    # error by age
    Error %>%
        do.call(what = rbind) %>%
        mutate(
            error = predicted - observed
            , ageGroup = factor(ageGroup, labels = config$ageGroupLabels)
        ) %>%
        ggplot(aes(x = ageGroup, y = error/observed, group = ageGroup))+
        geom_violin(fill = 'black')+
        ylab('Error : Observed')+
        xlab('Age Group')+
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

    # Error by year
    , Error %>%
        do.call(what = rbind) %>%
        mutate(
            error = predicted - observed
        ) %>%
        ggplot(aes(x = year, y = error/observed, group = year))+
        geom_violin(fill = 'black')+
        ylab('Error : Observed')+
        xlab('Year')+
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

    , nrow = 2
    , ncol = 1
    , labels = 'auto'
)
ggsave( g
    ,filename = paste0(plotdir,"/gof_bias.pdf")
    ,width = 6
    ,height = 8
)






# describe the reconstructed susceptibility states
    # : countrywide
    recon %>%
        do.call(what=rbind) %>%
        filter(Year %in% c(1981,2010)) %>%
        select(province, Year, S.x:S.z4) %>%
        gather(state, count, -Year, -province) %>%
        group_by(Year) %>%
        mutate( pop = sum(count) ) %>%
        group_by(Year, state) %>%
        summarize( prop = sum(count) / unique(pop) ) %>%
        spread(Year, prop) %>%
        mutate(
            fold.change = (`2010`-`1981`) / `1981`
        )
        
    # : average across provinces (just to ensure that they are similar)
    recon %>%
        do.call(what=rbind) %>%
        filter(Year %in% c(1981,2010)) %>%
        select(province, Year, S.x:S.z4) %>%
        gather(state, count, -Year, -province) %>%
        group_by(province, Year) %>%
        mutate(
            prop = count/sum(count)
        ) %>%
        group_by(Year, state) %>%
        summarize(meanProp = mean(prop)) %>%
        spread(Year, meanProp) %>%
        mutate(
            fold.change = (`2010`-`1981`) / `1981`
        )

# describe the reconstructed infectious fractions

    # : countrywide
    recon %>%
        do.call(what=rbind) %>%
        filter(Year %in% c(1981,2010)) %>%
        select(province, Year, I.x:I.z4) %>%
        gather(state, count, -Year, -province) %>%
        group_by(Year, state) %>%
        summarize(count = sum(count)) %>%
        # get the population size
        mutate(pop = sapply(Year, function(year){
            sapply(dat$pop, function(pop){
                pop[ , paste0('X',year)]
            }) %>% sum
        })) %>%
        mutate(prop = count/pop) %>%
        select(-count,-pop) %>%
        ungroup() %>%
        spread(Year, prop) %>%
        mutate(
            fold.change = (`2010`-`1981`) / `1981`
        )



# assess the trends in I(t)
recon.I = lapply(recon, function(x){
    x %>% 
        select(starts_with("I."), Year) %>%
        gather(state, count, -Year) %>%
        group_by(Year) %>%
        summarize( count = sum(count) ) %>%
        mutate( pop = c(dat$pop[[ x$province[1] ]] %>%
                colSums
            )[ paste0('X',Year) ]
            , I = count/pop
        )
})
lapply(recon.I, function(x){

    fit = glm( I ~ Year
        , data = x %>% filter(between(Year, 1981, 2004))
    ) %>% summary
    out = coef(fit)[2,1] + c(-2,0,2) * coef(fit)[2,2]
    c(
        include.zero = between(0, out[1], out[2])
        , below.zero = out[2] < 0
        , above.zero = out[1] > 0
    )  
}) %>%
do.call(what = rbind) %>%
summary


# plot: Beta
gBeta = recon %>%
    do.call(what=rbind) %>%
    ggplot(aes( x = Year, y = Beta))+
    geom_line(aes(group = province), alpha=.25)+
    geom_smooth(method = 'gam')+
    scale_x_continuous(expand=c(0,0))+
    scale_y_continuous( expression(beta(t))
    ) +
    theme(panel.grid = element_blank())

# plot: marginal lambda
gLambda = recon %>%
    do.call(what=rbind) %>%
    ggplot(aes( x = Year, y = lambda))+
    geom_line(aes(group = province), alpha=.15)+
    geom_smooth(method = 'gam')+
    scale_x_continuous(expand=c(0,0))+
    theme(panel.grid = element_blank())+
    ylab("Biweekly FoI")

# plot: log(Beta) x log(lambdas)
gBetaLambda = recon %>%
    do.call(what=rbind) %>%
    ggplot(aes( x = lambda, y = Beta))+
    geom_abline(slope = 1, intercept = -(20:0), col="#dddddd", lty=3)+
    geom_point(aes(col = Year), alpha=.2)+
    geom_smooth(method="gam")+
    scale_x_continuous('Biweekly FoI' #expression(bar(lambda)(t))
        , trans = "log10"
        , expand=c(0,0)
        , breaks = function(x){
            x = log10(x)
            10 ^ seq(floor(x[1]), ceiling(x[2]), by = 1)
        }
        , labels = scales::trans_format("log10", scales::math_format(10^.x))

    )+
    scale_y_continuous(expression(beta(t)))+
    annotation_logticks(sides = "tb"
        , mid = unit(0.1, "cm")
        , short = unit(0.1, "cm")
    )+
    scale_color_gradient2(
        low = "blue"
        ,mid = "grey"
        ,high = "red"
        ,midpoint = 2000
    )+
    theme(panel.grid = element_blank())


# plot: S trends
gS = recon %>%
    do.call(what=rbind) %>%
    select(province, Year, S.x:S.z4) %>%
    gather(state, count, -Year, -province) %>%
    group_by(province, Year) %>%
    mutate(
        prop = count/sum(count)
        ,state = sapply(state, function(x){
            paste0(
                "Susceptible to\n"
                ,4-sum(as.integer(gsub("[^0-9]","",x)),na.rm=T)
                ," serotype(s)"
            )
        })
        ,state = factor(state, levels=sort(unique(state), decreasing=T))
    ) %>%

    ggplot(aes( x = Year, y = prop, col = state, group = province))+
    geom_line(alpha=.2)+
    facet_wrap(~ state, scale="free_y", nrow=1)+
    scale_x_continuous(expand=c(0,0))+
    scale_y_continuous("Proportion\nof population", expand=c(0,0))+
    scale_color_manual(values=c(infectionCols,"black"))+
    guides(color=F)

# plot: I proportions x time
gIprop = recon %>%
    do.call(what=rbind) %>%
    select(province, Year, I.x:I.z4) %>%
    gather(state, count, -Year, -province) %>%
    group_by(province, Year) %>%
    mutate(
        prop = count/sum(count)
        ,state = sapply(state, function(x){
            paste("Infection",sum(as.integer(gsub("[^0-9]","",x)),na.rm=T)+1)
        })
    ) %>%
    filter( prop != 0 ) %>%

    ggplot(aes( x = Year, y = prop, col = state, group = province))+
    geom_line(alpha=.2)+
    facet_wrap(~ state, scale="free_y", nrow=1)+
    scale_x_continuous(expand=c(0,0))+
    scale_y_continuous("Proportion\nof infected", expand=c(0,0))+
    scale_color_manual(values=infectionCols)+
    guides(color=F)

# plot: I x time
gI = recon %>%
    do.call(what=rbind) %>%
    mutate( popsize = S.x + S.z1 + S.z2 + S.z3 + S.z4 ) %>% 
    select(province, Year, I.x:I.z3, popsize) %>%
    gather(state, count, -Year, -province, -popsize) %>%
    group_by(Year, state, popsize) %>%
    summarize( count = sum(count) ) %>%
    ungroup %>%
    mutate(
        state = sapply(state, function(x){
            paste("Infection",sum(as.integer(gsub("[^0-9]","",x)),na.rm=T)+1)
        })
        ,prop = count/popsize
    ) %>%
    mutate(state = factor(state, unique(state) %>% sort(decreasing=T))) %>%
    ggplot(aes( x = Year, y = prop, fill = state))+
    geom_col(width=1)+
    scale_x_continuous(expand=c(0,0))+
    scale_y_continuous(expand=c(0,0))+
    scale_fill_manual(values=infectionCols %>% rev)+
    guides(fill=F)

# Figure 4  
g = ggarrange(
    gS
    ,ggarrange(
        gIprop
        ,gI + ylab("Infected : Population")
        ,ncol = 2
        ,nrow = 1
        ,widths = c(4,1.2)
    )
    ,ggarrange( gLambda, gBeta, gBetaLambda
        ,ncol = 3
        ,nrow = 1
        ,widths = c(1,1,1.5)
        ,labels = c('c','d','e')
    )

    ,nrow = 3
    ,ncol = 1
    ,labels = c('a','b','')
)
ggsave(g
    ,filename = paste0(plotdir,"/states.pdf")
    ,width = 8
    ,height = 6
)




g = ggarrange(
    ggarrange(
        gS
        ,ggarrange(
            gIprop
            ,gI + ylab("Infected : Population")
            ,ncol = 2
            ,nrow = 1
            ,widths = c(4,1.2)
        )
        ,nrow = 2
        ,ncol = 1
        ,labels = c('a','b')
    )
    , ggarrange(gLambda, gBeta
        , nrow = 2
        , ncol = 1
        , align = 'v'
        , labels = c('c','d')
    )
    , ncol = 2
    , nrow = 1
    , widths = c(5,1.1)
)
ggsave(g
    ,filename = paste0(plotdir,"/states.pdf")
    ,width = 11
    ,height = 4
)



    #   Driver silencing
    #   ................
    setPlotdir('silencing')
    
# get mean age for all the configurations
fun$getMeanAgeForPlot = function(iProv){
    prov = elpd$province[iProv] 
    meanAge = list()


    fittedPar = list(
        exLambda    = elpd$expandedLambda[[iProv]]$mean
        ,exK        = elpd$expandedK[[iProv]]$mean
        ,exPhiA     = elpd$expandedPhiA[[iProv]]$mean
        ,exPhiT     = elpd$expandedPhiT[[iProv]]$mean
        ,exQ        = elpd$Qe[[iProv]]$mean
        ,pop        = dat$pop[[ elpd$province[iProv] ]]
    )
    mutatedPar = fittedPar

    silencedPar = list(
        exLambda    = c(mutatedPar$exLambda[1], rep(mean(mutatedPar$exLambda[-1]), length(mutatedPar$exLambda)-1))
        ,exK        = rep(mean(mutatedPar$exK), length(mutatedPar$exK))
        ,exPhiA     = rep(mean(mutatedPar$exPhiA), length(mutatedPar$exPhiA))
        ,exPhiT     = rep(mean(mutatedPar$exPhiT), length(mutatedPar$exPhiT))
        ,exQ        = fittedPar$exQ
        ,pop        = fittedPar$pop
    )

    # as fitted
    meanAge$fitted = fun$getMeanAgeDf(fittedPar,"fitted", prov)

    # silenced everything
    meanAge$silence.all = fun$getMeanAgeDf(silencedPar,"silenced_all", prov)


    # silenced: lambda(t)
    mutatedPar = fittedPar
    mutatedPar$exLambda = silencedPar$exLambda
    meanAge$silence.Lambda = fun$getMeanAgeDf(mutatedPar,"silenced_Lambda", prov)

    # silenced: K(a)
    mutatedPar = fittedPar
    mutatedPar$exK = silencedPar$exK
    meanAge$silence.K = fun$getMeanAgeDf(mutatedPar,"silenced_K", prov)

    # silenced: Phi(a)
    mutatedPar = fittedPar
    mutatedPar$exPhiA = silencedPar$exPhiA
    meanAge$silence.PhiA = fun$getMeanAgeDf(mutatedPar,"silenced_PhiA", prov)

    # silenced: Phi(t)
    mutatedPar = fittedPar
    mutatedPar$exPhiT = silencedPar$exPhiT
    meanAge$silence.PhiT = fun$getMeanAgeDf(mutatedPar,"silenced_PhiT", prov)



    # only: lambda(t)
    mutatedPar = silencedPar
    mutatedPar$exLambda = fittedPar$exLambda
    meanAge$only.Lambda = fun$getMeanAgeDf(mutatedPar,"only_Lambda", prov)

    # only: K(a)
    mutatedPar = silencedPar
    mutatedPar$exK = fittedPar$exK
    meanAge$only.K = fun$getMeanAgeDf(mutatedPar,"only_K", prov)

    # only: Phi(a)
    mutatedPar = silencedPar
    mutatedPar$exPhiA = fittedPar$exPhiA
    meanAge$only.PhiA = fun$getMeanAgeDf(mutatedPar,"only_PhiA", prov)

    # only: Phi(t)
    mutatedPar = silencedPar
    mutatedPar$exPhiT = fittedPar$exPhiT
    meanAge$only.PhiT = fun$getMeanAgeDf(mutatedPar,"only_PhiT", prov)

return(meanAge)
}

fun$getMeanAgeRelativeForPlot = function(ma, ref){
    lapply(ma, function(xma){
        xma %>%
            mutate(
                combined = combined - ma[[ref]]$combined
            )
    })
}

# get mean and 95%CI of the mean age increase
fun$getIncreasePerYear = function(dat){
    fitsum = glm(combined ~ year, data = dat) %>% summary

    list(
        # average annual increase
        rate = round(fitsum$coefficients[2,1], 2)
        # 95% CI
        , CI = round(fitsum$coefficients[2,1] + c(-2,2)*fitsum$coefficients[2,2], 2)
    )
}





    #   Compute/generate summaries

# recalculate if it doesn't exist
meanAgeRDS = paste0(plotdir,"/meanAge.RDS")
if(file.exists(meanAgeRDS)){
    meanAge = readRDS(meanAgeRDS)
} else {
    meanAge = lapply(seq_along(elpd$province), fun$getMeanAgeForPlot)
    saveRDS(meanAge, file=meanAgeRDS)
}
dMeanAge = do.call(rbind, do.call(c, meanAge))


# without any driver: purely demography
fun$getIncreasePerYear(dMeanAge %>% filter(ver == "silenced_all"))

# observed
fun$getIncreasePerYear(dMeanAge %>% filter(ver == "fitted"))
( dMeanAge %>%
        filter(ver == "fitted") %>%
        filter(year < 2000)
    ) %>%
    fun$getIncreasePerYear()
( dMeanAge %>%
        filter(ver == "fitted") %>%
        filter(year >= 2000)
    ) %>%
    fun$getIncreasePerYear()

fun$getIncreasePerYearTable = function(pname, eft){
    d = lapply(meanAge, fun$getMeanAgeRelativeForPlot
            , ref=c('only' = 'silence.all', 'silenced' = 'fitted')[eft]
        ) %>%
        do.call(what=c) %>%
        do.call(what=rbind) %>%
        filter(ver == paste0(eft,"_",pname)
        )

    c(
        'param' = pname
        , 'effectType' = eft
        , 'allTimes' = fun$getIncreasePerYear(d) %>%
            with(paste0(rate, ' (',paste(CI, collapse = ","),')'))
        , 'pre2000' = fun$getIncreasePerYear(d %>% filter(year < 2000)) %>%
            with(paste0(rate, ' (',paste(CI, collapse = ","),')'))
        , '2000onwards' = fun$getIncreasePerYear(d %>% filter(year >= 2000)) %>%
            with(paste0(rate, ' (',paste(CI, collapse = ","),')'))
    )
}


lapply(c('K','Lambda','PhiA','PhiT'), fun$getIncreasePerYearTable, eft = 'only') %>%
do.call(what = rbind)

lapply(c('K','Lambda','PhiA','PhiT'), fun$getIncreasePerYearTable, eft = 'silenced') %>%
do.call(what = rbind)





    #   Plot: variation silencing experiment
    #   (Figure 3)
    
library(latex2exp)

# mean age itself
gFitted = dMeanAge %>%
    filter(ver %in% c("fitted")) %>%
    ggplot(aes( x = year, y = combined, group = province))+
    geom_line(alpha=.2)+
    ylab("Mean age of reported cases")+
    facet_wrap(~ver, nrow=2, labeller = labeller(
        ver = c("fitted"="As fitted", "silenced_all" = "Removed all\nchanges")
    ))+
    theme(
        axis.text.x = element_text(angle=45,hjust=1)
        ,axis.title.x = element_blank()
        ,strip.background = element_rect(fill="black")
        ,strip.text = element_text(color="white")
        ,panel.border = element_rect(fill=NA, size=3)
    )+
    scale_x_continuous(expand=c(0,0))+
    ylim(c(5,35))+
    guides(color=F)

gSilenced = dMeanAge %>%
    filter(ver %in% c("silenced_all")) %>%
    ggplot(aes( x = year, y = combined, group = province))+
    geom_line(alpha=.2)+
    ylab("Mean age of reported cases")+
    facet_wrap(~ver, nrow=2, labeller = labeller(
        ver = c("fitted"="As fitted", "silenced_all" = "Removed all\nchanges")
    ))+
    theme(
        axis.text.x = element_text(angle=45,hjust=1)
        ,axis.title.x = element_blank()
        ,strip.background = element_rect(fill="black")
        ,strip.text = element_text(color="white")
        ,panel.border = element_rect(fill=NA, size=3)
    )+
    scale_x_continuous(expand=c(0,0))+
    ylim(c(5,35))+
    guides(color=F)


# mean age relative to "fitted"
gRefFitted = lapply(meanAge, fun$getMeanAgeRelativeForPlot, ref="fitted") %>%
    do.call(what=c) %>%
    do.call(what=rbind) %>%
    filter(grepl("^silence",ver), ver != "silenced_all") %>%
    mutate(
        ver = gsub("silenced_","", ver)
        , verPanel = factor( ver
            , levels = c('K','Lambda','PhiA','PhiT')
            , labels = c(
                TeX(r'(Without $\kappa(a)$)')
                , TeX(r'(Without $\bar{\tau}(t)$)')
                , TeX(r'(Without $\phi(a)$)')
                , TeX(r'(Without $\phi(t)$)')
            )
        )
    ) %>%
    ggplot(aes( x = year, y = combined))+
    geom_hline(yintercept=0, size=2, lty=1)+
    geom_line(aes(col = ver, group = province), alpha=.4)+
    ylab("Mean age of reported cases\nrelative to fitted")+
    facet_wrap(~verPanel, nrow=1, labeller = label_parsed)+
    theme(
        axis.text.x = element_text(angle=45,hjust=1)
        ,axis.title.x = element_blank()
        ,strip.background = element_rect(fill="#757575")
        ,strip.text = element_text(color="white")
        ,panel.border = element_rect(fill=NA, size=.5, color="#dddddd")
    )+
    scale_x_continuous(expand=c(0,0))+
    scale_color_manual(values = paramCols)+
    guides(color=F)

# mean age relative to "silenced_all"
gRefSilenced = lapply(meanAge, fun$getMeanAgeRelativeForPlot, ref="silence.all") %>%
    do.call(what=c) %>%
    do.call(what=rbind) %>%
    filter(grepl("^only",ver)) %>%
    mutate(
        ver = gsub("only_","", ver)
        , verPanel = factor( ver
            , levels = c('K','Lambda','PhiA','PhiT')
            , labels = c(
                TeX(r'(Solely $\kappa(a)$)')
                , TeX(r'(Solely $\bar{\tau}(t)$)')
                , TeX(r'(Solely $\phi(a)$)')
                , TeX(r'(Solely $\phi(t)$)')
            )
        )
    ) %>%
    ggplot(aes( x = year, y = combined))+
    geom_hline(yintercept=0, size=2, lty=1)+
    geom_line(aes(col = ver, group = province), alpha=.4)+
    ylab("Mean age of reported cases\nrelative to all changes removed")+
    facet_wrap(~verPanel, nrow=1, labeller = label_parsed)+
    theme(
        axis.text.x = element_text(angle=45,hjust=1)
        ,axis.title.x = element_blank()
        ,strip.background = element_rect(fill="#757575")
        ,strip.text = element_text(color="white")
        ,panel.border = element_rect(fill=NA, size=.5, color="#dddddd")
    )+
    scale_x_continuous(expand=c(0,0))+
    scale_color_manual(values = paramCols)+
    guides(color=F)
    
g = ggarrange(
    ggarrange( gFitted, gSilenced
        ,align = "hv"
        ,ncol = 1
        ,nrow = 2
        ,labels = c('a','b')
        ,vjust = 1
    )
    ,ggarrange( gRefFitted,  gRefSilenced
        ,align = "hv"
        ,ncol = 1
        ,nrow = 2
        ,labels = c('c','d')
        ,vjust = 1
    )
    ,ncol = 2
    ,nrow = 1
    ,widths = c(1,4)
    ,align = "h"
)
ggsave( g
    ,filename = paste0(plotdir,"/meanAgePanel_colored.pdf")
    ,width = 7
    ,height = 5
)  




    # Age increase explained by tau(t) + K(a) + Phi(a)
    
meanAge.driver = lapply(seq_along(elpd$province), function(iProv){
    prov = elpd$province[iProv] 

    fittedPar = list(
        exLambda    = elpd$expandedLambda[[iProv]]$mean
        ,exK        = elpd$expandedK[[iProv]]$mean
        ,exPhiA     = elpd$expandedPhiA[[iProv]]$mean
        ,exPhiT     = elpd$expandedPhiT[[iProv]]$mean
        ,exQ        = elpd$Qe[[iProv]]$mean
        ,pop        = dat$pop[[ elpd$province[iProv] ]]
    )
    mutatedPar = fittedPar

    silencedPar = list(
        exLambda    = c(mutatedPar$exLambda[1], rep(mean(mutatedPar$exLambda[-1]), length(mutatedPar$exLambda)-1))
        ,exK        = rep(mean(mutatedPar$exK), length(mutatedPar$exK))
        ,exPhiA     = rep(mean(mutatedPar$exPhiA), length(mutatedPar$exPhiA))
        ,exPhiT     = rep(mean(mutatedPar$exPhiT), length(mutatedPar$exPhiT))
        ,exQ        = fittedPar$exQ
        ,pop        = fittedPar$pop
    )

    mutatedPar = fittedPar
    mutatedPar$exPhiT = silencedPar$exQ
    mutatedPar$exPhiT = silencedPar$exPhiT
    fun$getMeanAgeDf(mutatedPar,"major_drivers", prov)
}) %>%
do.call(what = rbind)

fun$getIncreasePerYear(meanAge.driver)







        #   sliding window averages of marginal FoI
        #   (Figure S7)
        #   .......................................

g = lapply(seq(1,7, by = 2), function(window_size){
    xx = lapply(recon, function(x, wsize = window_size){
            oh = (wsize-1)/2 # overhang
            iYears = (1+oh):(nrow(x)-oh)
            data.frame(
                province = x$province[1]
                , Year = x$Year[iYears]
                , avg = sapply(iYears, function(iYear){
                    x$lambda[iYear + (-oh:oh)] %>% mean
                })
            )
        }) %>%
        do.call(what = rbind) %>%
        mutate(avg = avg * 26) # convert to per year

    rate = lm( avg ~ Year, data = xx %>% filter(Year <= 2005))$coefficients[2] %>%
        round(4) %>%
        {function(x){ -x }}() %>%
        format(digits = 4, scientific = F)
    xx %>%
        ggplot(aes(x = Year, y = avg))+
        geom_line(aes(group = province), alpha = 0.3)+
        geom_smooth(method="gam")+
        annotate("text"
            , x = 2017
            , y = 0.28
            , hjust = 1
            , vjust = 1
            , label = paste("Rate of decline up till 2005:", rate, "per year")
        )+
        scale_y_continuous("Smoothed annual FoI")+
        coord_cartesian(ylim = c(0,0.3))+
        scale_x_continuous(limits = c(1981,2017))+
        ggtitle(paste("Sliding window size:", window_size, "year"))
})
g = ggarrange( plotlist = g
    , ncol = 2
    , nrow = 2
    , align = "hv"
)
ggsave(g
    ,filename = paste0(plotdir,"/lambdaDeclineRate.pdf")
    ,width = 8
    ,height = 6
)


# 95%IQR for the decline in FoI up to 2005
sapply(recon, function(xx){
    glm( lambda ~ Year, data = xx %>% 
        filter(Year <= 2005) %>%
        mutate(lambda = lambda * 26)
    )$coefficients[2]
}) %>%
{function(x){ -x }}() %>%
quantile(probs = c(0.025, 0.5, 0.975)) %>%
round(3) %>%
format(digits = 4, scientific = F)




    #   Mechanistic explanation of FoI decline
    #   ......................................
    setPlotdir('reconstruction')

# assess the trends in beta

    # countrywide
    # : up to 1990
    fit = glm( Beta ~ Year
        , data = recon %>% 
            do.call(what = rbind) %>% 
            filter(Year <= 1990)
    ) %>% summary
    round(coef(fit)[2,1] + c(-2,0,2) * coef(fit)[2,2], 2)

    # : post-1990
    fit = glm( Beta ~ Year
        , data = recon %>% 
            do.call(what = rbind) %>% 
            filter(Year > 1990)
    ) %>% summary
    round(coef(fit)[2,1] + c(-2,0,2) * coef(fit)[2,2], 2)

    # by province
    # : up to 1990
    lapply(recon, function(x){
        fit = glm( Beta ~ Year
            , data = x %>% 
                filter(Year <= 1990)
        ) %>% summary
        out = round(coef(fit)[2,1] + c(-2,0,2) * coef(fit)[2,2], 2)
        c( out
            , include.zero = between(0, out[1], out[2])
            , below.zero = out[2] < 0
            , above.zero = out[1] > 0
        )
    }) %>%
    do.call(what = rbind) %>%
    summary

    # : post-1990
    lapply(recon, function(x){
        fit = glm( Beta ~ Year
            , data = x %>% 
                filter(Year > 1990)
        ) %>% summary
        out = round(coef(fit)[2,1] + c(-2,0,2) * coef(fit)[2,2], 2)
        c( out
            , include.zero = between(0, out[1], out[2])
            , below.zero = out[2] < 0
            , above.zero = out[1] > 0
        )
    }) %>%
    do.call(what = rbind) %>%
    summary



# Delineate the contributions of Beta and I, on lambda
lambdaContributions = function(yearRange){
    lapply(recon, function(x){
            fit = glm( Beta ~ Year
                , data = x %>% filter(between(Year, yearRange[1], yearRange[2]))
            ) %>%
            predict(newdata = data.frame(Year = yearRange))
            fold.change.Beta = (fit[2]-fit[1]) / fit[1]
            fold.Beta = fit[2]/fit[1]

            fit = glm( lambda ~ Year
                , data = x %>% filter(between(Year, yearRange[1], yearRange[2]))
            ) %>%
            predict(newdata = data.frame(Year = yearRange))
            fold.change.lambda = (fit[2]-fit[1]) / fit[1]
            fold.lambda = fit[2]/fit[1]

            data.frame(province = x$province[1]
                , fold.change.lambda
                , fold.change.Beta
                , fold.lambda
                , fold.Beta
            )
        }) %>%
        do.call(what = rbind) %>%
        mutate(prop.Beta = fold.change.Beta / fold.change.lambda) 
}

dFold = lambdaContributions(c(1981,2017))

dFold %>%
    select(-province) %>%
    lapply(quantile, probs = c(0.025, 0.5, 0.975))

# number of provinces with beta in the same direction as lambda
with(dFold, sign(fold.change.lambda) == sign(fold.change.Beta)) %>% sum


# Figure S9: plot in space
sp$map = sp$provincePolygons %>%
    left_join(sp$provinceCentroids, by = 'id') %>%
    rename(province = province72)

plotMap = function(xvar, fill.label, dFold){
    sp$map %>%
        left_join(dFold
            , by = "province"
        ) %>%
        mutate_at(xvar, function(x) x*100) %>%
        ggplot(aes(long, lat, group = group))+
        geom_polygon(aes_string(fill = xvar), color = "white", size = 0.3)+
        # scale_fill_gradient2(fill.label
        #     , low = "blue"
        #     , high = "red"
        #     , mid = "#dddddd"
        #     , midpoint = 0
        # )+
        scale_fill_gradientn(fill.label
            , colors = c("blue", "#dddddd", "red")
            , limits = c(-100, 100)
        )+

        coord_quickmap()+
        theme_void()+
        theme(
            legend.position = 'bottom'
            , legend.title = element_text(hjust = 1, vjust = 0.5, margin = margin(r = 10))
        )+
        facet_grid(. ~ Years)
}

dFold = lapply(list(
        c(1981, 2017)
        , c(1981, 1990)
        , c(1991, 2017)
    ), function(yearRange){
        lambdaContributions(yearRange) %>%
        mutate(Years = paste(yearRange, collapse="-"))
    }) %>%
    do.call(what = rbind) %>%
    mutate(Years = factor(Years, levels = unique(Years)))

gMap = ggarrange(plotlist = 
    mapply( plotMap
        , xvar = c(
            'fold.change.lambda'
            , 'fold.change.Beta'
        )
        , fill.label = c(
            'Percentage change in\nforce of infection (FoI)'
            , 'Percentage change in\ntransmission efficiency'
        )
        , dFold = list(dFold)
        , SIMPLIFY = FALSE
    )
    , nrow = 2
    , ncol = 1
)

g = ggarrange(
    dFold %>%
        ggplot(aes(x = fold.change.Beta * 100, y = fold.change.lambda * 100))+
        geom_vline(xintercept = 0, linetype = 2, size = 0.1)+
        geom_hline(yintercept = 0, linetype = 2, size = 0.1)+
        geom_abline(slope = 1, linetype = 2, size = 0.2, color = 'red')+
        geom_point(shape = 1, size = 1.5)+
        theme_classic()+
        ylab('Percentage change in FoI')+
        xlab('Percentage change in\ntransmission efficiency')+
        coord_fixed(ratio = 1)+
        facet_grid(. ~ Years)
    , gMap
    , ncol = 1
    , nrow = 2
    , heights = c(1,2)
    , labels = 'auto'
)
ggsave(g
    , file = file.path(outdir,'foldChange.pdf')
    , width = 7
    , height = 12.5
)














    # Beta

# Median, 95%IQR of Beta
recon %>%
    do.call(what = rbind) %>%
    filter(Year %in% seq(1981, 2017, by = 9)) %>%
    group_by(Year) %>%
    summarize(
        Median = median(Beta)
        , lower = quantile(Beta, 0.025)
        , upper = quantile(Beta, 0.975)
    )
    
# fit and summarize the linear regressions
fitsum.Beta = lapply(recon, function(x){
    lm(log2(Beta) ~ Year, data = x) %>% summary
})


    # I
    
fitsums = lapply(recon, function(x){
    x = x %>%
        mutate(
            infected = I.x + I.z1 + I.z2 + I.z3
            , popsize = S.x + S.z1 + S.z2 + S.z3 + S.z4
            , prop = infected/popsize
        ) %>% 
        filter(between(Year, 1988, 2004))
    lm( prop ~ Year, data = x) %>% summary
})
# rate of change excludes zero
sapply(fitsums, function(x){
    x$coefficients[2,4] <= 0.05
}) %>% summary

# rate of change
I.changerate = sapply(fitsums, function(x){
    x$coefficients[2,1]
}) 
I.changerate %>% quantile(c(0.025, 0.5, 0.975))


