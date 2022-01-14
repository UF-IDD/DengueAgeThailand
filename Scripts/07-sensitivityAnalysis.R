
# Assess impacts of
# 1) Informative vs wide priors on Q(i,t)
# 2) Less degrees of freedom on tau(t): 37 vs 16


library(tidyverse); theme_set(theme_classic())
library(ggpubr)


outdir = "07-sensitivityAnalysis"
dir.create(outdir, recursive = T)



    #   1) Informative vs wide priors on Q(i,t)
    #   ............................

Params = c(
    'lambda' = 'Per-serotype\ninfection intensity'
    , 'binnedK' = 'Age-modifier for\ninfection intensity'
    , 'binnedPhiT' = 'Time-varying\nreporting rate'
    , 'binnedPhiA' = 'Age-specific\nreporting rate'
    , 'binnedQe' = 'Clinical\ndetectability'
)

dat = lapply(c("_widePrior", ""), function(ver){
    infiles = paste0(
            "06-assessFits"
            , ver
            , "/DHF_add/summary_converged/m3h.rtaQeT"        
        ) %>%
        list.files(pattern = "\\.csv", full.names = T)
    infiles %>%
        lapply(function(f){
            read.csv(f, stringsAsFactors = F) %>%
            mutate(province = gsub('\\.csv$', '', basename(f)))
        }) %>%
        do.call(what = rbind) %>%
        
        filter(var %in% names(Params)) %>%
        select(mean, lower, upper, var, ivar, province) %>%
        rename_at(vars(mean, lower, upper), function(x){
            paste0(x, ver)
        })
}) %>%
plyr::join_all(by = c('province','var', 'ivar'))

Epoches = paste0(
    seq(1981,2012,by=9)
    ,"-"
    ,c(seq(1981+8,2012,by=9),2017)
)

g = ggarrange(
    dat %>%
        filter(var != 'binnedQe') %>%
        mutate(
            var = factor(var, levels = names(Params), labels = Params)
            , var = droplevels(var)
        ) %>%
        # mutate(var = Params[var]) %>%
        ggplot(aes(x = mean, y = mean_widePrior))+
        geom_abline(slope = 1, color = 'red', size = 0.2)+
        geom_point(shape = 1, stroke = 0.3)+
        facet_wrap(~ var, nrow = 1)+
        coord_fixed(ratio = 1)+
        theme_bw()+
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
        xlab('Informative priors')+
        ylab('Weakly informative priors')
    , dat %>%
        filter(var == 'binnedQe') %>%
        mutate(
            epoch = ceiling(ivar/4)
            , epoch = factor(epoch, labels = Epoches)
            , ivar = (ivar %% 4)
            , ivar = ivar + (ivar==0)*4
            , ivar = paste('Infection', ivar)
        ) %>%
        ggplot(aes(x = mean, y = mean_widePrior))+
        geom_abline(slope = 1, color = 'red', size = 0.2)+
        geom_point(shape = 1, stroke = 0.3)+
        facet_grid(ivar ~ epoch)+
        coord_fixed(ratio = 1)+
        theme_bw()+
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
        xlab('Informative priors')+
        ylab('Weakly informative priors')
    , ncol = 1
    , nrow = 2
    , heights = c(2,4)
    , labels = 'auto'
)
ggsave(g
    , filename = file.path(outdir, 'priors_effects.pdf')
    , width = 6
    , height = 8
)





    #   2) Less degrees of freedom on tau(t): 37 vs 16
    #   ............................

dat = lapply(c("_16tau", ""), function(ver){
    infiles = paste0(
            "06-assessFits"
            , ver
            , "/DHF_add/summary_converged/m3h.rtaQeT"        
        ) %>%
        list.files(pattern = "\\.csv", full.names = T)
    dat = infiles %>%
        lapply(function(f){
            read.csv(f, stringsAsFactors = F) %>%
            mutate(province = gsub('\\.csv$', '', basename(f)))
        }) %>%
        do.call(what = rbind) %>%
        filter(var %in% names(Params)) %>%
        select(mean, lower, upper, var, ivar, province) %>%
        rename_at(vars(mean, lower, upper), function(x){
            paste0(x, ver)
        })
    if(ver == "_16tau"){
        out = 1
        out = append(out, rep((1:8) + length(unique(out)), each = 2))
        out = append(out, rep((1:7) + length(unique(out)), each = 3))
        out = split(seq_along(out), out)
        lambda.scaler = 3
        dat = dat %>%
            filter(var == 'lambda') %>%
            mutate(ivar = out[ivar]) %>%
            unnest(cols = ivar) %>%
            mutate(mean_16tau = mean_16tau * lambda.scaler) %>%
            rbind(dat %>%
                filter(var != 'lambda') %>%
                mutate(mean_16tau = mean_16tau / ifelse(var=='binnedK', lambda.scaler, 1))
            )
    }
    return(dat)
}) %>%
plyr::join_all(by = c('province','var', 'ivar'))

g = dat %>%
    mutate(
        var = factor(var, levels = names(Params), labels = Params)
        , var = droplevels(var)
    ) %>%
    ggplot(aes(x = mean, y = mean_16tau))+
    geom_abline(slope = 1, color = 'red', size = 0.2)+
    geom_point(shape = 1, stroke = 0.3)+
    facet_wrap(~ var, nrow = 1)+
    coord_fixed(ratio = 1)+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
    xlab('Annual per-serotype infection intensity')+
    ylab('16 bins of per-serotype infection intensity')
ggsave(g
    , filename = file.path(outdir, 'tauBinning.pdf')
    , width = 8
    , height = 5
)




# performance: compare ELPD
library(loo)
getPath = function(ver = ""){
    paste0("06-assessFits", ver,"/DHF_add/elpd_converged/loo")
}
# for each province
Loo = lapply(list.files(getPath()), function(f){
    Vers = c("Annual" = "", "Binned" = "_16tau")
    Loo = lapply(Vers, function(ver){
        Loo = file.path(getPath(ver), f) %>% readRDS
        meta = strsplit(names(Loo),"/") %>%
            do.call(what=rbind)
        names(Loo) = meta[ ,3]
        Loo
    })
    # for each model
    lapply(names(Loo[[1]]), function(model.name){
        x = lapply(Loo, function(x) x[[model.name]])
        x = x[!sapply(x, is.null)]
        if(length(x)==0){
            return(NULL)
        } else if(length(x)<2){
            out = data.frame(
                elpd_diff = 0
                ,se_diff = 0
            )
            rownames(out) = names(x)
        } else {
            out = loo_compare(x = x) %>%
                as.data.frame %>%
                select(elpd_diff, se_diff)
        }
        out$rank = with(out, (se_diff * 2) <= -elpd_diff ) %>% cumsum
        out %>%
            mutate(
                tauVer = rownames(out)
                , province = gsub('\\.RDS$', '', f)
                , model = model.name
            )
    }) %>%
    do.call(what = rbind)
}) %>%
do.call(what = rbind)

# How often does binnedTau has non-inferior fit?
convergence = by(Loo, Loo$tauVer, function(x){
        x %>%
            group_by(model) %>%
            summarize(
                Best = sum(rank==1)
                , Converged = n()
            )
    }) %>%
    do.call(what = cbind) %>%
    select(-Binned.model) %>%
    rename(Model = Annual.model) %>%
    arrange(Binned.Converged)

convergence %>%
    mutate(delta = Binned.Converged - Annual.Converged) %>%
    summarize(
        median(delta)
        , quantile(delta, 0.025)
        , quantile(delta, 0.975)
        , sum(delta < 0)
    )


convergence %>%
    write_csv(file.path(outdir, 'tauBinning.csv'))

