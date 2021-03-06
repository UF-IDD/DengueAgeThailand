#!/usr/bin/env Rscript --vanilla

library(argparse)
library(tidyverse)
library(rstan)

parser = ArgumentParser(description = "Plot priors against posteriors.")


parser$add_argument('dataver', default = "DHF_add")
parser$add_argument('-stanroot', default = "05-estimates")
parser$add_argument('-od', default = "06-posteriors")
parser$add_argument('-model.include', nargs = "*", default = "m3h.rtaQeT")

# parse arguments
inputArg = parser$parse_args()
# inputArg = parser$parse_args('DHF_add')


# function to set plot directory
outdir = with(inputArg, file.path(od, dataver))
setPlotdir = function(...){
    plotdir <<- file.path(outdir, ...)
    message(paste("Set plot directory:",plotdir))
    dir.create(plotdir, recursive=T)
}
# create root for outputs
setPlotdir()



#   For each province,
#   1) import fits
#   2) summarize posterior by parameter
x = list.files( with(inputArg, file.path(stanroot, dataver))
        , pattern = "\\.RDS$"
        , recursive=T
    ) %>%
    strsplit("/") %>%
    do.call(what=rbind) %>%
    as.data.frame
    colnames(x) = c('model','province')
x$province = gsub('\\.RDS$','',x$province)
x$rds = list.files( with(inputArg, file.path(stanroot, dataver))
    , pattern = "\\.RDS$"
    , recursive=T
    , full.names = T
)

x = x %>% filter(model %in% inputArg$model.include)





    #   Plot density of Phi(a), Phi(t), one panel per age/time
    #   ...............................

ageBinLabels = function(x = seq(0, 60, by = 4)){
    paste0(x, c(paste0('-',x[-1]-1), "+"))
}
timeBinLabels = function(x = seq(1981, 2011, by = 2)){
    paste0(x, c(paste0('-',x[-1]-1), "+"))
}


plotDensityPanel = function(varname, xval, yval.func, yLimit, Title, panel.labels){
    # posterior
    # : process one province at a time
    xx = lapply(x$rds, function(rds){
        readRDS(rds) %>%
            as.matrix(pars = varname) %>%
            apply(2, density)
    })
    par(oma = c(1,1,1.5,0), mar = c(2,2,3,1), mfrow = c(4,4))
    lapply(seq_along(xx[[1]]), function(i){
        plot(x=range(xval), y = c(0,yLimit), type = 'n'
            , xlab = NA
            , ylab = NA
            , main = panel.labels[i]
        )
        lapply(xx, function(xxprov){
            with(xxprov[[i]], lines(x,y, lwd = 0.5, col = "#00000020"))
        })
        # prior
        # lines(xval, yval.func(xval), col = 'blue')
    }) %>%
    invisible
    mtext(Title, side = 3, line = 0, outer = TRUE)
    mtext("Density", side = 2, line = 0, outer = TRUE, cex = 0.9)
    mtext("Value", side = 1, line = 0, outer = TRUE, cex = 0.9)
}


pdf(file = file.path(plotdir, "PhiT.pdf")
    , width = 6
    , height = 6
)
    plotDensityPanel("binnedPhiT"
        , xval = seq(0, 1, by = 0.01)
        , yval.func = function(x){ dnorm(x, 0.5, 1) }
        , yLimit = 8
        , Title = "Time-varying reporting rate"
        , panel.labels = timeBinLabels()
    )
dev.off()


pdf(file = file.path(plotdir, "PhiA.pdf")
    , width = 6
    , height = 6
)
    plotDensityPanel("binnedPhiA"
        , xval = seq(0, 1, by = 0.01)
        , yval.func = function(x){ dnorm(x, 0.5, 1) }
        , yLimit = 10
        , Title = "Age-specific reporting rate"
        , panel.labels = paste('Age', ageBinLabels())
    )
dev.off()


q(save = 'n')



    #   One parameter per panel
    #   .......................

plotDensity = function(varname, xval, yval.func, yLimit, Title){
    plot(x=range(xval), y = c(0,yLimit), type = 'n'
        , xlab = "Value"
        , ylab = "Density"
        , main = Title
    )
    # posterior
    # : process one province at a time
    lapply(x$rds, function(rds){
        xx = readRDS(rds) %>%
            as.matrix(pars = varname) %>%
            apply(2, density)
        lapply(xx, function(a){
            with(a, lines(x,y, lwd = 0.5, col = "#00000020"))
        }) %>%
        invisible
    }) %>%
    invisible
    # prior
    lines(xval, yval.func(xval), col = 'blue')
}



pdf(file = file.path(plotdir, "Param.pdf")
    , width = 10
    , height = 5
)
par(mfrow = c(1,4))
    plotDensity("lambda"
        , xval = seq(0, 1, by = 0.01)
        , yval.func = function(x){ dnorm(x, 0.05, 0.5) }
        , yLimit = 20
        , Title = "Time-varying\ninfection intensity"
    )

    plotDensity("binnedK"
        , xval = seq(0, 1, by = 0.01)
        , yval.func = function(x){ dnorm(x, 1, 3) }
        , yLimit = 20
        , Title = "Age-specific modifer\nfor infection intensity"
    )

    plotDensity("binnedPhiT"
        , xval = seq(0, 1, by = 0.01)
        , yval.func = function(x){ dnorm(x, 0.5, 1) }
        , yLimit = 8
        , Title = "Time-varying\nreporting rate"
    )

    plotDensity("binnedPhiA"
        , xval = seq(0, 1, by = 0.01)
        , yval.func = function(x){ dnorm(x, 0.5, 1) }
        , yLimit = 10
        , Title = "Age-specific\nreporting rate"
    )
dev.off()
 


plotDensity.Qe = function(qi, xval = seq(0, 1, by = 0.01), yval.func, yLimit){
    plot(x=range(xval), y = c(0,yLimit), type = 'n'
        , xlab = "Value"
        , ylab = "Density"
        , main = paste("Clinical detectability\nof infection", qi)
    )
    # posterior
    # : process one province at a time
    lapply(x$rds, function(rds){
        xx = (readRDS(rds) %>%
            as.matrix(pars = "binnedQe"))[ , (1:4) + (qi-1)*4] %>%
            apply(2, density)
        lapply(xx, function(a){
            with(a, lines(x,y, col = "#00000020"))
        }) %>%
        invisible
    }) %>%
    invisible
    # prior
    lines(xval, yval.func(xval), col = 'blue')
}


pdf(file = file.path(plotdir, "Qe.pdf")
    , width = 10
    , height = 3
)
par(mfrow = c(1,4))
    plotDensity.Qe(1
        , yval.func = function(x){ dbeta(x, 1, 9) }
        , yLimit = 35
    )
    plotDensity.Qe(2
        , yval.func = function(x){ dbeta(x, 19, 1) }
        , yLimit = 20
    )
    plotDensity.Qe(3
        , yval.func = function(x){ dbeta(x, 1, 19) }
        , yLimit = 35
    )
    plotDensity.Qe(4
        , yval.func = function(x){ dbeta(x, 1, 19) }
        , yLimit = 35
    )
dev.off()


