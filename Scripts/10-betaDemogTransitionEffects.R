#!/usr/bin/env Rscript --vanilla
rm(list=ls())


outdir = "10-demogTransitionEffects"
dir.create(outdir, recursive=T)


library(deSolve)

library(dplyr)
library(tidyr)
library(ggplot2); theme_set(theme_bw())
library(ggpubr)

# define population parameters
popsize = 10^6     # starting population size

t.constant = 400    # length of time in which popsize is constant
t.transit = 240     # length transition time
t.transit.lag = 170 # time lag between death rate decline to birth rate decline
t.total = 900       # total time



    # define time points
    times = seq(1/26,t.total,by=1/26)
    # define ages available in the population
    ages = seq(0,100-1/26, by=1/26)


getPop = function(
    mu      # initial birth rate, death rate
    ,mu2    # birth/death at new equilibrium
  ){

    # convert from yearly to biweekly timescale
    mu = mu/26
    mu2 = mu2/26
  
    muDeath = c(
        rep(mu, t.constant*26)
        ,seq(mu, mu2, length.out = t.transit*26)
        ,rep(mu2, (t.total - t.constant - t.transit) * 26)
    )
    muBirth = c(
        rep(mu, (t.constant + t.transit.lag)*26)   # birth doesn't decline until later
        ,seq(mu, mu2, length.out = t.transit*26)
        ,rep(mu2, (t.total - t.constant - t.transit.lag - t.transit)*26)
    )[1:(t.total*26)]
    
    # get total population sizes
    popsize.t = matrix(NA, length(times), 6)
    for(it in seq_along(times)){
        time = times[it]
        dnBirth = muBirth[time*26]
        dnDeath = muDeath[time*26]
        
        pop1 = ifelse(
            it==1
            ,popsize
            ,popsize.t[it-1,2]
        )

        out = c(
            time = time
            ,popsize = pop1 * (1 + dnBirth - dnDeath)
            ,nBirth = dnBirth * pop1
            ,nDeath = dnDeath * pop1
            ,muBirth = dnBirth
            ,muDeath = dnDeath
        )
        popsize.t[it, ] = out
    }
    colnames(popsize.t) = names(out)
    plot(popsize.t[,'popsize'], type='l')

    return(popsize.t %>% as.data.frame)
}

popsize.t = getPop(1/40,1/80)




    #   Mechanistic simulations of lambda
    #   .................................

# define transmission parameters at introduction
t.intro = 100       # introduce dengue long before demographic transient

getLambdaTs = function(Beta = 2.3, t.lambda.fix = Inf){   # use value ~ 1 to 10

    S = matrix(c(1,0,0,0,0), 5, length(ages)) %>% t
    I = 1/popsize.t[[t.intro*26 ,"popsize"]]

    timepoints = times[-(1:(26*t.intro - 1))] 
    xx = matrix(NA, length(timepoints), 11)

    for(it in seq_along(timepoints)){

            t = timepoints[it]

            # pop count of each cohort at time t
            popAge = exp( -
                # sum up death rates encountered by each cohort
                c(0, popsize.t[t*26 - (2:2600) + 1, "muDeath"] %>% cumsum)
            ) * popsize.t[ t*26 - (1:2600) + 1, "nBirth" ]

            
            if(t < t.lambda.fix){
                # update lambda based on I(t)            
                lambda = Beta * I/4 # per serotype lambda (there are 4 serotypes)
            } else {
                # fix the lambda and find the beta that 
                # would keep the transmission intensity the same
                Beta = lambda * 4 / I
            }
            
            # log the FoI
            foi = sum((S %*% ((4:0)*lambda)) * popAge) / sum(popAge)
            
            p = exp( -lambda )    # proportion that survived a particular serotype at this timepoint
            # fraction of individuals that got infected in this time step
            # : age x 4
            p = sapply(4:1, function(i){
                # for an individual still susceptible to i serotypes
                i * (1-p)^(1) * (p)^(i-1)
            })
            I.age = S[ ,-5] %*% diag(p)

            # update I to infect people at next time step
            I = sum(I.age * popAge)/sum(popAge)

            # update S due to infections at time t
            S = S +
                 cbind( 0, I.age) -         # inflow
                 cbind(I.age, 0)            # outflow            

            # update S due to birth and aging
            S = rbind(c(1,0,0,0,0), S[-nrow(S), ])

        out = c(
            "time" = t
            ,"I" = I
            ,"lambda" = lambda
            ,"beta" = Beta
            ,"foi" = foi
            ,"meanAge" = colSums(I.age * popAge * ages)/colSums(I.age * popAge)
            ,"meanAge.combined" = sum(I.age * popAge * ages)/sum(I.age * popAge)
            ,"meanAge.pop" = sum(popAge * ages)/sum(popAge)
        )
        xx[it, ] = out
    }
    colnames(xx) = names(out)
    as.data.frame(xx) %>%
        mutate(Year = time - t.constant + 1800) %>%
        filter(Year >= year.plot.start) 
}

plotAgeTs = function(xxZoom){

    g = ggplot( mapping = aes(x = Year)
        )+
        geom_vline(data = pop.points, aes(xintercept = Year), color = "#dddddd", lty = 2)+
        geom_line(data = xxZoom %>%
            select(Year, meanAge1:meanAge4) %>%
            gather(infection, meanAge, -Year)
            , aes( y = meanAge, col = infection )
            , size = 0.5
        )+
        geom_line(data = xxZoom
            , aes( y = meanAge.combined )
            , size = 1.2
            , col = "black"
        )+
        geom_line(data = xxZoom
            , aes( y = meanAge.pop )
            , size = 1.2
            , lty = "22"
            , col = "grey"
        )+
        guides(color=F)+
        geom_text( data = xxZoom %>% 
            filter(Year == min(Year)) %>%
            select(meanAge1:meanAge.pop, Year) %>%
            gather(infection,age,-Year)
            ,aes( x = Year + 20, y = age + 0.2, label = infection)
            ,label = c(
                paste("Infection",1:4)
                ,"All infections"
                ,"Population"
            )
            ,hjust = 0
            ,vjust = 0
            ,color = c(infectionCols,"black","grey")
        )+
        scale_x_continuous(expand=c(0,0), position="top") +
        scale_color_manual( values = infectionCols)+
        theme(
            axis.title.x = element_blank()
            ,panel.grid = element_blank()
            ,plot.margin = unit(c(0,1,1,1) ,"lines")
        )+
        ylab("Mean Age (years)")
    return(g)
}

plotTs = function(xxZoom, ts.param, fill.color, y.label){
    xxZoom %>%
    ggplot(aes_(x = ~Year, ymax = as.name(ts.param), ymin = 0))+
        geom_ribbon(fill = fill.color)+
        geom_vline(data = pop.points, aes(xintercept = Year), col = "white")+
        scale_x_continuous(expand=c(0,0))+
        ylab(y.label)+
        theme(
            axis.text.x = element_blank()
            ,axis.title.x = element_blank()
            ,panel.grid = element_blank()
            ,plot.margin = unit(c(0,1,0,1) ,"lines")
        )
}

infectionCols = c(
    "#edd312"
    ,"#ed5412"
    ,"#4a7500"
    ,"#6a12ed"
)
simCols = c(
    fixedBeta = "#363636"
    , fixedLambda = "#ff9063"
)

pop.points = data.frame(
        time = c(t.constant
            , t.constant + t.transit
            , t.constant + t.transit.lag
            , t.constant + t.transit.lag + t.transit
        )
        ,event = c(
            "Death starts\ndeclining"
            ,"Death stops\ndeclining"
            ,"Birth starts\ndeclining"
            ,"Birth stops\ndeclining"
        )
    ) %>%
    left_join( popsize.t %>%
            mutate(time = round(time,2))
        ,by = "time"
    ) %>%
    mutate(Year = time + 1800 - t.constant)


# simulate
year.plot.start = 1700
sim.fixedBeta = getLambdaTs()
sim.fixedLambda = getLambdaTs( t.lambda.fix = t.constant + t.transit.lag)



# plot
library(latex2exp)
g.fixedBeta = ggarrange(
    # population size
    ggplot( popsize.t %>% 
            mutate(
                Year = time + 1800 - t.constant
            ) %>%
            filter(Year >= year.plot.start)
        , aes(x = Year, y = popsize/(10^6)))+
        geom_ribbon(aes(
            ymax = popsize/(10^6)
            ,ymin = 0
        ), alpha = .3)+
        geom_vline(data = pop.points, aes(xintercept = Year), col = "white")+
        geom_point( data = pop.points
            , size = 3
        )+
        geom_text( data = pop.points %>% filter(time < max(time))
            , aes(
                label = event
                , x = Year - 40
                , y = popsize/(10^6) + .1
            )
            , angle = 90
            , hjust = 0
            , lineheight = 0.8
        )+
        geom_text( data = pop.points %>% filter(time == max(time))
            , aes(
                label = event
                , x = Year - 40
                , y = popsize/(10^6) - .1
            )
            , angle = 90
            , hjust = 1
            , lineheight = 0.8
        )+
        scale_x_continuous(
            expand=c(0,0)
        )+
        scale_y_continuous(
            "Population size\n(million)"
            ,expand=c(0,.15)
        )+
        theme(
            axis.text.x = element_blank()
            ,axis.title.x = element_blank()
            ,panel.grid = element_blank()
            ,plot.margin = unit(c(1,1,0,1) ,"lines")
        )

    # Beta
    , plotTs(sim.fixedBeta
        , "beta", simCols['fixedBeta']
        , TeX(r'(Biweekly $\beta(t)$)')
    )

    # Lambda
    , plotTs(sim.fixedBeta
        , "lambda", simCols['fixedBeta']
        , TeX(r'(\overset{Biweekly $\bar{\tau}(t)$,}{per serotype})')
    )
        
    # FOI
    , plotTs(sim.fixedBeta, "foi", simCols['fixedBeta'], "Biweekly FoI")

    # Age
    , plotAgeTs(sim.fixedBeta)
    
    ,nrow = 5
    ,ncol = 1
    ,heights = c(1,0.7,0.7,0.7,1.3)
    ,align = "v"
    ,common.legend = T
    ,labels = c('a','b','','')
)




g.fixedLambda = ggarrange(
    ggplot(mapping = aes(x = Year, ymax = I, ymin = 0))+
        geom_ribbon(data = sim.fixedLambda, fill = "#ff9063")+
        geom_ribbon(data = sim.fixedBeta, fill = "#363636")+
        geom_vline(data = pop.points, aes(xintercept = Year), col = "white")+
        scale_x_continuous(expand=c(0,0))+
        ylab("Infectious fraction")+
        theme(
            axis.text.x = element_blank()
            ,axis.title.x = element_blank()
            ,panel.grid = element_blank()
            ,plot.margin = unit(c(1,1,0,1) ,"lines")
        )

    # Beta
    , plotTs(sim.fixedLambda
        , "beta", simCols['fixedLambda']
        , TeX(r'(Biweekly $\beta(t)$)')
    )
        
    # Lambda
    , plotTs(sim.fixedLambda
        , "lambda", simCols['fixedLambda']
        , TeX(r'(\overset{Biweekly $\bar{\tau}(t)$,}{per serotype})')
    )
        
    # FOI
    , plotTs(sim.fixedLambda, "foi", simCols['fixedLambda'], "Biweekly FoI")


    ,plotAgeTs(sim.fixedLambda)
    ,nrow = 5
    ,ncol = 1
    ,heights = c(1,0.7,0.7,0.7,1.3)
    ,align = "v"
    ,common.legend = T
    ,labels = c('c','d','','')
)


g = ggarrange( g.fixedBeta, g.fixedLambda
    , nrow = 1
    , ncol = 2
)
ggsave( g
    ,filename = paste0(outdir,"/lambda-beta-meanAge.pdf")
    ,width = 8
    ,height = 9
)

