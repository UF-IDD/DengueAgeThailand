
library(tidyverse); theme_set(theme_bw())
library(ggpubr)

datadir = '04-dataForStan'
casever = 'DHF_add'
config = readRDS( paste0(datadir,"/config.RDS") )

plotdir = '04-describeData'
dir.create(plotdir, recursive = T)

    #   Read in data
    #   ............
    
fun = new.env()
fun$readDataCsv = function(x){
    as.matrix(read.csv(x, stringsAsFactors=F))
}

places = list.files(
    file.path(datadir,casever)
    , pattern = '\\.csv$'
) %>%
gsub(pattern = '\\.csv$', replacement = '')

case = file.path(datadir,casever, paste0(places,'.csv')) %>%
    lapply(fun$readDataCsv)
pop = file.path(datadir,'pop', paste0(places,'.csv')) %>%
    lapply(fun$readDataCsv)




    #   Plot time-series
    #   ................

# start and end of age groups
ageBound = with(config, {
    apply(ageGroups, 1, function(x){
        x = ages[x==1]
        c(min(x), max(x)+1)
    })
})

# midpoints of the age strata
ageMid = colMeans(ageBound)

# age of cohort at each time point in the time-series
cohortAge = with(config, outer(cohorts, times, function(h,t) t-h))


Ts = mapply( function(prov,ca,pa){
        data.frame(
            province = prov
            , year = config$times
            , meanPopAge = pa
            , meanCaseAge = ca
        ) %>%
        gather(ageType, meanAge, -province, -year)
    }
    , prov = places
    , ca = lapply(case, function(casep){
        apply(casep, 2, function(x){
            sum((x * ageMid)[x != -1]) / sum(x[x != -1])
        })
    })
    , pa = lapply(pop, function(popp){
        include = apply(cohortAge, 2, function(x) x %in% config$ages)
        colSums(popp * cohortAge * include) /
        colSums(popp * include)
    })
    , SIMPLIFY = FALSE
    , USE.NAMES = FALSE
) %>%
do.call(what = rbind)


gTs = Ts %>%
    ggplot(aes(x = year, y = meanAge, col = ageType))+
    geom_line(aes(group=paste(province,ageType)), alpha=.1)+
    geom_point(alpha=.1, size=.5)+
    theme_minimal()+
    theme(
        axis.text.x = element_text(angle=0, hjust=.5)
        ,legend.position = c(0.8, 0.15)
        ,legend.title = element_blank()
    )+
    scale_x_continuous('Year', expand=c(.01,.01))+
    scale_y_continuous('Mean age (yrs)')+
    scale_color_manual(
        values = c("red","black")
        ,labels = c("DHF cases","Population")
    )+
    guides(color = guide_legend(override.aes = list(size = 2, alpha=1)))
    


# case counts
dCase = lapply(case, function(casep){
    mapply( function(x,y){ x %>% mutate(year = y) }
        , x = apply(casep, 2, function(x){
            out = as.data.frame(t(ageBound[ ,x != -1]))
            colnames(out) = c('ageStart','ageEnd')
            out %>%
                mutate(count = x[x != -1])        
        })
        , y = config$times
        , SIMPLIFY = FALSE
        , USE.NAMES = FALSE
    ) %>%
    do.call(what = rbind)
}) %>%
Reduce(f=rbind) %>%
group_by(ageStart, ageEnd, year) %>%
summarize(count = sum(count)) %>%
group_by(year) %>%
mutate(denseCase = (count/(ageEnd-ageStart))/sum(count))

# popsize
dPop = lapply(pop, function(popp){
    include = apply(cohortAge, 2, function(x) x %in% config$ages)
    lapply(1:ncol(include), function(T){
        data.frame(
            year = config$times[T]
            , popsize = popp[include[ ,T], T]
            , age = cohortAge[include[ ,T], T]
        )
    }) %>%
    Reduce(f=rbind)
}) %>%
Reduce(f=rbind) %>%
group_by(age, year) %>%
summarize(popsize = sum(popsize)) 

# density of case age and pop age
ageDist = by(dCase, dCase %>% select(ageStart,ageEnd,year), function(x){
    x %>%
        mutate(popsize = x %>%
            mutate(age = mapply( function(a,b) a:(b-1)
                , a = ageStart
                , b = ageEnd
                , SIMPLIFY = FALSE
            )) %>%
            unnest(cols = age) %>%
            left_join(dPop, by = c('age','year')) %>%
            with(sum(popsize))
        )
}) %>%
Reduce(f=rbind) %>%
group_by(year) %>%
mutate(densePop = (popsize/(ageEnd-ageStart))/sum(popsize))


plotAgeDist = function(x){
    ggplot(x)+
    geom_rect(aes(
        ymin = ageStart
        ,ymax = ageEnd
        ,xmin = 0
        ,xmax = densePop / max(ageDist$densePop)
    ), fill="#dddddd", size = 0.3)+

    geom_rect(aes(
        ymin = ageStart
        ,ymax = ageEnd
        ,xmin = 0
        ,xmax = denseCase / max(ageDist$denseCase)
    ), fill='red', size = 0.3)+
    facet_wrap(~year, nrow=1)+
    scale_y_continuous('Age', limits=c(0,65), breaks=c(0,25,50), expand=c(0,0))+
    scale_x_continuous('Density (DHF cases)'
        , expand=c(0,.1)
        , sec.axis = sec_axis( ~., name = "Density (population)")
    )+
    theme(
        panel.grid.major = element_blank()
        ,panel.grid.minor = element_blank()
        ,panel.border = element_rect(color="#bbbbbb")
        ,strip.background = element_rect(fill="white", color="#bbbbbb")
        ,axis.text.y = element_text(size=rel(.8), color="#bbbbbb")
        ,axis.ticks.y = element_line(color="#bbbbbb")
        ,axis.text.x = element_blank()
        ,axis.ticks.x = element_blank()
        ,axis.title.x.bottom = element_text(color = "red")
    )
}
    

# Figure 1
g = ggarrange(  gTs
    , ageDist %>%
        filter(year %in%  seq(1981,2017,by=9)) %>%
        plotAgeDist()
    ,nrow = 2
    ,ncol = 1
    ,heights = c(2,1)
    ,labels = "auto"
)
ggsave(g
    ,filename = paste0(plotdir,"/meanAge_dist.pdf")
    ,width = 3
    ,height = 4
)

# Figure S2
g = ageDist %>%
    ggplot()+
    geom_rect(aes(
        xmin = ageStart
        ,xmax = ageEnd
        ,ymin = 0
        ,ymax = densePop / max(ageDist$densePop)
    ), fill="#dddddd")+

    geom_rect(aes(
        xmin = ageStart
        ,xmax = ageEnd
        ,ymin = 0
        ,ymax = denseCase / max(ageDist$denseCase)
    ), fill='red')+
    facet_wrap(~year, nrow=4)+
    scale_x_continuous('Age', limits=c(0,65), breaks=c(0,25,50), expand=c(0,0))+
    scale_y_continuous('Density of population'
        , expand=c(0,.1)
        , sec.axis = sec_axis( ~., name = "Density of DHF cases")
    )+
    theme(
        panel.grid.major = element_blank()
        ,panel.grid.minor = element_blank()
        ,panel.border = element_rect(color="#bbbbbb")
        ,strip.background = element_rect(fill="white", color="#bbbbbb")
        ,axis.text.x = element_text(size=rel(.8), color="#bbbbbb")
        ,axis.ticks.x = element_line(color="#bbbbbb")
        ,axis.text.y = element_blank()
        ,axis.ticks.y = element_blank()
        ,axis.title.y.right = element_text(color = "red")
    )
ggsave(g
    ,filename = paste0(plotdir,"/meanAge_dist_all.pdf")
    ,width = 6
    ,height = 4
)
