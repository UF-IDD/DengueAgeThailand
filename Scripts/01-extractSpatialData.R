#!/usr/bin/env Rscript --vanilla
# rm(list=ls())


spatialdir = "01-extractSpatialData"

if(FALSE){

    indir = "00-ProvinceShapefiles"

    library(rgdal)
    library(rgeos)
    library(dplyr)
    library(ggplot2)

    # import mapping
    source('Scripts/importMapping.R')

    # import shapefile of Thai 77 provinces
    dataProjected <- readOGR(
        dsn = indir
        ,layer = gsub('\\.shp$','',list.files(indir, pattern='\\.shp$'))
    )
    # add to data a new column termed "id" composed of the rownames of data
    dataProjected@data$id <- rownames(dataProjected@data)

    # create dataframe of province centroids
    provinceCentroids <- dataProjected@data %>%
        cbind(
            # polygon centroids
            gCentroid(dataProjected,byid=TRUE) %>%
            as.data.frame %>%
            rename(cenX = x, cenY = y) %>%
            mutate(area = sapply(dataProjected@polygons, function(x){
                slot(x,'area')
            }))
        ) %>%
        left_join(mapping, by=c("PROV_NAMT"="province")) %>%
        select(id, region, province72, province77, cenX, cenY, area) %>%
        arrange(cenY)

    # simplify the polygons to save plotting time
    dataProjected <- gSimplify(dataProjected, 0.05, topologyPreserve=TRUE)
    
    # create polygon dataframe for ggplot
    provincePolygons <- fortify(dataProjected, region = "id")
  
    # save data to ease importation
    dir.create(spatialdir, recursive=T)
    save(
        provinceCentroids, provincePolygons
        ,file = paste0(spatialdir,'/provinceSpatialData.Rdata')
    )
} else {
    load(paste0(spatialdir,'/provinceSpatialData.Rdata'))
}

rm(spatialdir)




