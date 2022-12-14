### Visualise Your Spatial Data ###
#Written by Hannah Wauchope, last edited Tuesday 13th December, 2022

#### Overview ####
#This script contains variation on a theme of occasions where you want to view maps of things: either points, polygons, rasters, or any combination of the three. 
#I refer to those three as PPR in this script for simplicity

#First, there is a function to get (fairly basic) basemap polygons, you can choose either the whole world, continental regions, or specific countries
#This is "GetThatBasemap"
#You can use this as the polygon in the next functions

#Then, there is a script to map any of your PPR. You can specify which of PPR you are mapping, and make a bunch of aesthetic decisions)
#This is "MakeThatMap"

#Next, there are two functions to build animations. 

#The first, "CreateThoseImages" uses MakeThatMap to make the map image "frames" that will then be stacked together into an animated gif. It'll also add a timescale bar if you want.

#Finally, "MakeThatAnimation" stitches the images together and saves an animation. 

#### Initialise ####
library(tidyverse)
library(tidyterra)
library(terra)
library(viridis)
library(gtools)
library(magick)
library(pbapply)
library(pbmcapply)
library(RColorBrewer)
library(cowplot)
library(plyr)

WGSCRS <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
MollCRS <- "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"

#### GetThatBasemap ####

##BEFORE YOU RUN THIS FUNCTION##
#You need to download the following shapefile and save it somewhere:
# TM_WORLD_BORDERS-0.3.zip from
# https://thematicmapping.org/downloads/world_borders.php

#This function is a fairly janky function to get background polygon.

#ARGUMENTS:
#"WorldPolyFP" needs to be the FP to where you saved the TMWorldBorders shapefile
#"Region" can either be:
### a continental region (out of World, Europe, Russia, Africa, Amercicas, North America, South America, Asia Middle East, Antarctica)
### A vector of multiple continental regions (E.g. c("Europe", "Russia", "North America"))
### A vector of country names (can be found by running "ReturnCountryNames=TRUE"), that you want. 
### IF YOU CHOOSE THIS, make "CountryNames=TRUE"
#"CropCoords": a vector of decimal degree coordinates for if you want to crop the polygon. Give as: c(xmin, xmax, ymin, ymax)

GetThatBasemap <- function(WorldPolyFP, Region, CropCoords=NULL, CountryNames=FALSE, ReturnCountryNames=FALSE){
  WorldOutline <- vect(WorldPolyFP)
  if(ReturnCountryNames==TRUE){
    return(unique(WorldOutline$NAME))
  }
  WorldOutline <- project(WorldOutline, WGSCRS)
  RegionPoly <- NULL
  if(grepl("World", Region)){
    RegionPoly <- WorldOutline
  }
  if(grepl("Europe", Region)){
    RegionPoly <- WorldOutline[WorldOutline$REGION==150 & WorldOutline$NAME!="Russia",]
  }
  if(grepl("Russia", Region)){
    if(is.null(RegionPoly)){
      RegionPoly <- WorldOutline[WorldOutline$NAME=="Russia",]
    } else {
      RegionPoly <- union(RegionPoly, WorldOutline[WorldOutline$NAME=="Russia",])
    }
  }
  if(grepl("Africa", Region)){
    if(is.null(RegionPoly)){
      RegionPoly <- WorldOutline[WorldOutline$REGION==2,]
    } else {
      RegionPoly <- union(RegionPoly, WorldOutline[WorldOutline$REGION==2,])
    }
  }
  if(grepl("America", Region)){
    Americas <- WorldOutline[WorldOutline$REGION==19,]
    Americas <- crop(Americas, ext(-179, -10, -60, 90))
    if(grepl("Americas", Region)){
      if(is.null(RegionPoly)){
        RegionPoly <- Americas
      } else {
        RegionPoly <- union(RegionPoly, Americas)
      }
    }
    SouthAmerica <- crop(Americas, ext(-179, -10, -60, 12))
    SouthAmerica <- SouthAmerica[!SouthAmerica$NAME %in% c("Nicaragua", "Costa Rica", "Panama")]
    if(grepl("South America", Region)){
      if(is.null(RegionPoly)){
        RegionPoly <- SouthAmerica
      } else {
        RegionPoly <- union(RegionPoly, SouthAmerica)
      }
    }
    if(grepl("North America", Region)){
      if(is.null(RegionPoly)){
        RegionPoly <- Americas[!Americas$NAME %in% SouthAmerica$NAME,]
      } else {
        RegionPoly <- union(RegionPoly, Americas[!Americas$NAME %in% SouthAmerica$NAME,])
      }
    }
  }
  if(grepl("Asia Middle East", Region)){
    if(is.null(RegionPoly)){
      RegionPoly <- WorldOutline[WorldOutline$REGION==142,]
    } else {
      RegionPoly <- union(RegionPoly, WorldOutline[WorldOutline$REGION==142,])
    }
  }
  if(grepl("Australasia", Region)){
    Australasia <- WorldOutline[WorldOutline$REGION==9,]
    if(is.null(RegionPoly)){
      RegionPoly <- project(Australaisa, "+proj=eqearth +lon_0=150 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +type=crs")
    } else {
      RegionPoly <- union(RegionPoly, Australasia)
    }
  }
  if(grepl("Antarctica", Region)){
    Antarctica <- WorldOutline[WorldOutline$REGION==0,]
    Antarctica <- crop(Antarctica, ext(-180, 180, -90, -60))
    if(is.null(RegionPoly)){
      RegionPoly <- project(Antarctica, "+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +type=crs")
    } else {
      RegionPoly <- union(RegionPoly, Antarctica)
    }
  }
  if(CountryNames==TRUE){
    RegionPoly <- WorldOutline[WorldOutline$NAME %in% Region,]
  }
  if(is.null(RegionPoly)){
    stop("Given Region(s) not recognised")
  }
  if(!is.null(CropCoords)){
    RegionPoly <- crop(RegionPoly, ext(CropCoords))
  }
  return(RegionPoly)
}

#### MakeThatMap ####
#This function makes a map of points, polygons, rasters, or a combination ('PPR')
#All of them need to be in terra format (either SpatVect for points and polygons, or SpatRast for Rasters), and can be any projection you like as long as they're consistent 
#Alternatively, you can give points as a dataframe (with coordinates named as either x,y; lon,lat; longitude, latitude; or Longitude, Latitude), but these must be in decimal degrees, not any other projection,as must be the raster or polygon if you include it

#ARGUMENTS
#MapPoints = Your point data
#MapPolys = Your polygon data
#MapRaster = Your raster data
#MapOrder = The order you want them in, as a character vector including "Rast", "Poly", "Points" (defaults to raster on base, then polygon, then points on top)
## If you've only provided, e.g., points and raster, then this should be c("Rast", "Points")

#Now we get in to a bunch of aesthetic choices you may wanna make. There are many more that could be made, I may edit this more at some point
#PointColVal is the name of a column in your dataframe that you want to use to colour your points. AT THE MOMENT THIS MUST BE A CHARACTER, NOT NUMERIC!!
#Point Col Name is want you want PointColVal to be called in the map legend
#PointCols is a vector of colours if you want to manually specify point colour (OR "viridis" to use the viridis colour palette)
#PointSize = size of your points

#RastColName = name of the raster legend
#RastMin = Min val in your rast colourscale (if null R will calculate this for you)
#RastMax = Max val in your rast colourscale (if null R will calculate this for you)

#PolyFillCol = Fill colour of the polygon
#PolyLineCol = Line colour of the polygon
#PolyAlpha = alpha of the polygon


MakeThatMap <- function(MapPoint=NULL, MapPoly=NULL, MapRast=NULL, MapOrder=c("Rast", "Poly", "Points"),
                        PointColVal = "NumAxesOutsideRange", 
                        PointColName = "Number of PCA's outside Range",
                        PointCols = brewer.pal(8, "RdYlBu")[4:1],
                        PointSize=1.5,
                        RastColName = "Probability of Occurrence",
                        RastMin = NULL,
                        RastMax = NULL,
                        PolyFillCol = "gray70",
                        PolyLineCol = NA,
                        PolyAlpha = 0.5){
  if(class(MapPoint)!="SpatVector"){
    GeomLon <- names(MapPoint)[names(MapPoint) %in% c("x", "lon", "longitude", "Longitude", "Lon")]
    GeomLat <- names(MapPoint)[names(MapPoint) %in% c("y", "lat", "latitude", "Latitude", "Lat")]
    MapPoint <- vect(MapPoint, geom=c(GeomLon, GeomLat), crs=WGSCRS)
  }
  names(MapPoint)[names(MapPoint) == PointColVal] <- "PointCol"
  if(class(MapPoly)!="SpatVector"){
    MapPoly <- vect(MapPoly)
  }
  
  if(class(MapRast)!="SpatRaster"){
    MapRast <- rast(MapRast)
  }
  
  while(length(MapOrder) != 3){
    MapOrder <- c(MapOrder, "NULL")
  }
  
  MapImage <- ggplot()+
    theme(#aspect.ratio = ifelse(SDMS==TRUE, 0.38, 0.8),
      panel.grid = element_blank(), 
      axis.ticks = element_blank(),
      axis.text = element_blank(),
      axis.title = element_blank(), 
      plot.background=element_blank(),
      panel.background = element_blank(),
      legend.key=element_blank(),
      legend.text = element_text(size=6),
      legend.justification =  "top",
      legend.spacing.x = unit(0.1, 'cm'),
      plot.title = element_text(size=10))
  
  if(MapOrder[[1]]=="Rast"){
    MapImage <- MapImage + geom_spatraster(data=MapRast)
  } else if (MapOrder[[1]]=="Poly"){
    MapImage <- MapImage + geom_spatvector(data=MapPoly, colour=PolyLineCol, fill=PolyFillCol, alpha=PolyAlpha)
  } else if (MapOrder[[1]] == "Points"){
    MapImage <- MapImage + geom_spatvector(data=MapPoint, aes(colour=PointCol), size=PointSize)
  }
  
  if(MapOrder[[2]]=="Rast"){
    MapImage <- MapImage + geom_spatraster(data=MapRast)
  } else if (MapOrder[[2]]=="Poly"){
    MapImage <- MapImage + geom_spatvector(data=MapPoly, colour=PolyLineCol, fill=PolyFillCol, alpha=PolyAlpha)
  } else if (MapOrder[[2]] == "Points"){
    MapImage <- MapImage + geom_spatvector(data=MapPoint, aes(colour=PointCol), size=PointSize)
  }
  
  if(MapOrder[[3]]=="Rast"){
    MapImage <- MapImage + geom_spatraster(data=MapRast)
  } else if (MapOrder[[3]]=="Poly"){
    MapImage <- MapImage + geom_spatvector(data=MapPoly, colour=PolyLineCol, fill=PolyFillCol, alpha=PolyAlpha)
  } else if (MapOrder[[3]] == "Points"){
    MapImage <- MapImage + geom_spatvector(data=MapPoint, aes(colour=PointCol), size=PointSize)
  }
  
  if("Points" %in% MapOrder & !is.null(PointCols)){
    if(PointCols[[1]]=="viridis"){
      MapImage <- MapImage + scale_colour_viridis(na.value=NA, direction=1, name=PointColName, discrete=TRUE, option="A")
    } else {
      MapImage <- MapImage + scale_colour_manual(values=PointCols, name=PointColName, drop=FALSE)
    }
  }
  
  if(is.null(RastMin)){RastMin <- minmax(MapRast)[[1]]}
  if(is.null(RastMax)){RastMax <- minmax(MapRast)[[2]]}
  
  if("Rast" %in% MapOrder){
    MapImage <- MapImage + scale_fill_viridis(na.value=NA, direction=1, name=RastColName, limits=c(RastMin, RastMax))
  }
  
  MapImage <- MapImage + guides(fill = guide_legend(order = 1), 
                                  colour = guide_legend(order = 2))
  
  return(MapImage)
}


#### CreateThoseImages ####
#This functions creates the images that will then be stacked together to make an animation
#All this code is structured around maps of animations through time, so, you know, bear that in mind.
#You can either just have the maps as an animation, or you can add a timeline bar at the moment with a dot that shows where in time you are ("Time Scale")
#I've set it up so you can provide any combination of location points, polygons, or rasters. So you could have just rasters, rasters over polygons, points of polyons, points over rasters over polygons etc etc. 
#I'm going to refer to the Points/Polygons/Rasters as PPR as shorthand
#PPR can either be provide as one vector/raster, or as a list of them. If just one, it'll be static throughout the animation, if many then they'll change
#For any PPRs provided as a list, the list must be the same length for each one AND the same length as TimeSteps

### ARGUMENTS
#AnimPoints = An individual points vector or list of vectors
#AnimPolys = An individual polygon or list of polygons
#AnimRast = An individual raster or list of rasters
#TimeSteps is a character vector of the timesteps you're using (e.g. c(1940,1950,1960, 1970))
#MapOrder specifies the order you want PPR in (from base layer to top layer). 
#AnimFP is a path to a folder where you want the animation to be created. 
#TimeScale = a logical for whether you want a time bar at the bottom that labels where in time you are
#Overwrite = logical for whether you want to erase old images from a previous run
#TimeStepTitle = logical for whether you want the name of the time step to be in the top right corner
#TimeScaleLab = The x axis label for the time scale
#SaveReps = getting sneaky! If you want certain map images to hang around for longer (e.g. the last image), run this function separately for that/those timestep(s) and up the number of save reps.
##Then this image will make up multiple frames in the final animation

#### Finally, there are aesthetics (these are the SAME as in "MakeThatMap", and are passed straight onto it)
#PointColVal is the name of a column in your dataframe that you want to use to colour your points. AT THE MOMENT THIS MUST BE A CHARACTER, NOT NUMERIC!!
#Point Col Name is want you want PointColVal to be called in the map legend
#PointCols is a vector of colours if you want to manually specify point colour (OR "viridis" to use the viridis colour palette)
#PointSize = size of your points

#RastColName = name of the raster legend

#PolyFillCol = Fill colour of the polygon
#PolyLineCol = Line colour of the polygon
#PolyAlpha = alpha of the polygon

CreateThoseImages <- function(AnimPoints=NULL, AnimPolys=NULL, AnimRasts=NULL, TimeSteps, MapOrder=c("Rast", "Poly", "Points"),
                                  AnimFP = paste0(FP, "/Animations/", Spec, "/"),
                                  overwrite=TRUE, TimeScale=TRUE,
                                  TimeStepTitle=FALSE, SaveReps=1,
                                  TimeScaleLab = "Thousand Years Before Present",
                                  PointColVal = "NumAxesOutsideRange", 
                                  PointColName = "Number of PCA's outside Range",
                                  PointCols = brewer.pal(8, "RdYlBu")[4:1],
                                  PointSize=1.5,
                                  RastColName = "Probability of Occurrence",
                                  PolyFillCol = "gray70",
                                  PolyLineCol = NA,
                                  PolyAlpha = 0.5){
  dir.create(paste0(AnimFP, "Images/"), showWarnings=FALSE)
  if(overwrite==TRUE){
    sapply(list.files(paste0(AnimFP, "Images/"), full.names=TRUE), unlink)
  }
  
  if(!is.null(AnimRasts)){
    RastsMinMax <- lapply(AnimRasts, minmax)
    RastMin <- min(sapply(RastsMinMax, function(x) x[[1]]))
    RastMax <- max(sapply(RastsMinMax, function(x) x[[2]])) 
  }

  while(length(MapOrder) != 3){
    MapOrder <- c(MapOrder, "NULL")
  }
  
  pblapply(TimeSteps, function(TS){
    if(class(AnimPoints)=="list"){
      AnimPoint <- AnimPoints[[TS]]
    } else {AnimPoint <- AnimPoints}
    if(class(AnimRasts)=="list"){
      AnimRast <- AnimRasts[[TS]]
    } else {AnimRast <- AnimRasts}
    if(class(AnimPolys)=="list"){
      AnimPoly <- AnimPolys[[TS]]
    } else {AnimPoly <- AnimPolys}
    
    AnimImage <- MakeThatMap(AnimPoint, AnimPoly, AnimRast, 
                             PointColVal = PointColVal, 
                             PointColName = PointColName,
                             PointCols = PointCols,
                             PointSize=PointSize,
                             RastColName = RastColName,
                             RastMin = RastMin,
                             RastMax = RastMax,
                             PolyFillCol = PolyFillCol,
                             PolyLineCol = PolyLineCol,
                             PolyAlpha = PolyAlpha)
    if(TimeStepTitle){
      AnimImage <- AnimImage+labs(title=TS)
    }
    
    if(TimeScale == TRUE){
      ScalePlot <- ggplot()+
        geom_point(aes(x=as.numeric(TS), y=0.01), size=4, colour=ifelse(PointCols[[1]]=="viridis", "#5ec962", PointCols[[5]]))+
        geom_line(aes(x=c(as.numeric(tail(TimeSteps,1)), as.numeric(TimeSteps[[1]])), y=0))+
        scale_x_reverse(expand=c(0.1,0.1), breaks=rev(seq(floor(as.numeric(tail(TimeSteps,1))), ceiling(as.numeric(TimeSteps[[1]])), round((as.numeric(TimeSteps[[1]])-as.numeric(tail(TimeSteps,1)))/10))))+
        scale_y_continuous(expand=c(0,0), limits=c(0,1))+
        xlab(TimeScaleLab)+
        theme(axis.line.y=element_blank(), axis.text.y = element_blank(), axis.ticks.y=element_blank(), axis.title.y = element_blank(),
              plot.background = element_blank(), panel.background = element_blank(), panel.grid = element_blank(),
              axis.text.x = element_text(size=8, colour="black"), axis.title.x = element_text(size=10, colour="black"))
      
      AnimImage <- plot_grid(AnimImage, ScalePlot, ncol=1, rel_heights = c(0.9, 0.2)) +
        theme(plot.background = element_rect(fill = "white", colour = NA))
    }
    for(i in c(1:SaveReps)){
      ggsave(paste0(AnimFP, "Images/", TS, "_", i, ".png"), AnimImage, height=150, width=250, units="mm", dpi=500)
    }
  })
  return("Images Created")
}

#### MakeThatAnimation ####
#Finally, function creates the animation
#NOTE. This reads in all the files from the "Images" folder (within the AnimFP you've specified) ***in order they were written***
#This means if you wanted to make an animation with various data sources, run "CreateThoseImages" with your data sources in the order you want them to appear in the animation (just make sure to set overwrite to false so it doesn't erase each time!)
#You can play around with the number of SaveReps you use in "CreateThoseImages" to make certain maps last for more or less time in the animation

##ARGUMENTS 
#AnimFP = same as the filepath in CreateAnimationImages
#Frames per second = how fast you want the animation to be. More = more frames per second!
#Save name is what you want the animation gif to be saved as 

MakeThatAnimation <- function(AnimFP, FramesPerSecond=2, SaveName){
  print("compile animation")
  timestamps <- data.frame(filenames=list.files(paste0(AnimFP, "Images"), full.names=TRUE), creationdate=file.info(list.files(paste0(AnimFP, "Images"), full.names=TRUE))$mtime)
  timestamps <- timestamps[order(timestamps$creationdate),]
  timestamps <- sapply(timestamps$filenames, function (x) image_read(x))
  myanimation <- image_animate(image_join(timestamps), fps=FramesPerSecond)
  print("save animation")
  image_write(myanimation, paste0(AnimFP, SaveName,".gif"))
}
