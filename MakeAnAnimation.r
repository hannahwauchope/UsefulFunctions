FP <- "/Users/hannahwauchope/Dropbox/Work/Projects/PaleoSDMs/"
DataFP <- "/Users/hannahwauchope/Dropbox/Work/Data/"

source("~/Documents/GitHub/UsefulFunctions/VisualiseYourSpatialData.R")

#Get a polygon
#Download the shapefile, save it somewhere. 
WorldPolyFP <- paste0(DataFP, "GISLayers/CountryContinentWorldSHPS/World/TM_WORLD_BORDERS-0.3-SSudan.shp")

#### EXAMPLE RUN WITH MY DATA ####
# Create animations of SDM output, with fossil locations overlayed for the past, and GBIF locations for the present
#Running this for many species
PolSDMs3 <- read.csv(paste0(FP, "SDMs/Output/SDMSampleOutput.csv"))
load(file=paste0(FP, "GBIFDownloads/GBIFReduce.RData"))

SpecToModel <- c("Populus.tremuloides", "Betula.pendula", "Juniperus.communis", "Pinus.contorta", "Alnus.rubra", "Centaurea.scabiosa")
SpecToModel <- gsub("[.]", " ", SpecToModel)
SpecToModel <- data.frame(Species=SpecToModel, Region=c("North America", "Europe Russia Asia Middle East", "Europe Russia North America", "North America", "North America", "Europe Russia Asia Middle East"))

pblapply(1:nrow(SpecToModel), function(STM){
  Spec <- SpecToModel[STM,]$Species
  Region <- SpecToModel[STM,]$Region
  print(Spec)
  #Get the region for the polygon (we're just using one static one)
  if(Region=="Europe Russia Asia Middle East"){
    AnimPolys <- GetThatBasemap(WorldPolyFP, Region, c(-30, 180, -11, 82))
  } else {
    AnimPolys <- GetThatBasemap(WorldPolyFP, Region)
  }
  
  #Get the points for the relevant species, do a bit of cleaning
  AnimPoints <- subset(PolSDMs3, Species==Spec & AgeScale==1)
  AnimPoints <- vect(AnimPoints, geom=c("x", "y"), crs=MollCRS)
  AnimPoints <- project(AnimPoints, WGSCRS)
  AnimPoints <- crop(AnimPoints, AnimPolys)
  AnimPoints <- AnimPoints[!is.na(AnimPoints$SDMScale),]
  names(AnimPoints)[names(AnimPoints)=="Clamp"] <- "NumAxesOutsideRange"
  AnimPoints$NumAxesOutsideRange <- as.factor(AnimPoints$NumAxesOutsideRange)
  AnimPoints$PointVal <- factor(round_any(AnimPoints$SDMScale, 0.25), levels=c(0,0.25, 0.5, 0.75, 1))
  
  #Get the time steps
  TimeSteps <- c(rev(sort(as.numeric(unique(AnimPoints[AnimPoints$Dataset_Sample!="GBIF"]$Ages)))), unique(AnimPoints[AnimPoints$Dataset_Sample=="GBIF"]$Ages))
  TimeSteps <- as.character(TimeSteps)
  
  #Split points vector into a list of vectors based on timesteps. 
  AnimPoints <- pblapply(TimeSteps, function(TS) AnimPoints[AnimPoints$Ages == TS,])
  names(AnimPoints) <- TimeSteps
  
  print("Compile Rasters")

  AnimRasts <- pblapply(TimeSteps, function(TS){
    rasty <- rast(paste0(FP, "/SDMs/", gsub(" ", ".", Spec), "/HindcastRasters/Hindcast_", TS, ".tif")) #Read in raster
    rasty <- project(rasty, WGSCRS) #Project
    if(!is.null(AnimPolys)){ #Crop to same region as the polygon
      rasty <- crop(rasty, AnimPolys)
      rasty <- mask(rasty, AnimPolys)
    }
    rasty[!is.na(rasty)] <- round_any(rasty[!is.na(rasty)], 100)/1000 #Convert raster values into 10 breaks (0,1,0.1)
    return(rasty)
  }) #, mc.cores=ncores
  names(AnimRasts) <- TimeSteps
  TimeSteps <- as.character(TimeSteps)
  
  print("Create Images")
  AnimPointsPast <- AnimPoints[!names(AnimPoints) %in% unique(PolSDMs3[PolSDMs3$Dataset_Sample=="GBIF",]$Ages)]
  AnimPointsPresent<- AnimPoints[names(AnimPoints) %in% unique(PolSDMs3[PolSDMs3$Dataset_Sample=="GBIF",]$Ages)]
  AnimRastsPast <- AnimRasts[!names(AnimRasts) %in% unique(PolSDMs3[PolSDMs3$Dataset_Sample=="GBIF",]$Ages)]
  AnimRastsPresent<- AnimRasts[names(AnimRasts) %in% unique(PolSDMs3[PolSDMs3$Dataset_Sample=="GBIF",]$Ages)]
  TimeStepsPast <- TimeSteps[!TimeSteps %in% unique(PolSDMs3[PolSDMs3$Dataset_Sample=="GBIF",]$Ages)]
  TimeStepsPresent <- TimeSteps[TimeSteps %in% unique(PolSDMs3[PolSDMs3$Dataset_Sample=="GBIF",]$Ages)]
  
  names(AnimPointsPast) <- as.character(as.numeric(names(AnimPointsPast))/1000)
  names(AnimRastsPast) <- as.character(as.numeric(names(AnimRastsPast))/1000)
  TimeStepsPast <- as.character(as.numeric(TimeStepsPast)/1000)
  
  CreateThoseImages(AnimPointsPast, AnimPolys, AnimRastsPast, TimeStepsPast, MapOrder=c("Poly", "Rast", "Points"), AnimFP = paste0(FP, "/Animations/", Spec, "/"),
                        PointColVal = "PointVal", 
                        PointColName = "Probability of Occurrence",
                        PointCols = brewer.pal(5, "YlOrRd"))
  CreateThoseImages(AnimPointsPresent, AnimPolys, AnimRastsPresent, TimeStepsPresent, MapOrder=c("Poly", "Rast", "Points"), AnimFP = paste0(FP, "/Animations/", Spec, "/"),
                        PointColVal = "PointVal", 
                        PointColName = "Probability of Occurrence",
                        PointCols = brewer.pal(5, "YlOrRd"),
                        PointSize = 0.5,
                        overwrite=FALSE, 
                        TimeStepTitle=TRUE,
                        SaveReps=5)#
  MakeThatAnimation(AnimFP = paste0(FP, "/Animations/", Spec, "/"), FramesPerSecond=5, SaveName=Spec)
}) #, mc.cores=6


#### GRAVEYARD ####
#This function is a fairly janky function to get background polygon.
#WorldPolyFP needs to be the FP to where you saved the WorldPoly
#Region can either be:
### a continental region (out of World, Europe, Russia, Africa, Amercicas, North America, South America, Asia Middle East, Antarctica)
### A vector of multiple continental regions (E.g. c("Europe", "Russia", "North America"))
### A vector of country names (can be found by running "ReturnCountryNames=TRUE"), that you want. IF YOU CHOOSE THIS, name CountryNames=TRUE
#Crop coords: a vector of decimal degree coordinates for if you want to crop the polygon. Give as: c(xmin, xmax, ymin, ymax)

GetBackground <- function(WorldPolyFP, Region, CropCoords=NULL, CountryNames=FALSE, ReturnCountryNames=FALSE){
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

#This functions creates the images that will then be stacked together to make an animation 
#All this code is structured around maps of animations through time, so, you know, bear that in mind. 
#I've set it up so you can provide any combination of location points, polygons, or rasters. So you could have just rasters, rasters over polygons, points of polyons, points over rasters over polygons etc etc. 
#I'm going to refer to the Points/Polygons/Rasters as PPR as shorthand
#PPR can either be provide as one vector/raster, or as a list of them. If just one, it'll be static throughout the animation, if many then they'll change
#For any PPRs provided as a list, the list must be the same length for each one AND the same length as TimeSteps
#TimeSteps is a character vector of the timesteps you're using (e.g. c(1940,1950,1960, 1970))
#MapOrder specifies the order you want PPR in (from base layer to top layer). 
#AnimFP is a path to a folder where you want the animation to be created. 
#TimeScale = a logical for whether you want a time bar at the bottom that labels where in time you are

#Now we get in to a bunch of aesthetic choices you may wanna make. There are many more that could be made, I may edit this more at some point
#PointColVal is the name of a column in your dataframe that you want to use to colour your points
#Point Col Name is want you want PointColVal to be called in the map legend
#PointCols is a vector of colours if you want to manually specify point colour
#PointSize = size of your points
#RastViridis = a logical, if TRUE the raster colour scheme will be viridis. 
#RastColName = name of the raster legend
#Overwrite = logical for whether you want to erase old images from a previous run
#TimeStepTitle = logical for whether you want the name of the time step to be in the top right corner
#Save reps = getting sneaky! If you want certain map images to hang around for longer (e.g. the last image), run this function separately for that/those timestep(s) and up the number of save reps.
##Then this image will make up multiple frames in the final animation

CreateAnimationImages <- function(AnimPoints=NULL, AnimPolys=NULL, AnimRasts=NULL, TimeSteps, MapOrder=c("Rast", "Poly", "Points"),
                                  AnimFP = paste0(FP, "/Animations/", Spec, "/"),
                                  PointColVal = "NumAxesOutsideRange", 
                                  PointColName = "Number of PCA's outside Range",
                                  PointCols = brewer.pal(8, "RdYlBu")[3:1],
                                  PointSize=1.5,
                                  RastViridis = TRUE,
                                  RastColName = "Probability of Occurrence",
                                  overwrite=TRUE, TimeScale=TRUE,
                                  TimeStepTitle=FALSE,
                                  SaveReps=1){
  dir.create(paste0(AnimFP, "Images/"), showWarnings=FALSE)
  if(overwrite==TRUE){
    sapply(list.files(paste0(AnimFP, "Images/"), full.names=TRUE), unlink)
  }
  
  if(!is.null(AnimRasts)){
    RastsMinMax <- lapply(AnimRasts, minmax)
    RastMin <- min(sapply(RastsMinMax, function(x) x[[1]]))
    RastMax <- max(sapply(RastsMinMax, function(x) x[[2]])) 
  }
  
  if(!is.null(AnimPoints)){
    if(class(AnimPoints)=="list"){
      AnimPoints <- lapply(AnimPoints, function(x){
        names(x)[names(x) == PointColVal] <- "PointCol"
        return(x)
      })
    } else {
      names(AnimPoints)[names(AnimPoints) == PointColVal] <- "PointCol"
    }
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
    
    AnimImage <- ggplot()+
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
      AnimImage <- AnimImage + geom_spatraster(data=AnimRast)
    } else if (MapOrder[[1]]=="Poly"){
      AnimImage <- AnimImage + geom_spatvector(data=AnimPoly, colour=NA, fill="gray70", alpha=0.5)
    } else if (MapOrder[[1]] == "Points"){
      AnimImage <- AnimImage + geom_spatvector(data=AnimPoint, aes(colour=PointCol), size=PointSize)
    }
    
    if(MapOrder[[2]]=="Rast"){
      AnimImage <- AnimImage + geom_spatraster(data=AnimRast)
    } else if (MapOrder[[2]]=="Poly"){
      AnimImage <- AnimImage + geom_spatvector(data=AnimPoly, colour=NA, fill="gray70", alpha=0.5)
    } else if (MapOrder[[2]] == "Points"){
      AnimImage <- AnimImage + geom_spatvector(data=AnimPoint, aes(colour=PointCol), size=PointSize)
    }
    
    if(MapOrder[[3]]=="Rast"){
      AnimImage <- AnimImage + geom_spatraster(data=AnimRast)
    } else if (MapOrder[[3]]=="Poly"){
      AnimImage <- AnimImage + geom_spatvector(data=AnimPoly, colour=NA, fill="gray70", alpha=0.5)
    } else if (MapOrder[[3]] == "Points"){
      AnimImage <- AnimImage + geom_spatvector(data=AnimPoint, aes(colour=PointCol), size=PointSize)
    }
    
    if("Points" %in% MapOrder & !is.null(PointCols)){
      if(PointCols[[1]]=="viridis"){
        AnimImage <- AnimImage + scale_colour_viridis(na.value=NA, direction=1, name=PointColName, discrete=TRUE, option="A")
      } else {
        AnimImage <- AnimImage + scale_colour_manual(values=PointCols, name=PointColName, drop=FALSE)
      }
    }
    
    if("Rast" %in% MapOrder & RastViridis==TRUE){
      AnimImage <- AnimImage + scale_fill_viridis(na.value=NA, direction=1, name=RastColName, limits=c(RastMin, RastMax))
    }
    
    AnimImage <- AnimImage + guides(fill = guide_legend(order = 1), 
                                    colour = guide_legend(order = 2))
    
    if(TimeStepTitle){
      AnimImage <- AnimImage+labs(title=TS)
    }
    
    if(TimeScale == TRUE){
      ScalePlot <- ggplot()+
        geom_point(aes(x=as.numeric(TS)/1000, y=0.01), size=4, colour=ifelse(PointCols[[1]]=="viridis", "#5ec962", PointCols[[5]]))+
        geom_line(aes(x=c(0, 22), y=0))+
        scale_x_reverse(expand=c(0.1,0.1), breaks=rev(seq(0, 22, 2)))+
        scale_y_continuous(expand=c(0,0), limits=c(0,1))+
        xlab("Thousand years before present")+
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

#This function creates the animation
#AnimFP = same as the filepath in CreateAnimationImages
#Frames per second = how fast you want the animation to be. More = more frames per second!
#Save name is what you want the animation gif to be saved as 

CreateAnimation <- function(AnimFP, FramesPerSecond=2, SaveName){
  print("compile animation")
  timestamps <- data.frame(filenames=list.files(paste0(AnimFP, "Images"), full.names=TRUE), creationdate=file.info(list.files(paste0(AnimFP, "Images"), full.names=TRUE))$mtime)
  timestamps <- timestamps[order(timestamps$creationdate),]
  timestamps <- sapply(timestamps$filenames, function (x) image_read(x))
  myanimation <- image_animate(image_join(timestamps), fps=FramesPerSecond)
  print("save animation")
  image_write(myanimation, paste0(AnimFP, SaveName,".gif"))
}
