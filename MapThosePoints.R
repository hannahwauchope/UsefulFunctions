#A simple function to plot your points on a world map!
#If doing this for a figure *highly* recommend getting your own better shapefile (if for the world ideally projecting it in Robinson or something)

library(ggplot2)
library(tidyterra)
SamplePoints <- data.frame(lon=c(130,49,-53,74,167,-102), lat=c(-30,-50,20,60,50,0))
SamplePoints$Color <- c(1,5,3,3,5,2)

#Spec points = dataframe with a row per point and lat longs
#coords = the column names for long (x) and lat (y)
#Colourscale = true if you want a column for colour
#ColourID = if colourscale is TRUE then this is the column name for the colour variable
MapPointsSimple <- function(SpecPoints, coords = c("x", "y"), colourscale=FALSE, colourid){
  WorldMap <- borders(database="world", colour="grey70", fill="grey70")
  SpecPoints <- as.data.frame(SpecPoints)
  SpecPoints$x <- SpecPoints[,c(coords[[1]])]
  SpecPoints$y <- SpecPoints[,c(coords[[2]])]
  
  if(colourscale==TRUE){
    SpecPoints$ColorScale <- SpecPoints[,c(colourid)]
  }
  PastPoints  <- ggplot(SpecPoints) + WorldMap+
    theme(panel.grid = element_blank(), 
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank(),
          plot.background=element_blank(),
          panel.background = element_rect(fill="white"),
          legend.key=element_blank(),
          legend.text = element_text(size=6),
          legend.justification =  "top",
          legend.spacing.x = unit(0.1, 'cm'),
          legend.position=NULL)
  if(colourscale==TRUE){
    PastPoints <- PastPoints +  geom_point(aes(x=x, y=y, colour=ColorScale), size=1, stroke=0.8)+
      scale_colour_viridis()
  } else {
    PastPoints <- PastPoints +  geom_point(aes(x=x, y=y), size=1, stroke=0.8)
  }
  return(PastPoints)
}
MapPointsSimple(SamplePoints, coords=c("lon", "lat"), colourscale=TRUE, colourid="Color")

MapPointsWithPolygon <- function(SpecPoints, BasePoly, coords = c("x", "y"), colourscale=FALSE, colourid, PolyFillCol="grey70", PolyLineCol="grey70", MapProj = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"){
  SpecPoints <- as.data.frame(SpecPoints)
  SpecPoints$x <- SpecPoints[,c(coords[[1]])]
  SpecPoints$y <- SpecPoints[,c(coords[[2]])]
  if(colourscale==TRUE){
    SpecPoints$ColorScale <- SpecPoints[,c(colourid)]
  }
  SpecPoints <- vect(SpecPoints, geom=c("x", "y"), crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
  if(class(BasePoly)!="SpatVector"){
    BasePoly <- vect(BasePoly)
  }
  SpecPoints <- project(SpecPoints, BasePoly)
  
  PastPoints  <- ggplot() + 
    geom_spatvector(data=BasePoly, colour=PolyLineCol, fill=PolyFillCol)+
    theme(panel.grid = element_blank(), 
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank(),
          plot.background=element_blank(),
          panel.background = element_rect(fill="white"),
          legend.key=element_blank(),
          legend.text = element_text(size=6),
          legend.justification =  "top",
          legend.spacing.x = unit(0.1, 'cm'),
          legend.position=NULL)
  if(colourscale==TRUE){
    PastPoints <- PastPoints + geom_spatvector(data=SpecPoints, aes(colour=ColorScale), size=1, stroke=0.8)+
      scale_colour_viridis()
  } else {
    PastPoints <- PastPoints + geom_spatvector(data=SpecPoints, size=1, stroke=0.8)
  }
  return(PastPoints)
}


ModernTemps <- rast(paste0(DataFP, "/Climate_CRU/BioClim/", "1901_1950", "BioClim.grd"))[[1]]
crs(ModernTemps) <- WGSCRS
ModernTemps <- project(ModernTemps, MollCRS)
BaseRast <- ModernTemps

MapPointsWithRaster(SamplePoints, BaseRast, coords=c("lon", "lat"), colourscale=TRUE, colourid="Color")

MapPointsWithRaster <- function(SpecPoints, BaseRast, coords = c("x", "y"), colourscale=FALSE, colourid,
                                PointLegend = "none", RastLegend = "none"){
  SpecPoints <- as.data.frame(SpecPoints)
  SpecPoints$x <- SpecPoints[,c(coords[[1]])]
  SpecPoints$y <- SpecPoints[,c(coords[[2]])]
  if(colourscale==TRUE){
    SpecPoints$ColorScale <- SpecPoints[,c(colourid)]
  }
  SpecPoints <- vect(SpecPoints, geom=c("x", "y"), crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
  if(class(BaseRast)!="SpatRaster"){
    BaseRast <- rast(BaseRast)
  }
  SpecPoints <- project(SpecPoints, BaseRast)
  
  PastPoints  <- ggplot() + 
    geom_spatraster(data=BaseRast)+
    scale_fill_viridis(na.value="white", option="magma")+
    theme(panel.grid = element_blank(), 
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank(),
          plot.background=element_blank(),
          panel.background = element_rect(fill="white"),
          legend.key=element_blank(),
          legend.text = element_text(size=6),
          legend.justification =  "top",
          legend.spacing.x = unit(0.1, 'cm'),
          legend.position=RastLegend)
  if(colourscale==TRUE){
    PastPoints <- PastPoints + geom_spatvector(data=SpecPoints, aes(colour=ColorScale), size=1, stroke=0.8)+
      scale_colour_viridis()+
      theme(legend.position=PointLegend)
  } else {
    PastPoints <- PastPoints + geom_spatvector(data=SpecPoints, size=1, stroke=0.8)
  }
  return(PastPoints)
}

MapPointsPolygonRaster <- function(SpecPoints, BaseRast, coords = c("x", "y"), colourscale=FALSE, colourid,
                                   PointLegend = "none", RastLegend = "none"){
  SpecPoints <- as.data.frame(SpecPoints)
  SpecPoints$x <- SpecPoints[,c(coords[[1]])]
  SpecPoints$y <- SpecPoints[,c(coords[[2]])]
  if(colourscale==TRUE){
    SpecPoints$ColorScale <- SpecPoints[,c(colourid)]
  }
  SpecPoints <- vect(SpecPoints, geom=c("x", "y"), crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
  if(class(BaseRast)!="SpatRaster"){
    BaseRast <- rast(BaseRast)
  }
  SpecPoints <- project(SpecPoints, BaseRast)
  
  PastPoints  <- ggplot() + 
    geom_spatraster(data=BaseRast)+
    scale_fill_viridis(na.value="white", option="magma")+
    theme(panel.grid = element_blank(), 
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank(),
          plot.background=element_blank(),
          panel.background = element_rect(fill="white"),
          legend.key=element_blank(),
          legend.text = element_text(size=6),
          legend.justification =  "top",
          legend.spacing.x = unit(0.1, 'cm'),
          legend.position=RastLegend)
  if(colourscale==TRUE){
    PastPoints <- PastPoints + geom_spatvector(data=SpecPoints, aes(colour=ColorScale), size=1, stroke=0.8)+
      scale_colour_viridis()+
      theme(legend.position=PointLegend)
  } else {
    PastPoints <- PastPoints + geom_spatvector(data=SpecPoints, size=1, stroke=0.8)
  }
  return(PastPoints)
}





Poly <- vect("/Users/hannahwauchope/Dropbox/Work/Data/GISLayers/TM_WORLD_BORDERS_SIMPL-0/TM_WORLD_BORDERS_SIMPL-0.3.shp")
Poly <- Poly[Poly$REGION == 142,]
Poly <- project(Poly, "+proj=robin +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")

SamplePoints <- data.frame(lon=c(130,49,-53,74,167,-102), lat=c(-30,-50,20,60,50,0))
SamplePoints$Color <- c(1,5,3,3,5,2)
SamplePoints$PointCol <- as.character(SamplePoints$Color)

#SamplePoints <- vect(SamplePoints, geom=c("lon", "lat"), crs=WGSCRS)
ModernTemps <- rast(paste0(DataFP, "/Climate_CRU/BioClim/", "1901_1950", "BioClim.grd"))[[1]]
crs(ModernTemps) <- WGSCRS
AnimRast <- ModernTemps
AnimPoint <- SamplePoints
AnimPoly <- Poly

MakeThatMap <- function(AnimPoints=NULL, AnimPolys=NULL, AnimRasts=NULL, MapOrder=c("Rast", "Poly", "Points"),
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
  if(class(AnimPoint)!="SpatVector"){
    GeomLon <- names(AnimPoint)[names(AnimPoint) %in% c("x", "lon", "longitude", "Longitude", "Lon")]
    GeomLat <- names(AnimPoint)[names(AnimPoint) %in% c("y", "lat", "latitude", "Latitude", "Lat")]
    AnimPoint <- vect(AnimPoint, geom=c(GeomLon, GeomLat), crs=WGSCRS)
  }
  
  if(class(AnimPoly)!="SpatVector"){
    AnimPoly <- vect(AnimPoly)
  }
  
  if(class(AnimRast)!="SpatRast"){
    AnimRast <- rast(AnimRast)
  }
  
  while(length(MapOrder) != 3){
    MapOrder <- c(MapOrder, "NULL")
  }

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
    AnimImage <- AnimImage + geom_spatvector(data=AnimPoly, colour=PolyLineCol, fill=PolyFillCol, alpha=PolyAlpha)
  } else if (MapOrder[[1]] == "Points"){
    AnimImage <- AnimImage + geom_spatvector(data=AnimPoint, aes(colour=PointCol), size=PointSize)
  }
  
  if(MapOrder[[2]]=="Rast"){
    AnimImage <- AnimImage + geom_spatraster(data=AnimRast)
  } else if (MapOrder[[2]]=="Poly"){
    AnimImage <- AnimImage + geom_spatvector(data=AnimPoly, colour=PolyLineCol, fill=PolyFillCol, alpha=PolyAlpha)
  } else if (MapOrder[[2]] == "Points"){
    AnimImage <- AnimImage + geom_spatvector(data=AnimPoint, aes(colour=PointCol), size=PointSize)
  }
  
  if(MapOrder[[3]]=="Rast"){
    AnimImage <- AnimImage + geom_spatraster(data=AnimRast)
  } else if (MapOrder[[3]]=="Poly"){
    AnimImage <- AnimImage + geom_spatvector(data=AnimPoly, colour=PolyLineCol, fill=PolyFillCol, alpha=PolyAlpha)
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
  
  if(is.null(RastMin)){RastMin <- minmax(AnimRast)[[1]]}
  if(is.null(RastMax)){RastMax <- minmax(AnimRast)[[2]]}
  
  if("Rast" %in% MapOrder){
    AnimImage <- AnimImage + scale_fill_viridis(na.value=NA, direction=1, name=RastColName, limits=c(RastMin, RastMax))
  }
  
  AnimImage <- AnimImage + guides(fill = guide_legend(order = 1), 
                                  colour = guide_legend(order = 2))
   
  return(AnimImage)
}











