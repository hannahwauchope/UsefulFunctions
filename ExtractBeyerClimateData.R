### A function to extract rasters of bioclimatic climate data from Beyer et al "High-resolution terrestrial climate, bioclimate and vegetation for the last 120,000 years" https://www.nature.com/articles/s41597-020-0552-1
#It will return a list of raster stacks, one stack for every year of data requested. Each stack contains a raster for each bioclimatic variable 

#The dataset offers numerous climate options but this function just returns bioclimatic variables

#If BBox is used, this should be an extent object - it will cut the rasters down from global to the extent of BBox
#If you set "TellMeYearOptions" to TRUE, the function will just output which years you can get data for - this is useful to know what to input to reqyears
#use reqyears to specify which years you want. Either "All" for all years, or eg c("0", "1000", "2000") for thos years
#Use EqualAr if you want the raster converted to mollweide equal area projection - defaults to false

#First, set the file path for the netcdffile
library(ncdf4)


GetPastClimateBeyer <- function(NetCDFFilePath, BBox = NULL, reqyears = "All", EqualAr = FALSE, TellMeYearOptions = F){
  WGSCRS <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
  MollCRS <- CRS("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
  RobCRS <- CRS("+proj=robin +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
  
  env_nc      <- ncdf4::nc_open(NetCDFFilePath) #open the netcdf
  longitude   <- ncdf4::ncvar_get(env_nc, "longitude") #get coordinate data
  latitude    <- ncdf4::ncvar_get(env_nc, "latitude") #get coordinate data
  years <- ncdf4::ncvar_get(env_nc, "time") #get the avaliable years
  
  if(reqyears[[1]]=="All"){ #if all years required, specify
    reqyears <- years
  }
  
  if(TellMeYearOptions == TRUE){ #Report the years available
    ncdf4::nc_close(env_nc)
    return(years)
  }
  
  BioVars <- names(env_nc$var)[11:length(names(env_nc$var))] #Just get bioclim variables
  
  PalBIO <- pblapply(BioVars, function(BV){
    BIO <- ncdf4::ncvar_get(env_nc, BV)
    return(BIO)
  }) #extract bioclim data for all years
  ncdf4::nc_close(env_nc)
  
  lonlat <- as.matrix(expand.grid(longitude,latitude)) #Make a matrix of coordinates
  
  PalRasters <- pblapply(reqyears, function(YR){ #For each year in required years
    BioRas <- lapply(c(1:length(PalBIO)), function(BIO){ #For each bioclim variable
      BioYear <- PalBIO[[BIO]]
      BioYear <- BioYear[,,which.min(abs(years - YR))] #get the year
      BioYear <- as.vector(BioYear) 
      BioYear_df <- data.frame(cbind(lonlat,BioYear)) #get biodata
      names(BioYear_df) <- c("lon","lat","BIO1")
      BIO_ras <- rasterFromXYZ(BioYear_df)  #Convert first two columns as lon-lat and third as value                
      if(!is.null(BBox)){ #Crop (if bbox given)
        BIO_ras <- crop(BIO_ras, BBox)
      }
      crs(BIO_ras) <- WGSCRS
      if(EqualAr == T){ #Make equal area if desired
        BIO_ras <- projectRaster(BIO_ras, crs=MollCRS)
        BIO_ras <- squareRastCells(BIO_ras)
      }
      return(BIO_ras)
    })
    names(BioRas) <- BioVars
    return(stack(BioRas))
  })
  names(PalRasters) <- reqyears
  return(PalRasters)
}

NetCDFFilePath <- ("/Users/hannahwauchope/Dropbox/Work/Data/LateQuaternaryClimate_BeyerManica/LateQuaternary_Environment.nc")


#e.g. this returns just what years are available
GetPastClimateBeyer(NetCDFFilePath, TellMeYearOptions = T)
BeyerClim <- GetPastClimateBeyer(NetCDFFilePath, BBox = extent(ArcticBiomes), reqyears = "All", EqualAr = FALSE, TellMeYearOptions = F)




