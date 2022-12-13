#Get rasters of recent-historical to present day climate data from CRU TS4

#NOTE, THIS CODE IS UPDATED TO WORK WITH THE TERRA SPATIAL PACKAGE, *NOT* THE RASTER PACKAGE (which is going to die soon)

#Before you begin, you need to go to the CRU website, and download the big netcdf files:

#STEP 1.
#Go to the CRU TS4 website, currently hosted on CEDA, and find the latest version (As of writing [Nov 2022] it's CRU TS4.06, which goes up to 2021. It'll be 4.07 sometime in 2023, with data going up to 2022 etc)
#https://catalogue.ceda.ac.uk/uuid/e0b4e1e56c1c4460b796073a31366980 #This is currently the link that works

#STEP 2. register a free account

#STEP 3. Hit the download button, review the readme and release notes if you need, then click on the "data" folder
#You'll find a bunch of folders representing different climate variables here:
  #cloud cover (cld)
  #diurnal temperature range (dtr)
  #frost day frequency (frs)
  #precipitation (pre)
  #daily mean temperature (tmp)
  #monthly average daily maximum (tmx)
  #minimum (tmn) temperature
  #vapour pressure (vap)
  #Potential Evapo-transpiration (pet)
  #wet day frequency (wet).

#STEP 4. For all variables that you are interested in, click on the folder and download the BIG NETCDF file
#Unless they change the filename structure, it'll be called e.g. "cru_ts4.06.1901.2021.cld.dat.nc.gz" (if you've got the 2021 release, and want cloud cover (cld) data)
  #The important thing is its the 1901-202X data, and that the file name ends in dat.nc.gz

#STEP 5. Put the (one or multiple) netcdf files in a folder

#STEP 6. Run the below function to extract your data

library(ncdf4)
library(chron)
library(terra)

DataFP <- "/Users/hannahwauchope/Dropbox/Work/Data/"
filepath <- paste0(DataFP, "/Climate_CRU/ClimateData")

#Filepath is the path to where you netcdf files are saved
#months is a vector of the months you want data from, expressed as "01" (for January) to "12" (for December)
#years is the years you want data from, expressed as e.g. "1901" (if you want a range, enter it as "as.character(2001:2020)")
#climtype is the three letter shortenings of the climate type as per the file name (see above, e.g. "tmp" for daily mean temperature)
##To avoid the output getting very confusing, you need to run this function separately for each climtype you want data for

#The function will return all combos of the months and years you provide (e.g. if you give c("01", "02", "03") and c(as.character(2000-2010)), you'll get rasters for those three months for each of those years)
#UNLESS you make SpecifyCombos = TRUE, then months and years need to be the same length, and it will only give the combos as you provide them (so if you put c("01", "02", "03") and c("2001", "2002", "2003") you'd get Jan 2001, Feb 2002 and March 2003)

#You'll receive a list of rasters, where each raster (and list item, to be safe) is named as "month_year_climtype"

rastmonyear <- function(filepath, months, years, climtype, SpecifyCombos=FALSE){
  climate_output <- nc_open(list.files(path=filepath,pattern=paste0("*", climtype, ".dat.nc"), full.names=TRUE))
  
  #Get data
  lon <- ncvar_get(climate_output,"lon")
  lat <- ncvar_get(climate_output,"lat",verbose=F)
  time <- ncvar_get(climate_output,"time")
  tunits <- ncatt_get(climate_output,"time","units")
  fillvalue <- ncatt_get(climate_output, climtype,"_FillValue")
  clim_array <- ncvar_get(climate_output,climtype)
  clim_array[clim_array==fillvalue$value] <- NA #Change NAs to appropriate thing
  dlname <- ncatt_get(climate_output,climtype,"long_name")
  dunits <- ncatt_get(climate_output, climtype,"units")
  
  #Get metadata
  title <- ncatt_get(climate_output,0,"title")
  institution <- ncatt_get(climate_output,0,"institution")
  datasource <- ncatt_get(climate_output,0,"source")
  references <- ncatt_get(climate_output,0,"references")
  history <- ncatt_get(climate_output,0,"history")
  Conventions <- ncatt_get(climate_output,0,"Conventions")
  nc_close(climate_output)
  
  #Get the right times
  tustr <- strsplit(tunits$value, " ")
  tdstr <- strsplit(unlist(tustr)[3], "-")
  tmonth <- as.integer(unlist(tdstr)[2])
  tday <- as.integer(unlist(tdstr)[3])
  tyear <- as.integer(unlist(tdstr)[1])
  dates <- chron(time,origin=c(tmonth, tday, tyear))
  dates <- do.call("rbind", lapply(1:length(dates), function(x) as.data.frame(t(strsplit(as.character(dates[x]), "/")[[1]][c(1,3)]))))
  names(dates) <- c("Month", "Year")
  dates$RasNum <- rownames(dates)
  
  #Fix up the dates (they only use the second two numbers of the date, meaning currently 1901 and 2001 are the same. Add the 19/20 back in)
  dates[1:(99*12),]$Year <- paste0("19", dates[1:(99*12),]$Year)
  dates[((99*12)+1):nrow(dates),]$Year <- paste0("20", dates[((99*12)+1):nrow(dates),]$Year)
  
  #Extract raster of designated month/year
  if(SpecifyCombos == FALSE){
    ClimCombos <- expand.grid(months, years)
  } else {
    if(length(months)!=length(years)){stop("Length of months and years need to match if using SpecifyCombos=TRUE")}
    ClimCombos <- data.frame(Var1 = months, Var2 = years)
  }
  
  ClimRasters <- lapply(1:nrow(ClimCombos), function(x){
    month <- ClimCombos[x,1]
    year <- ClimCombos[x,2]
    m <- as.numeric(as.character(dates[dates$Month==month & dates$Year==year,]$RasNum))
    clim_slice <- clim_array[,,m]
    #Ok so that gets alll our data out
    # matrix (nlon*nlon rows by 2 cols) of lons and lats
    lonlat <- as.matrix(expand.grid(lon,lat))
    # vector of values
    clim_vec <- as.vector(clim_slice)
    clim_df <- data.frame(cbind(lonlat,clim_vec))
    names(clim_df) <- c("lon","lat",paste(climtype,as.character(m), sep="_"))
    clim_ras <- rast(clim_df, type="xyz")  #Convert first two columns as lon-lat and third as value
    names(clim_ras) <- paste(month, year, climtype, sep="_")
    crs(clim_ras) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
    return(clim_ras)
  })
  names(ClimRasters) <- sapply(1:nrow(ClimCombos), function(x) paste(ClimCombos[x,1], ClimCombos[x,2], climtype, sep="_"))
  return(ClimRasters)
}




