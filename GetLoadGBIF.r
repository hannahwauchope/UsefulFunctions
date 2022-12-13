#Required libraries
library(rgbif)
library(pbapply)

#This function downloads species occurrence data from GBIF, and then reads it in ready for use in R

#The first function, "GetGBIFData", creates a download request, and then checks every 3 minutes to see if GBIF has procesed the request
#Once it has, it downloads the data and saves it

#You first need to create a gbif account at gbif.org, and then enter your username, password, and email into the function

#SpecNames should be a vector of species names (formatted as, e.g. "homo sapiens" or "Homo sapiens", no underscore separating the words)
#SaveFP should be where downloaded records should be saved
#Overwrite is for whether you want already downloaded records to be redownloaded

#If the species name doesn't exist in GBIF, it'll tell you the species isn't found in GBIF. 

#The function will return a dataframe, with a row for every species in SpecNames, and a column detailing what happened either "Downloaded", "Already Downloaded", or "Not found in GBIF"

GetGBIFData <- function(SpecNames, SaveFP, GBIFuser, GBIFpwd, GBIFemail, ReDownload=FALSE, cluster=FALSE){
  return(do.call("rbind", pblapply(SpecNames, function(i){
    if(ReDownload==FALSE){
      if(length(list.files(paste0(SaveFP, i), full.names=TRUE, recursive=TRUE))>0){
        return(data.frame(Species = i, Status = "Already downloaded"))
      }
    }
    if(cluster==TRUE){
      dir.create(paste0(SaveFP, "LogBook/"), showWarnings = FALSE)
      if(file.exists(paste0(SaveFP, "LogBook/", i, " Begin.csv"))){
        return(data.frame(Species = i, Status = "Another core's on it"))
      } else {
        write.csv(NULL, paste0(SaveFP, "LogBook/", i, " Begin.csv"))
      }
    }
    dir.create(paste0(SaveFP, i), showWarnings = FALSE)
    unlink(list.files(paste0(SaveFP, i), full.names=TRUE, recursive = TRUE), recursive = TRUE)
    TaxKey <- name_suggest(q=i, rank="species")$data$key #First get the "TaxonKey" (each species has its own key id)
    if(is.null(TaxKey)){
      write.csv(NULL, paste0(SaveFP, "LogBook/", i, " NotFoundinGBIF.csv"))
      return(data.frame(Species = i, Status = "Not found in GBIF"))
    }
    GBIFDownload <- occ_download(pred_in("taxonKey", TaxKey), user=GBIFuser, pwd=GBIFpwd, email=GBIFemail)
    print(paste0("Begin downloading ", i))
    while(tryCatch(occ_download_get(GBIFDownload, path = paste0(SaveFP, i), overwrite = TRUE), error=function(e) "ERROR") == "ERROR"){
      Sys.sleep(180)
    }
    print(paste0("Downloaded ", i))
    write.csv(NULL, paste0(SaveFP, "LogBook/", i, " Done.csv"))
    return(data.frame(Species = i, Status = "Downloaded"))
  })))
}
#GetGBIFData(SpecNames, SaveFP, GBIFuser, GBIFpwd, GBIFemail)
## This function loads in the downloaded GBIF data and cleans it up a bit.
#It needs SpecNames and SaveFP as per above
#Columns to keep should be a vector of which column names you want. Either "LocationAndTaxa" (default, the ones I tend to keep), "All", or your own defined vector
  ### The downloaded data has 257 columns (!) (at time of writing, Nov 22). An option to figure out which ones you want is to just enter one species as your "SpecName" with ColumnsToKeep="All", just to examine the columns and then decide
  ### It's all a bit confusing to get info on the columns, are bunch of them are defined here under queries (https://www.gbif.org/developer/occurrence#p_basisOfRecord), otherwise google gbif API and the column name and you should find it eventually
#RemoveLocationsWithNoCoordinates ditches locations that are false for "hasCoordinate" and true for "hasGeospatialIssues"

LoadGBIFData <- function(SpecNames, SaveFP, GBIFuser, GBIFpwd, GBIFemail, ColumnsToKeep = "LocationAndTaxa", RemoveLocationsWithNoCoordinates=TRUE, redownload = TRUE){
  ReadGBIFDat <- pblapply(SpecNames, function(i){
    print(paste0("Attempt to read ", i))
    if(length(list.files(paste0(SaveFP, i)))<2){
      zipfile <- list.files(paste0(SaveFP, i), pattern="*.zip", full.names=TRUE, recursive=TRUE)
      if(length(zipfile)==0){
        print(paste0(i, " has not been downloaded from GBIF"))
        return(NULL)
      }
      UnzipTry <- tryCatch(unzip(zipfile, exdir=paste0(SaveFP, i)), warning=function(e) e)
      if(UnzipTry[[1]] == "error 1 in extracting from zip file" & redownload == TRUE){
        print(paste0(i, "'s file is corrupted, attempting redownload"))
        unlink(zipfile)
        GetGBIFData(i, SaveFP, GBIFuser, GBIFpwd, GBIFemail)
        zipfile <- list.files(paste0(SaveFP, i), pattern="*.zip", full.names=TRUE, recursive=TRUE)
        UnzipTry <- tryCatch(unzip(zipfile, exdir=paste0(SaveFP, i)), warning=function(e) e)
        if(UnzipTry[[1]] == "error 1 in extracting from zip file"){
          return(paste0("something's fucked with ", i))
        }
      } else if (UnzipTry[[1]] == "error 1 in extracting from zip file" & redownload == FALSE){
        print(paste0(i, "'s file is corrupted and redownload is set to false"))
        return(NULL)
      }
    }
    GBIFDat <- read.csv(paste0(SaveFP, i, "/occurrence.txt"), header = TRUE, sep="\t", fill = TRUE, quote = "")
    if(ColumnsToKeep[[1]] == "LocationAndTaxa"){
      ColumnsToKeep <- c("gbifID", "spatial", "temporal", "basisOfRecord", "occurrenceID", "sex", "lifeStage", "occurrenceStatus", "occurrenceRemarks", "organismID", 
      "year", "month", "day", "verbatimEventDate", "locality", "verbatimLocality", "decimalLatitude", "decimalLongitude", "coordinateUncertaintyInMeters", "coordinatePrecision", "verbatimCoordinateSystem", "taxonID",
      "scientificNameID", "scientificName", "kingdom", "phylum", "class", "order", "family", "genus", "subgenus", "specificEpithet", "infraspecificEpithet", "hasGeospatialIssues", "hasCoordinate", "species", "genericName","verbatimScientificName", "iucnRedListCategory")
    }
    if(ColumnsToKeep[[1]] == "All"){
      ColumnsToKeep <- c(1:ncol(GBIFDat))
    }
    if(length(ColumnsToKeep[!ColumnsToKeep %in% colnames(GBIFDat)]) != 0){
      stop("Some of the names in ColumnsToKeep don't exist as column names in the download")
    }
    GBIFDat <- GBIFDat[,c(ColumnsToKeep)]
    if(RemoveLocationsWithNoCoordinates==TRUE){
      #For some dumb reason GBIF has started quoting it's true/false logicals as characters, but I don't know if they'll switch back, so this is clunky to make sure the function works
      if(class(unique(GBIFDat$hasGeospatialIssues)[[1]])=="logical"){
        GBIFDat <- GBIFDat[GBIFDat$hasCoordinate,]
        GBIFDat <- GBIFDat[!GBIFDat$hasGeospatialIssues,]
      } else {
        GBIFDat <- subset(GBIFDat, hasCoordinate=="true")
        GBIFDat <- subset(GBIFDat, hasGeospatialIssues=="false")
      }
    }
    if(nrow(GBIFDat)==0){
      print(paste0(i, " has no data after subsetting"))
      return(NULL)
    }
    GBIFDat$GivenSpec <- i
    print(paste0(i, " has been read in, success!"))
    return(GBIFDat)
  }) #
  return(do.call("rbind", ReadGBIFDat))
}
#LoadGBIFData(SpecNames, SaveFP)

