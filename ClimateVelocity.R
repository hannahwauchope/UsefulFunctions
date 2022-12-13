library(CircStats)

#Calculate climate velocity. The concept of this all is taken from Burrows et al. 2011. The pace of shifting climate in marine and terrestrial ecosystems. Science, 334, 652-655.
#{http://science.sciencemag.org/content/334/6056/652}
#Would HIGHLY recommend reading this before proceeding! 
#Basically, spatial velocity calculates velocity in little blocks of 9 raster cells, and calculates how far you'd have to travel to remain in the same climate (by deg C)
#Temporal velocity calculates velocity by running many simple linear models through time for each grid cell in a raster stack, and calculates change per year. 
#Climate Velocity combines the two

### These functions are taken from package VoCC. Have modified one of the functions, so found it easier to just paste the relevant functions here. 
#BUT VoCC should always be cited for these: Garcia Molinos, J et al (2019) VoCC: An r package for calculating the velocity of climate change and related climatic metrics. Methods in Ecology and Evolution
#All except "tempTrendPerYear" are taken verbatim from VoCC 1.0.0. 
#Alternatively, just install VoCC for these (but note that TempTrendPerYear is still required from this code)

library(devtools)
devtools::install_github("JorGarMol/VoCC", dependencies = TRUE, build_vignettes = TRUE)

#The reason TempTrendPerYear is modified is that tempTrend gives you an estimate of change per time step (e.g. if you provide 5 rasters, it'll give you change per raster)
#Another way of putting this is that it assumes the rasters you provide are in consecutive years, e.g. years 1,2,3,4,5. 
#I modified this because I wanted to estimate *yearly* climate velocity, but from 5 rasters spanning 100 years. So I modified the formuala to give a model estimate per year, you just need to specify what the time steps between your rasters are
#Possibly all you'd need to do is divide the output of the original function by however many years each raster represents? But I can't be bothered to figure out the maths for that and this works. 

#spatGrad calculates the 'spatial velocity' (i.e. spatial gradient) and requires the function angulo
#tempTrend/tempTrendPerYear calculates the temporal velocity
#gVoCC brings the two together to give the overall velocity. It returns two rasters - one for magnitude and one for direction (as a value out of 360deg)

#Internal function, ignore
angulo <- function(dx, dy){
  d <- cbind(dx, dy)
  angline <- function(rw){
    angle <- ifelse(rw[2] < 0, 180 + CircStats::deg(atan(rw[1]/rw[2])),
                    ifelse(rw[1] < 0, 360 + CircStats::deg(atan(rw[1]/rw[2])), CircStats::deg(atan(rw[1]/rw[2]))))
    return(angle)
  }
  return(apply(d, 1, angline))
}

##Calculate spatial velocity
#r is the raster you want to use (e.g. if you're calculating velocity over the next 100 years it'd typically be the current raster)
#th is used to indicate a lower threshold to truncate the spatial gradient with. Use -Inf (default) if no threshold required.
#Projected indicates whether the raster is in a projected coordinate system or not (if not, it'll correct for latitudinal distortion)
#Output: a raster stack with the magnitude of spatial gradient (degC per km for unprojected rasters, or degC per spatial unit for projected) and angle (angle in degrees)
spatGrad <- function(r, th = -Inf, projected = FALSE){
  if(raster::nlayers(r) > 1){r <- raster::calc(r, mean, na.rm=T)}
  # get resolution of the raster
  re <- raster::res(r)
  # Create a columns for focal and each of its 8 adjacent cells
  y <- data.table(raster::adjacent(r, 1:ncell(r), directions = 8, pairs = TRUE))
  y <- na.omit(y[, climFocal := getValues(r)[from]][order(from, to)])   # Get value for focal cell, order the table by raster sequence and omit NAs (land cells)
  y[, clim := getValues(r)[to]] # Insert values for adjacent cells
  y[, sy := rowFromCell(r, from)-rowFromCell(r, to)]  # Column to identify rows in the raster (N = 1, mid = 0, S = -1)
  y[, sx := colFromCell(r, to)-colFromCell(r, from)]  # Same for columns (E = 1, mid = 0, W = -1)
  y[sx > 1, sx := -1]   # Sort out the W-E wrap at the dateline, part I
  y[sx < -1, sx := 1]   # Sort out the W-E wrap at the dateline, part II
  y[, code := paste0(sx, sy)]    # Make a unique code for each of the eight neighbouring cells
  # Code cells with positions
  y[.(code = c("10","-10","-11","-1-1","11","1-1","01","0-1"), to = c("climE","climW","climNW","climSW","climNE","climSE","climN","climS")), on = "code", code := i.to]
  y <- dcast(y[,c("from","code","clim")], from ~ code, value.var = "clim")
  y[, climFocal := getValues(r)[from]]   # Put climFocal back in
  y[, LAT := yFromCell(r, from)]         # Add focal cell latitude
  
  # Calculate individual spatial temperature gradients: grads (degC per km)
  # WE gradients difference in temperatures for each western and eastern pairs divided by the distance between the cells in each pair (corrected for  latitudinal distortion if unprojected)
  # Positive values indicate an increase in clim from W to E (i.e., in line with the Cartesian x axis)
  
  ifelse(projected == TRUE, d <- 1, d <- 111.325)
  ifelse(projected == TRUE, co <- 0, co <- 1)
  
  y[, gradWE1 := (climN-climNW)/(cos(co*CircStats::rad(LAT+re[2]))*(d*re[1]))]
  y[, gradWE2 := (climFocal-climW)/(cos(co*CircStats::rad(LAT))*(d*re[1]))]
  y[, gradWE3 := (climS-climSW)/(cos(co*CircStats::rad(LAT-re[2]))*(d*re[1]))]
  y[, gradWE4 := (climNE-climN)/(cos(co*CircStats::rad(LAT+re[2]))*(d*re[1]))]
  y[, gradWE5 := (climE-climFocal)/(cos(co*CircStats::rad(LAT))*(d*re[1]))]
  y[, gradWE6 := (climSE-climS)/(cos(co*CircStats::rad(LAT-re[2]))*(d*re[1]))]
  
  # NS gradients difference in temperatures for each northern and southern pairs divided by the distance between them (111.325 km per degC *re[2] degC)
  # Positive values indicate an increase in sst from S to N (i.e., in line with the Cartesian y axis)
  y[, gradNS1 := (climNW-climW)/(d*re[2])]
  y[, gradNS2 := (climN-climFocal)/(d*re[2])]
  y[, gradNS3 := (climNE-climE)/(d*re[2])]
  y[, gradNS4 := (climW-climSW)/(d*re[2])]
  y[, gradNS5 := (climFocal-climS)/(d*re[2])]
  y[, gradNS6 := (climE-climSE)/(d*re[2])]
  
  # Calulate NS and WE gradients. NOTE: for angles to work (at least using simple positive and negative values on Cartesian axes), S-N & W-E gradients need to be positive)
  y[, WEgrad := apply(.SD, 1, function(x) stats::weighted.mean(x, c(1,2,1,1,2,1), na.rm = T)), .SDcols = 12:17]
  y[, NSgrad := apply(.SD, 1, function(x) stats::weighted.mean(x, c(1,2,1,1,2,1), na.rm = T)), .SDcols = 18:23]
  y[is.na(WEgrad) & !is.na(NSgrad), WEgrad := 0L]     # Where NSgrad does not exist, but WEgrad does, make NSgrad 0
  y[!is.na(WEgrad) & is.na(NSgrad), NSgrad := 0L]     # same the other way around
  
  # Calculate angles of gradients (degrees) - adjusted for quadrant (0 deg is North)
  y[, angle := angulo(.SD$WEgrad, .SD$NSgrad), .SDcols = c("WEgrad", "NSgrad")]
  
  # Calculate the vector sum of gradients (C/km)
  y[, Grad := sqrt(apply(cbind((y$WEgrad^2), (y$NSgrad^2)), 1, sum, na.rm = TRUE))]
  
  # Merge the reduced file back into the main file to undo the initial na.omit
  from <- data.table(1:ncell(r)) # Make ordered from cells
  y <- y[from]   # merge both
  
  rAng <- rGrad <- raster::raster(r)
  rAng[y$from] <- y$angle
  rGrad[y$from] <- y$Grad
  rGrad[rGrad[] < th] <- th
  output <- raster::stack(rGrad,rAng)
  names(output) <- c("Grad", "Ang")
  return(output)
  
} #CHECK NOTES ON THE THRESHOLD

##Calculate temporal velocity
#r is a raster stack of the set of temp variables you want to use (e.g. for years 1, 2, 3)
#th let's you set the minimum number of observations in any given cell (e.g. if year 2 had an NA you could set it so that 3 so that all 3 years are required to calculate velocity)
#Output: a raster stack containing trend in degC per year per grid cell, plus standard error and statistical significance
tempTrend <- function(r, th) {
  y <- getValues(r)
  ocean <- which(rowSums(is.na(y))!= ncol(y))    # remove land cells
  y <- t(y[ocean, ])
  N <- apply(y, 2, function(x) sum(!is.na(x)))
  ind <- which(N >= th)
  y <- y[,ind]  # drop cells with less than th observations
  N <- apply(y, 2, function(x) sum(!is.na(x)))
  x <- matrix(nrow = nlayers(r), ncol = ncol(y))
  x[] <- 1:nlayers(r)
  # put NA values into the x values so they correspond with y
  x1 <- y
  x1[!is.na(x1)] <- 1
  x <- x*x1
  # calculate the sum terms
  sx <- apply(x, 2, sum, na.rm = T)
  sy <- apply(y, 2, sum, na.rm = T)
  sxx <- apply(x, 2, function(x) sum(x^2, na.rm = T))
  syy <- apply(y, 2, function(x) sum(x^2, na.rm = T))
  xy <- x*y
  sxy <- apply(xy, 2, sum, na.rm = T)
  # Estimate slope coefficients and associated standard errors and p-values
  slope <- (N*sxy-(sx*sy))/(N*sxx-sx^2)
  sres <- (N*syy-sy^2-slope^2*(N*sxx-sx^2))/(N*(N-2))
  SE <- suppressWarnings(sqrt((N*sres)/(N*sxx-sx^2)))
  Test <- slope/SE
  p <- mapply(function(x,y) (2*pt(abs(x), df = y-2, lower.tail = FALSE)), x = Test, y = N)
  
  slpTrends <- sigTrends <- seTrends <- raster(r[[1]])
  slpTrends[ocean[ind]] <- slope
  seTrends[ocean[ind]] <- SE
  sigTrends[ocean[ind]] <- p
  output <- stack(slpTrends,seTrends,sigTrends)
  names(output) <- c("slpTrends", "seTrends", "sigTrends")
  return(output)
}

#This is my modified function, with one extra parameter:
#xvals is a vector of the same length as the raster stack, giving the year values of those rasters
#e.g. if you provide 5 rasters for 2000,2020,2040,2060,2080 xvals would be c(0,20,40,60,80)
tempTrendPerYear <- function(r, th, xvals) {
  y <- getValues(r)
  ocean <- which(rowSums(is.na(y))!= ncol(y))    # remove land cells
  y <- t(y[ocean, ])
  N <- apply(y, 2, function(x) sum(!is.na(x)))
  ind <- which(N >= th)
  y <- y[,ind]  # drop cells with less than th observations
  N <- apply(y, 2, function(x) sum(!is.na(x)))
  x <- matrix(nrow = nlayers(r), ncol = ncol(y))
  if(nlayers(r)!=length(xvals)){stop("number of xvals and r layers are different")}
  x[] <- xvals
  # put NA values into the x values so they correspond with y
  x1 <- y
  x1[!is.na(x1)] <- 1
  x <- x*x1
  # calculate the sum terms
  sx <- apply(x, 2, sum, na.rm = T)
  sy <- apply(y, 2, sum, na.rm = T)
  sxx <- apply(x, 2, function(x) sum(x^2, na.rm = T))
  syy <- apply(y, 2, function(x) sum(x^2, na.rm = T))
  xy <- x*y
  sxy <- apply(xy, 2, sum, na.rm = T)
  # Estimate slope coefficients and associated standard errors and p-values
  slope <- (N*sxy-(sx*sy))/(N*sxx-sx^2)
  sres <- (N*syy-sy^2-slope^2*(N*sxx-sx^2))/(N*(N-2))
  SE <- suppressWarnings(sqrt((N*sres)/(N*sxx-sx^2)))
  Test <- slope/SE
  p <- mapply(function(x,y) (2*pt(abs(x), df = y-2, lower.tail = FALSE)), x = Test, y = N)
  
  slpTrends <- sigTrends <- seTrends <- raster(r[[1]])
  slpTrends[ocean[ind]] <- slope
  seTrends[ocean[ind]] <- SE
  sigTrends[ocean[ind]] <- p
  output <- stack(slpTrends,seTrends,sigTrends)
  names(output) <- c("slpTrends", "seTrends", "sigTrends")
  return(output)
}

##Calculate climate velocity
#This function takes the output from spatGrad and tempTrend (or tempTrendPerYear)
#It outputs a stack with the velocity (magnitude) and angle (degrees/360)
gVoCC <- function(tempTrend, spatGrad){
  VoCC <- tempTrend[[1]]/spatGrad[[1]]
  # velocity angles have opposite direction to the spatial climatic gradient if warming and same direction (cold to warm) if cooling
  ind <- which(getValues(VoCC) > 0)
  VoCCang <- spatGrad[[2]]
  VoCCang[ind] <- spatGrad[[2]][ind] + 180
  VoCCang[] <- ifelse(VoCCang[] >= 360, VoCCang[] - 360, VoCCang[])
  output <- stack(VoCC,VoCCang)
  names(output) <- c("voccMag", "voccAng")
  return(output)
}


