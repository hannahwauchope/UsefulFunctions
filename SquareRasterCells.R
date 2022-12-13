#Square rast cells
#With thanks to adamlilith/enmSdm https://rdrr.io/github/adamlilith/enmSdm/src/R/squareRastCells.r
squareRastCells <- function(x, keepWidth = TRUE, ...) {
  
  OGcrs <- terra::crs(x)
  OGext <- terra::ext(x)
  OGext <- c(OGext[1], OGext[2], OGext[3], OGext[4])
  
  if (keepWidth) {
    
    # pad extent so it can accommodate all necessary rows
    width <- res(x)[1]
    northSouth <- OGext[4] - OGext[3]
    ncols <- ncol(x)
    
    nrows <- ceiling(northSouth / width)
    crammedRows <- northSouth / width
    pad <- (width * (nrows - crammedRows)) / 2
    OGext <- c(
      OGext[1],
      OGext[2],
      OGext[3] - pad,
      OGext[4] + pad
    )
    
  } else if (!keepWidth) {
    
    # pad extent so it can accommodate all necessary columns
    height <- res(x)[2]
    eastWest <- OGext[2] - OGext[1]
    nrows <- nrow(x)
    
    ncols <- ceiling(eastWest / height)
    crammedCols <- eastWest / height
    pad <- (height * (ncols - crammedCols)) / 2
    OGext <- c(
      OGext[1] - pad,
      OGext[2] + pad,
      OGext[3],
      OGext[4]
    )
    
  }
  
  NEWext <- terra::ext(OGext)
  
  template <- terra::rast(NEWext, nrows=nrows, ncols=ncols, crs=OGcrs)
  rast <- terra::resample(x, template, ...)
  rast
  
} 
