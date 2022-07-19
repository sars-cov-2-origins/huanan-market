# https://doi.org/10.5281/zenodo.6786454

library(MASS)
library(ks)
library(sp)
library(graphics)
library(scales)
library(sf)                                                                                                                               

# get a list of polygons that make up a particular probability contour
getPolygons <- function(kde,probability) {
  contourLevel <- contourLevels(kde, prob=(1-probability))
  contourLines <- contourLines(kde$eval.points[[1]], kde$eval.points[[2]], kde$estimate, level=contourLevel)
  polygons <- list()
  for (j in 1:length(contourLines)){
    polygons[[j]] <- Polygon(cbind(contourLines[[j]]$x,contourLines[[j]]$y))
  }
  return(polygons)
}

# get a specific probability contour from a kde
getContour <- function(kde, probability) {
  polygons <- getPolygons(kde, probability) 
  return(SpatialPolygons(list(Polygons(polygons,1))))
}

# iterate over probability list and create contour geojson files
writeContours <- function(kde, probabilities, filename) {
  for (p in probabilities) {
    contour <- getContour(kde,p)
    spdf <- SpatialPolygonsDataFrame(contour, data.frame(ID=1:length(contour)))
    #plot(spdf, border="gray30", add=T) # quick visual check
    sf <- st_as_sf(spdf)
    st_write(sf, paste(filename, "KDE", p, "geojson", sep = "."), append=FALSE, quiet=TRUE, driver="geoJson")
  }
}

# start of main 

# location of Huanan Market
market.x <- 114.25738
market.y <- 30.61955

# size of KDE grid
grid_size = 5000

# process December case list
filename <- "who_cases_dec-2019" 

cases <- read.csv(paste(filename, "csv", sep = "."),header=T)
cases.linked <- cases[cases$huanan_linked == "TRUE",]
cases.notLinked <- cases[cases$huanan_linked == "FALSE",]

longLat.all <- cbind(cases$longitude,cases$latitude)
longLat.linked <- cbind(cases.linked$longitude,cases.linked$latitude)
longLat.notLinked <- cbind(cases.notLinked$longitude,cases.notLinked$latitude)

casesets <- list("all", "linked", "notLinked")
probabilities <- list(0.50, 0.25, 0.1, 0.05, 0.01)

for (set in casesets) {
  if (set == "all") {
    cases.loc <- longLat.all
  } else if (set == "linked") {
    cases.loc <- longLat.linked
  } else if (set == "notLinked") {
    cases.loc <- longLat.notLinked
  }
  print(paste(filename, " dataset: ", set))

  kde <- kde(cases.loc, H=Hpi(cases.loc), compute.cont=T, gridsize=c(grid_size,grid_size))

  writeContours(kde, probabilities, paste(filename, set, sep = "."))
}

