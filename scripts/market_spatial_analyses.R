# https://doi.org/10.5281/zenodo.6786454

library(ks)
library(sp)
library(sf)                                                                                                                               
library(geosphere)
library(raster)
library(ggplot2)
library(stats)

######## Functions ##########
getPolygons <- function(points,probability) {
  kde = kde(points, H=Hpi(points), compute.cont=T, gridsize=c(1000,1000))
  contourLevel = contourLevels(kde, prob=(1-probability))
  contourLines = contourLines(kde$eval.points[[1]], kde$eval.points[[2]], kde$estimate, level=contourLevel)
  polygons = list()
  for (j in 1:length(contourLines)){
    polygons[[j]] = Polygon(cbind(contourLines[[j]]$x,contourLines[[j]]$y))
  }
  return(polygons)
}

getPvalue <- function(metric,nullDis) {
  pValue <- sum(nullDis < metric)
  return(pValue/length(nullDis))
}

create.spPolygon <- function(x){
  p = SpatialPolygons(list(Polygons(list(Polygon(x)),1)))
  projection(p) <- CRS("+proj=longlat +datum=WGS84")
  p
}

create.circle.spPolygon <- function(long, lat, dist.m){
  bearing <- c(1:360)
  circle <- data.frame(bearing)
  
  circle[,2:3] <- destPoint(c(long, lat), circle$bearing, dist.m)
  circle$bearing <- NULL
  names(circle)[1] <- "lon"
  names(circle)[2] <- "lat"
  
  circle.p <- create.spPolygon(circle)
  circle.p
}

#fixed variables
market.lat <- 30.61955
market.long <- 114.25738
radius <- 1000
#radius <- 500
numberOfGhosts = 9
numberOfReps <- 1000
numberOfSubsampled <- 104
unsampled_KDEprobabily <- 0.25

# process December case list
data.path <- "../data/"
#filename1 <- paste(data.path,"who_cases_dec-2019" 
filename1 <- paste(data.path,"who_cases_dec-2019", sep = "")
cases = read.csv(paste(filename1, "csv", sep = "."),header=T)
cases.linked <- cases[cases$huanan_linked == "TRUE",]
cases.notLinked <- cases[cases$huanan_linked == "FALSE",]
cases.linked.B <- cases.linked[cases.linked$lineage == "B",]
cases.notLinked.B <- cases.notLinked[cases.notLinked$lineage == "B",]
cases.notLinked.A <- cases.notLinked[cases.notLinked$lineage == "A",]
cases.notLinked.A1 <- cases.notLinked.A
cases.B <- cases[cases$lineage == "B",]
cases.confirmed <- cases[cases$confirmed == "TRUE" & !is.na(cases$confirmed),]
cases.clinDiag <- cases[cases$confirmed == "FALSE" & !is.na(cases$confirmed),]
#adding one case not linked to the market known to have been infected with lineage A and known to have stayed in a nearby hotel for the 5 days before symptom onset
additionalCase <- c(NA, NA, NA, NA, NA, FALSE, FALSE, "A", NA, "additional case not linked to the market known to have been infected with lineage A and known to have stayed in a nearby hotel")
cases.notLinked.A[nrow(cases.notLinked.A) + 1,] <- additionalCase

#read in null distributions
median.distances.all.weibo <- read.csv2(paste(data.path,"median_distance_weibo_null_155.csv", sep = ""),header=T, sep = ",")
median.distances.notLinked.weibo <- read.csv2(paste(data.path,"median_distance_weibo_null_120.csv", sep = ""),header=T, sep = ",")
median.distances.linked.weibo <- read.csv2(paste(data.path,"median_distance_weibo_null_35.csv", sep = ""),header=T, sep = ",")
median.distances.linked.B.weibo <- read.csv2(paste(data.path,"median_distance_weibo_null_10.csv", sep = ""),header=T, sep = ",")
median.distances.notLinked.B.weibo <- read.csv2(paste(data.path,"median_distance_weibo_null_1.csv", sep = ""),header=T, sep = ",")
median.distances.notLinked.A1.weibo <- median.distances.notLinked.B.weibo
median.distances.notLinked.A.weibo <- read.csv2(paste(data.path,"median_distance_weibo_null_2.csv", sep = ""),header=T, sep = ",")
median.distances.B.weibo <- read.csv2(paste(data.path,"median_distance_weibo_null_11.csv", sep = ""),header=T, sep = ",")
median.distances.confirmed.weibo <- read.csv2(paste(data.path,"median_distance_weibo_null_65.csv", sep = ""),header=T, sep = ",")
median.distances.clinDiag.weibo <- read.csv2(paste(data.path,"median_distance_weibo_null_79.csv", sep = ""),header=T, sep = ",")

median.distances.notLinked.sub.weibo <- read.csv2(paste(data.path,"median_distance_weibo_null_104.csv", sep = ""),header=T, sep = ",")
median.distances.all.ghosts.weibo <- read.csv2(paste(data.path,"median_distance_weibo_null_164.csv", sep = ""),header=T, sep = ",")
median.distances.notLinked.ghosts.weibo <- read.csv2(paste(data.path,"median_distance_weibo_null_129.csv", sep = ""),header=T, sep = ",")
median.distances.linked.ghosts.weibo <- read.csv2(paste(data.path,"median_distance_weibo_null_44.csv", sep = ""),header=T, sep = ",")

centerpoint.distances.weibo <- read.csv2(paste(data.path,"Weibo_1000_points_sampled_with_replacement_distance_to_huanan.csv", sep = ""),header=T, sep = ",")

median.distances.all.worldpop.age <- read.csv2(paste(data.path,"distance_popdensityagegroups_null_155.csv", sep = ""),header=T, sep = ",")
median.distances.notLinked.worldpop.age <- read.csv2(paste(data.path,"distance_popdensityagegroups_null_120.csv", sep = ""),header=T, sep = ",")
median.distances.linked.worldpop.age <- read.csv2(paste(data.path,"distance_popdensityagegroups_null_35.csv", sep = ""),header=T, sep = ",")
median.distances.linked.B.worldpop.age <- read.csv2(paste(data.path,"distance_popdensityagegroups_null_10.csv", sep = ""),header=T, sep = ",")
median.distances.notLinked.B.worldpop.age <- read.csv2(paste(data.path,"distance_popdensityagegroups_null_1.csv", sep = ""),header=T, sep = ",")
median.distances.notLinked.A1.worldpop.age <- median.distances.notLinked.B.worldpop.age
median.distances.notLinked.A.worldpop.age <- read.csv2(paste(data.path,"distance_popdensityagegroups_null_2.csv", sep = ""),header=T, sep = ",")
median.distances.B.worldpop.age <- read.csv2(paste(data.path,"distance_popdensityagegroups_null_11.csv", sep = ""),header=T, sep = ",")
median.distances.confirmed.worldpop.age <- read.csv2(paste(data.path,"distance_popdensityagegroups_null_65.csv", sep = ""),header=T, sep = ",")
median.distances.clinDiag.worldpop.age <- read.csv2(paste(data.path,"distance_popdensityagegroups_null_79.csv", sep = ""),header=T, sep = ",")

median.distances.notLinked.sub.worldpop.age <- read.csv2(paste(data.path,"distance_popdensityagegroups_null_104.csv", sep = ""),header=T, sep = ",")
median.distances.all.ghosts.worldpop.age <- read.csv2(paste(data.path,"distance_popdensityagegroups_null_164.csv", sep = ""),header=T, sep = ",")
median.distances.notLinked.ghosts.worldpop.age <- read.csv2(paste(data.path,"distance_popdensityagegroups_null_129.csv", sep = ""),header=T, sep = ",")
median.distances.linked.ghosts.worldpop.age <- read.csv2(paste(data.path,"distance_popdensityagegroups_null_44.csv", sep = ""),header=T, sep = ",")

filename2 <- paste(data.path,"centroidnull_agegroups_1mill", sep = "") 
centerpoint.distances.worldpop.age = read.csv(paste(filename2, "csv", sep = "."),header=T)

######## median distance approach, all cases, not linked, linked, linked B, not linked A, not linked B ######## 
distances.toMarket <- vector(mode="numeric", length=length(cases[,1]))
distances.toMarket.notLinked <- vector(mode="numeric", length=length(cases.notLinked[,1]))
distances.toMarket.linked <- vector(mode="numeric", length=length(cases.linked[,1]))
distances.toMarket.linked.B <- vector(mode="numeric", length=length(cases.linked.B[,1]))
distances.toMarket.notLinked.B <- vector(mode="numeric", length=length(cases.notLinked.B[,1]))
distances.toMarket.notLinked.A1 <- vector(mode="numeric", length=length(cases.notLinked.A1[,1]))
distances.toMarket.notLinked.A <- vector(mode="numeric", length=length(cases.notLinked.A[,1]))
distances.toMarket.B <- vector(mode="numeric", length=length(cases.B[,1]))
distances.toMarket.confirmed <- vector(mode="numeric", length=length(cases.confirmed[,1]))
distances.toMarket.clinDiag <- vector(mode="numeric", length=length(cases.clinDiag[,1]))

notLinked.count <- 1
linked.count <- 1
linked.count.B <- 1
notLinked.count.B <- 1
notLinked.count.A <- 1
count.B <- 1
confirmed.count <- 1
clinDiag.count <- 1

for (i in 1:length(distances.toMarket)){
  dist<- (distVincentyEllipsoid(c(market.long,market.lat), c(cases$long[i],cases$lat[i])))/1000
  #print(dist)
  distances.toMarket[i] <- dist
  if (cases$huanan_linked[i] == FALSE){
    distances.toMarket.notLinked[notLinked.count] <- dist
    notLinked.count <- notLinked.count +1
    if (cases$lineage[i] == "A"){
      distances.toMarket.notLinked.A[notLinked.count.A] <- dist
      distances.toMarket.notLinked.A1[notLinked.count.A] <- dist
      notLinked.count.A <- notLinked.count.A +1
    } else if (cases$lineage[i] == "B"){
      distances.toMarket.notLinked.B[notLinked.count.B] <- dist
      notLinked.count.B <- notLinked.count.B +1
    }
  } else {
    distances.toMarket.linked[linked.count] <- dist
    linked.count <- linked.count +1
    if (cases$lineage[i] == "B"){
      distances.toMarket.linked.B[linked.count.B] <- dist
      linked.count.B <- linked.count.B +1
    }
  }
  if (cases$lineage[i] == "B"){
    distances.toMarket.B[count.B] <- dist
    count.B <- count.B + 1
  }
  if (cases$confirmed[i] == "TRUE" & !is.na(cases$confirmed[i])){
    distances.toMarket.confirmed[confirmed.count] <- dist
    confirmed.count <- confirmed.count + 1
  } else if (cases$confirmed[i] == "FALSE" & !is.na(cases$confirmed[i])){
    distances.toMarket.clinDiag[clinDiag.count] <- dist
    clinDiag.count <- clinDiag.count + 1
  }
}

#adding the distance for the one case not linked to the market known to have been infected with lineage A and known to have stayed in a nearby hotel. This hotel is much closer to the market than the residence of the other A case (=2.31km) 
dist.notLinked.A <- 2.31
distances.toMarket.notLinked.A[2] <- dist.notLinked.A

#t.test(distances.toMarket.notLinked, y = distances.toMarket.linked)
shapiro.test(distances.toMarket.notLinked)
shapiro.test(distances.toMarket.linked)
#unpaired two-samples Wilcoxon test 
wilcox.test(distances.toMarket.notLinked, distances.toMarket.linked, alternative = "two.sided")
wilcox.test(distances.toMarket.notLinked.A, distances.toMarket.B, alternative = "two.sided")


median.dist.all.pValue.weibo <- getPvalue(median(distances.toMarket),as.numeric(median.distances.all.weibo$Median_distance_huanan))
median.dist.notLinked.pValue.weibo <- getPvalue(median(distances.toMarket.notLinked),as.numeric(median.distances.notLinked.weibo$Median_distance_huanan))
median.dist.linked.pValue.weibo <- getPvalue(median(distances.toMarket.linked),as.numeric(median.distances.linked.weibo$Median_distance_huanan))
median.dist.linked.B.pValue.weibo <- getPvalue(median(distances.toMarket.linked.B),as.numeric(median.distances.linked.B.weibo$Median_distance_huanan))
median.dist.notLinked.B.pValue.weibo <- getPvalue(median(distances.toMarket.notLinked.B),as.numeric(median.distances.notLinked.B.weibo$Median_distance_huanan))
median.dist.notLinked.A1.pValue.weibo <- getPvalue(median(distances.toMarket.notLinked.A1),as.numeric(median.distances.notLinked.A1.weibo$Median_distance_huanan))
median.dist.notLinked.A.pValue.weibo <- getPvalue(median(distances.toMarket.notLinked.A),as.numeric(median.distances.notLinked.A.weibo$Median_distance_huanan))
median.dist.B.pValue.weibo <- getPvalue(median(distances.toMarket.B),as.numeric(median.distances.notLinked.A.weibo$Median_distance_huanan))
median.dist.confirmed.pValue.weibo <- getPvalue(median(distances.toMarket.confirmed),as.numeric(median.distances.confirmed.weibo$Median_distance_huanan))
median.dist.clinDiag.pValue.weibo <- getPvalue(median(distances.toMarket.clinDiag),as.numeric(median.distances.clinDiag.weibo$Median_distance_huanan))

median.dist.all.pValue.worldpop.age <- getPvalue(median(distances.toMarket),as.numeric(median.distances.all.worldpop.age$Median_distance_huanan))
median.dist.notLinked.pValue.worldpop.age <- getPvalue(median(distances.toMarket.notLinked),as.numeric(median.distances.notLinked.worldpop.age$Median_distance_huanan))
median.dist.linked.pValue.worldpop.age <- getPvalue(median(distances.toMarket.linked),as.numeric(median.distances.linked.worldpop.age$Median_distance_huanan))
median.dist.linked.B.pValue.worldpop.age <- getPvalue(median(distances.toMarket.linked.B),as.numeric(median.distances.linked.B.worldpop.age$Median_distance_huanan))
median.dist.notLinked.B.pValue.worldpop.age <- getPvalue(median(distances.toMarket.notLinked.B),as.numeric(median.distances.notLinked.B.worldpop.age$Median_distance_huanan))
median.dist.notLinked.A1.pValue.worldpop.age <- getPvalue(median(distances.toMarket.notLinked.A1),as.numeric(median.distances.notLinked.A1.worldpop.age$Median_distance_huanan))
median.dist.notLinked.A.pValue.worldpop.age <- getPvalue(median(distances.toMarket.notLinked.A),as.numeric(median.distances.notLinked.A.worldpop.age$Median_distance_huanan))
median.dist.B.pValue.worldpop.age <- getPvalue(median(distances.toMarket.B),as.numeric(median.distances.B.worldpop.age$Median_distance_huanan))
median.dist.confirmed.pValue.worldpop.age <- getPvalue(median(distances.toMarket.confirmed),as.numeric(median.distances.confirmed.worldpop.age$Median_distance_huanan))
median.dist.clinDiag.pValue.worldpop.age <- getPvalue(median(distances.toMarket.clinDiag),as.numeric(median.distances.clinDiag.worldpop.age$Median_distance_huanan))

######## centerpoint distance approach, 155 cases, not linked, linked, linked B, not linked A, not linked B ######## 
centerpoint.longitude <- median(cases$long)
centerpoint.latitude <- median(cases$lat)
centerpoint.longitude.notLinked <- median(cases.notLinked$long)
centerpoint.latitude.notLinked <- median(cases.notLinked$lat)
centerpoint.longitude.linked <- median(cases.linked$long)
centerpoint.latitude.linked <- median(cases.linked$lat)
centerpoint.longitude.linked.B <- median(cases.linked.B$long)
centerpoint.latitude.linked.B <- median(cases.linked.B$lat)
centerpoint.longitude.notLinked.B <- median(cases.notLinked.B$long)
centerpoint.latitude.notLinked.B <- median(cases.notLinked.B$lat)
centerpoint.longitude.notLinked.A <- median(as.numeric(cases.notLinked.A$long))
centerpoint.latitude.notLinked.A <- median(as.numeric(cases.notLinked.A$lat))
centerpoint.longitude.notLinked.A1 <- median(as.numeric(cases.notLinked.A1$long))
centerpoint.latitude.notLinked.A1 <- median(as.numeric(cases.notLinked.A1$lat))
centerpoint.longitude.B <- median(cases.B$long)
centerpoint.latitude.B <- median(cases.B$lat)
centerpoint.longitude.confirmed <- median(cases.confirmed$long)
centerpoint.latitude.confirmed <- median(cases.confirmed$lat)
centerpoint.longitude.clinDiag <- median(cases.clinDiag$long)
centerpoint.latitude.clinDiag <- median(cases.clinDiag$lat)

centerpoint.dist <- (distVincentyEllipsoid(c(market.long,market.lat), c(centerpoint.longitude,centerpoint.latitude)))/1000
centerpoint.dist.notLinked <- (distVincentyEllipsoid(c(market.long,market.lat), c(centerpoint.longitude.notLinked,centerpoint.latitude.notLinked)))/1000
centerpoint.dist.linked <- (distVincentyEllipsoid(c(market.long,market.lat), c(centerpoint.longitude.linked,centerpoint.latitude.linked)))/1000
centerpoint.dist.linked.B <- (distVincentyEllipsoid(c(market.long,market.lat), c(centerpoint.longitude.linked.B,centerpoint.latitude.linked.B)))/1000
centerpoint.dist.notLinked.B <- (distVincentyEllipsoid(c(market.long,market.lat), c(centerpoint.longitude.notLinked.B,centerpoint.latitude.notLinked.B)))/1000
centerpoint.dist.notLinked.A <- (distVincentyEllipsoid(c(market.long,market.lat), c(centerpoint.longitude.notLinked.A,centerpoint.latitude.notLinked.A)))/1000
centerpoint.dist.notLinked.A1 <- (distVincentyEllipsoid(c(market.long,market.lat), c(centerpoint.longitude.notLinked.A1,centerpoint.latitude.notLinked.A1)))/1000
centerpoint.dist.B <- (distVincentyEllipsoid(c(market.long,market.lat), c(centerpoint.longitude.B,centerpoint.latitude.B)))/1000
centerpoint.dist.confirmed <- (distVincentyEllipsoid(c(market.long,market.lat), c(centerpoint.longitude.confirmed,centerpoint.latitude.confirmed)))/1000
centerpoint.dist.clinDiag <- (distVincentyEllipsoid(c(market.long,market.lat), c(centerpoint.longitude.clinDiag,centerpoint.latitude.clinDiag)))/1000

centerpoint.dist.all.pValue.weibo <- getPvalue(centerpoint.dist,as.numeric(centerpoint.distances.weibo$distance_market_km))
centerpoint.dist.notLinked.pValue.weibo <- getPvalue(centerpoint.dist.notLinked,as.numeric(centerpoint.distances.weibo$distance_market_km))
centerpoint.dist.linked.pValue.weibo <- getPvalue(centerpoint.dist.linked,as.numeric(centerpoint.distances.weibo$distance_market_km))
centerpoint.dist.linked.B.pValue.weibo <- getPvalue(centerpoint.dist.linked.B,as.numeric(centerpoint.distances.weibo$distance_market_km))
centerpoint.dist.notLinked.B.pValue.weibo <- getPvalue(centerpoint.dist.notLinked.B,as.numeric(centerpoint.distances.weibo$distance_market_km))
centerpoint.dist.notLinked.A.pValue.weibo <- getPvalue(centerpoint.dist.notLinked.A,as.numeric(centerpoint.distances.weibo$distance_market_km))
centerpoint.dist.notLinked.A1.pValue.weibo <- getPvalue(centerpoint.dist.notLinked.A1,as.numeric(centerpoint.distances.weibo$distance_market_km))
centerpoint.dist.confirmed.pValue.weibo <- getPvalue(centerpoint.dist.confirmed,as.numeric(centerpoint.distances.weibo$distance_market_km))
centerpoint.dist.clinDiag.pValue.weibo <- getPvalue(centerpoint.dist.clinDiag,as.numeric(centerpoint.distances.weibo$distance_market_km))

centerpoint.dist.all.pValue.worldpop.age <- getPvalue(centerpoint.dist,centerpoint.distances.worldpop.age$distance_huanan)
centerpoint.dist.notLinked.pValue.worldpop.age <- getPvalue(centerpoint.dist.notLinked,centerpoint.distances.worldpop.age$distance_huanan)
centerpoint.dist.linked.pValue.worldpop.age <- getPvalue(centerpoint.dist.linked,centerpoint.distances.worldpop.age$distance_huanan)
centerpoint.dist.linked.B.pValue.worldpop.age <- getPvalue(centerpoint.dist.linked.B,centerpoint.distances.worldpop.age$distance_huanan)
centerpoint.dist.notLinked.B.pValue.worldpop.age <- getPvalue(centerpoint.dist.notLinked.B,centerpoint.distances.worldpop.age$distance_huanan)
centerpoint.dist.notLinked.A.pValue.worldpop.age <- getPvalue(centerpoint.dist.notLinked.A,centerpoint.distances.worldpop.age$distance_huanan)
centerpoint.dist.notLinked.A1.pValue.worldpop.age <- getPvalue(centerpoint.dist.notLinked.A1,centerpoint.distances.worldpop.age$distance_huanan)
centerpoint.dist.B.pValue.worldpop.age <- getPvalue(centerpoint.dist.B,centerpoint.distances.worldpop.age$distance_huanan)
centerpoint.dist.confirmed.pValue.worldpop.age <- getPvalue(centerpoint.dist.confirmed,centerpoint.distances.worldpop.age$distance_huanan)
centerpoint.dist.clinDiag.pValue.worldpop.age <- getPvalue(centerpoint.dist.clinDiag,centerpoint.distances.worldpop.age$distance_huanan)

######## write p-values for median and centroid distance metrics to file(all cases, not linked, linked) ######## 
sink('analysis-output.csv')
cat("median distance, all cases, ", median(distances.toMarket), ", weibo null, p = ", median.dist.all.pValue.weibo, "\n")
cat("median distance, not linked, ", median(distances.toMarket.notLinked), ", weibo null, p = ", median.dist.notLinked.pValue.weibo, "\n")
cat("median distance, linked, ", median(distances.toMarket.linked), ",  weibo null, p = ", median.dist.linked.pValue.weibo, "\n")
cat("median distance, linked B, ", median(distances.toMarket.linked.B), ",  weibo null, p = ", median.dist.linked.B.pValue.weibo, "\n")
cat("median distance, B, ", median(distances.toMarket.B), ",  weibo null, p = ", median.dist.B.pValue.weibo, "\n")
cat("median distance, not linked B, ", median(distances.toMarket.notLinked.B), ",  weibo null, p = ", median.dist.notLinked.B.pValue.weibo, "\n")
cat("median distance, not linked A (1 + hotel guest), ", median(distances.toMarket.notLinked.A), ",  weibo null, p = ", median.dist.notLinked.A.pValue.weibo, "\n")
cat("median distance, not linked A (1), ", median(distances.toMarket.notLinked.A1), ",  weibo null, p = ", median.dist.notLinked.A1.pValue.weibo, "\n")
cat("median distance, confirmed, ", median(distances.toMarket.confirmed), ",  weibo null, p = ", median.dist.confirmed.pValue.weibo, "\n")
cat("median distance, clinically Diagnosed, ", median(distances.toMarket.clinDiag), ",  weibo null, p = ", median.dist.clinDiag.pValue.weibo, "\n")
cat("\n")
cat("median distance, all cases, ", median(distances.toMarket), ", worldpop age null, p = ", median.dist.all.pValue.worldpop.age, "\n")
cat("median distance, not linked, ", median(distances.toMarket.notLinked), ", worldpop age null, p = ", median.dist.notLinked.pValue.worldpop.age, "\n")
cat("median distance, linked, ", median(distances.toMarket.linked), ",  worldpop age null, p = ", median.dist.linked.pValue.worldpop.age, "\n")
cat("median distance, linked B, ", median(distances.toMarket.linked.B), ",  worldpop age null, p = ", median.dist.linked.B.pValue.worldpop.age, "\n")
cat("median distance, B, ", median(distances.toMarket.B), ",  worldpop age null, p = ", median.dist.B.pValue.worldpop.age, "\n")
cat("median distance, not linked B, ", median(distances.toMarket.notLinked.B), ",  worldpop age null, p = ", median.dist.notLinked.B.pValue.worldpop.age, "\n")
cat("median distance, not linked A (1 + hotel guest), ", median(distances.toMarket.notLinked.A), ",  worldpop age null, p = ", median.dist.notLinked.A.pValue.worldpop.age, "\n")
cat("median distance, not linked A (1), ", median(distances.toMarket.notLinked.A1), ",  worldpop age null, p = ", median.dist.notLinked.A1.pValue.worldpop.age, "\n")
cat("median distance, confirmed, ", median(distances.toMarket.confirmed), ",  worldpop age null, p = ", median.dist.confirmed.pValue.worldpop.age, "\n")
cat("median distance, clinically diagnosed, ", median(distances.toMarket.clinDiag), ",  worldpop age null, p = ", median.dist.clinDiag.pValue.worldpop.age, "\n")
cat("\n")
cat("centerpoint distance, all cases, ",centerpoint.dist, ", weibo null, p = ", centerpoint.dist.all.pValue.weibo, "\n")
cat("centerpoint distance, not linked, ",centerpoint.dist.notLinked, ", weibo null, p = ", centerpoint.dist.notLinked.pValue.weibo, "\n")
cat("centerpoint distance, linked, ",centerpoint.dist.linked, ", weibo null, p = ", centerpoint.dist.linked.pValue.weibo, "\n")
cat("centerpoint distance, linked B, ",centerpoint.dist.linked.B, ", weibo null, p = ", centerpoint.dist.linked.B.pValue.weibo, "\n")
cat("centerpoint distance, B, ",centerpoint.dist.B, ", weibo null, p = ", centerpoint.dist.B.pValue.weibo, "\n")
cat("centerpoint distance, not linked B, ",centerpoint.dist.notLinked.B, ", weibo null, p = ", centerpoint.dist.notLinked.B.pValue.weibo, "\n")
#cat("centerpoint distance, not linked A (1 + hotel guest), ",centerpoint.dist.notLinked.A, ", weibo null, p = ", centerpoint.dist.notLinked.A.pValue.weibo, "\n")
cat("centerpoint distance, not linked A (1), ",centerpoint.dist.notLinked.A1, ", weibo null, p = ", centerpoint.dist.notLinked.A1.pValue.weibo, "\n")
cat("centerpoint distance, confirmed, ",centerpoint.dist.confirmed, ", weibo null, p = ", centerpoint.dist.confirmed.pValue.weibo, "\n")
cat("centerpoint distance, clinically diagnosed, ",centerpoint.dist.clinDiag, ", weibo null, p = ", centerpoint.dist.clinDiag.pValue.weibo, "\n")
cat("\n")
cat("centerpoint distance, all cases, ",centerpoint.dist, ", worldpop age null, p = ", centerpoint.dist.all.pValue.worldpop.age, "\n")
cat("centerpoint distance, not linked, ",centerpoint.dist.notLinked, ", worldpop age null, p = ", centerpoint.dist.notLinked.pValue.worldpop.age, "\n")
cat("centerpoint distance, linked, ",centerpoint.dist.linked, ", worldpop age null, p = ", centerpoint.dist.linked.pValue.worldpop.age, "\n")
cat("centerpoint distance, linked B, ",centerpoint.dist.linked.B, ", worldpop age null, p = ", centerpoint.dist.linked.B.pValue.worldpop.age, "\n")
cat("centerpoint distance, B, ",centerpoint.dist.B, ", worldpop age null, p = ", centerpoint.dist.B.pValue.worldpop.age, "\n")
cat("centerpoint distance, not linked B, ",centerpoint.dist.notLinked.B, ", worldpop age null, p = ", centerpoint.dist.notLinked.B.pValue.worldpop.age, "\n")
#cat("centerpoint distance, not linked A (1 + hotel guest), ",centerpoint.dist.notLinked.A, ", worldpop age null, p = ", centerpoint.dist.notLinked.A.pValue.worldpop.age, "\n")
cat("centerpoint distance, not linked A (1), ",centerpoint.dist.notLinked.A1, ", worldpop age null, p = ", centerpoint.dist.notLinked.A1.pValue.worldpop.age, "\n")
cat("centerpoint distance, confirmed, ",centerpoint.dist.confirmed, ", worldpop age null, p = ", centerpoint.dist.confirmed.pValue.worldpop.age, "\n")
cat("centerpoint distance, clinically diagnosed, ",centerpoint.dist.clinDiag, ", worldpop age null, p = ", centerpoint.dist.clinDiag.pValue.worldpop.age, "\n")
sink()


######## robustness wrt location noise, unsampled cases (ghosts) and mislabeling ######## 
# noise stats
median.distances.noise <- vector(mode="numeric", length=numberOfReps)
centerpoint.distances.noise <- vector(mode="numeric", length=numberOfReps)
median.distances.notLinked.noise <- vector(mode="numeric", length=numberOfReps)
centerpoint.distances.notLinked.noise <- vector(mode="numeric", length=numberOfReps)
median.distances.linked.noise <- vector(mode="numeric", length=numberOfReps)
centerpoint.distances.linked.noise <- vector(mode="numeric", length=numberOfReps)
# noise & ghosts stats
median.distances.noise.ghosts <- vector(mode="numeric", length=numberOfReps)
centerpoint.distances.noise.ghosts <- vector(mode="numeric", length=numberOfReps)
median.distances.notLinked.noise.ghosts <- vector(mode="numeric", length=numberOfReps)
centerpoint.distances.notLinked.noise.ghosts <- vector(mode="numeric", length=numberOfReps)
median.distances.linked.noise.ghosts <- vector(mode="numeric", length=numberOfReps)
centerpoint.distances.linked.noise.ghosts <- vector(mode="numeric", length=numberOfReps)
# subsampled not linked stats
median.distances.notLinked.subsampled <- vector(mode="numeric", length=numberOfReps)
centerpoint.distances.notLinked.subsampled <- vector(mode="numeric", length=numberOfReps)

for (j in 1:numberOfReps){
  print(j)
  ##### re-sampling = 'noise' #####
  distances.toMarket.withNoise <- vector(mode="numeric", length=length(cases[,1]))
  distances.toMarket.notLinked.withNoise <- vector(mode="numeric", length=length(cases.notLinked[,1]))
  distances.toMarket.linked.withNoise <- vector(mode="numeric", length=length(cases.linked[,1]))
  cases.withNoise.long <- vector(mode="numeric", length=length(cases[,1]))
  cases.withNoise.lat <- vector(mode="numeric", length=length(cases[,1]))
  cases.withNoise.notLinked.long <- vector(mode="numeric", length=length(cases.notLinked[,1]))
  cases.withNoise.notLinked.lat <- vector(mode="numeric", length=length(cases.notLinked[,1]))
  cases.withNoise.linked.long <- vector(mode="numeric", length=length(cases.linked[,1]))
  cases.withNoise.linked.lat <- vector(mode="numeric", length=length(cases.linked[,1]))
  notLinked.count <- 1
  linked.count <- 1
  for (k in 1:length(cases$id)){
    circle <- create.circle.spPolygon(cases[k,2],cases[k,3],radius)
    projection(circle) <- CRS("+proj=longlat +datum=WGS84")
    df <- as.data.frame(cbind("circle around data point"))
    spdf <- SpatialPolygonsDataFrame(circle, df, match.ID = FALSE)
    spdf.sf = st_as_sf(spdf)                                                                                                                 
    pointWithNoise = sf::st_sample(spdf.sf, size=1)
    loc <-st_coordinates(pointWithNoise)
    cases.withNoise.long[k] <- loc[1]
    cases.withNoise.lat[k] <- loc[2]
    dist.withNoise<- (distVincentyEllipsoid(c(market.long,market.lat), c(loc[1],loc[2])))/1000
    distances.toMarket.withNoise[k] <- dist.withNoise
    if (cases$huanan_linked[k] == "FALSE"){
      cases.withNoise.notLinked.long[notLinked.count] <- loc[1]
      cases.withNoise.notLinked.lat[notLinked.count] <- loc[2]
      distances.toMarket.notLinked.withNoise[notLinked.count] <- dist.withNoise
      notLinked.count <- notLinked.count +1 
    } else {
      cases.withNoise.linked.long[linked.count] <- loc[1]
      cases.withNoise.linked.lat[linked.count] <- loc[2]
      distances.toMarket.linked.withNoise[linked.count] <- dist.withNoise
      linked.count <- linked.count +1 
    }
  }
  
  median.distances.noise[j] <- median(distances.toMarket.withNoise)
  median.distances.notLinked.noise[j] <- median(distances.toMarket.notLinked.withNoise)
  median.distances.linked.noise[j] <- median(distances.toMarket.linked.withNoise)
  centerpoint.dist.noise <- (distVincentyEllipsoid(c(market.long,market.lat), c(median(cases.withNoise.long),median(cases.withNoise.lat))))/1000
  centerpoint.distances.noise[j] <- centerpoint.dist.noise
  centerpoint.dist.notLinked.noise <- (distVincentyEllipsoid(c(market.long,market.lat), c(median(cases.withNoise.notLinked.long),median(cases.withNoise.notLinked.lat))))/1000
  centerpoint.distances.notLinked.noise[j] <- centerpoint.dist.notLinked.noise
  centerpoint.dist.linked.noise <- (distVincentyEllipsoid(c(market.long,market.lat), c(median(cases.withNoise.linked.long),median(cases.withNoise.linked.lat))))/1000
  centerpoint.distances.linked.noise[j] <- centerpoint.dist.linked.noise
  
  ##### subsample not linked cases, always including 66 and 153 ######
  distances.toMarket.notLinked.subsampled <- vector(mode="numeric", length=numberOfSubsampled)
  cases.subsampled.notLinked.long <- vector(mode="numeric", length=numberOfSubsampled)
  cases.subsampled.notLinked.lat <- vector(mode="numeric", length=numberOfSubsampled)
  cases.notLinked.2fixed <- cases.notLinked[cases.notLinked$id == "66" | cases.notLinked$id == "153",]
  cases.notLinked.except.2fixed <- cases.notLinked[cases.notLinked$id != "66" & cases.notLinked$id != "153",]
  cases.subsampled.notLinked.long[1] <- cases.notLinked.2fixed$long[1] 
  cases.subsampled.notLinked.lat[1] <- cases.notLinked.2fixed$lat[1] 
  cases.subsampled.notLinked.long[2] <- cases.notLinked.2fixed$long[2] 
  cases.subsampled.notLinked.lat[2] <- cases.notLinked.2fixed$lat[2] 
  distances.toMarket.notLinked.subsampled[1] <-  (distVincentyEllipsoid(c(market.long,market.lat), c(cases.subsampled.notLinked.long[1],cases.subsampled.notLinked.lat[1])))/1000
  distances.toMarket.notLinked.subsampled[2] <-  (distVincentyEllipsoid(c(market.long,market.lat), c(cases.subsampled.notLinked.long[2],cases.subsampled.notLinked.lat[2])))/1000
  indices.rand <- sample(1:length(cases.notLinked.except.2fixed$id))

  for (a in 3:numberOfSubsampled){
    cases.subsampled.notLinked.long[a] <- cases.notLinked.except.2fixed$long[indices.rand[(a-2)]]
    cases.subsampled.notLinked.lat[a] <- cases.notLinked.except.2fixed$lat[indices.rand[(a-2)]]
    distances.toMarket.notLinked.subsampled[a] <- (distVincentyEllipsoid(c(market.long,market.lat), c(cases.subsampled.notLinked.long[a],cases.subsampled.notLinked.lat[a])))/1000
  }

  median.distances.notLinked.subsampled[j] <- median(distances.toMarket.notLinked.subsampled)
  centerpoint.dist.notLinked.subsampled <- (distVincentyEllipsoid(c(market.long,market.lat), c(median(cases.subsampled.notLinked.long),median(cases.subsampled.notLinked.lat))))/1000
  centerpoint.distances.notLinked.subsampled[j] <- centerpoint.dist.notLinked.subsampled

  ##### adding 'unsampled cases', a.k.a. ghosts #####
  distances.toMarket.withNoise.withGhosts <- distances.toMarket.withNoise
  distances.toMarket.notLinked.withNoise.withGhosts <- distances.toMarket.notLinked.withNoise
  distances.toMarket.linked.withNoise.withGhosts <- distances.toMarket.linked.withNoise
  cases.withNoise.withGhosts.long <- cases.withNoise.long
  cases.withNoise.withGhosts.lat <- cases.withNoise.lat
  cases.withNoise.notLinked.withGhosts.long <- cases.withNoise.notLinked.long
  cases.withNoise.notLinked.withGhosts.lat <- cases.withNoise.notLinked.lat
  cases.withNoise.linked.withGhosts.long <- cases.withNoise.linked.long
  cases.withNoise.linked.withGhosts.lat <- cases.withNoise.linked.lat
  
  market.polygons<-getPolygons(cbind(cases$longitude,cases$latitude),unsampled_KDEprobabily) 
  market.spatialPolygons = SpatialPolygons(list(Polygons(market.polygons,1)))
  market.spdf = SpatialPolygonsDataFrame(market.spatialPolygons, data.frame(ID=1:length(market.spatialPolygons)))
  market.spdf.sf = st_as_sf(market.spdf)                                                                                                                 
  ghostpoints = sf::st_sample(market.spdf.sf, size=numberOfGhosts)
  for (l in 1:numberOfGhosts){
    loc.ghost <-st_coordinates(ghostpoints[l])
    dist.ghost<- (distVincentyEllipsoid(c(market.long,market.lat), c(loc.ghost[1],loc.ghost[2])))/1000
    # all ghosts are considered linked, so added to linked only
    distances.toMarket.withNoise.withGhosts <- c(distances.toMarket.withNoise.withGhosts, dist.ghost)
    distances.toMarket.linked.withNoise.withGhosts <- c(distances.toMarket.withNoise.withGhosts, dist.ghost)
    cases.withNoise.withGhosts.long <- c(cases.withNoise.withGhosts.long, loc.ghost[1])
    cases.withNoise.withGhosts.lat <- c(cases.withNoise.withGhosts.lat, loc.ghost[2])
    cases.withNoise.linked.withGhosts.long <- c(cases.withNoise.linked.withGhosts.long, loc.ghost[1])
    cases.withNoise.linked.withGhosts.lat <- c(cases.withNoise.linked.withGhosts.lat , loc.ghost[2])
  }
  
  median.distances.noise.ghosts[j] <- median(distances.toMarket.withNoise.withGhosts)
  median.distances.notLinked.noise.ghosts[j] <- median(distances.toMarket.notLinked.withNoise.withGhosts)
  median.distances.linked.noise.ghosts[j] <- median(distances.toMarket.linked.withNoise.withGhosts)
  centerpoint.dist.noise.ghosts <- (distVincentyEllipsoid(c(market.long,market.lat), c(median(cases.withNoise.withGhosts.long),median(cases.withNoise.withGhosts.lat))))/1000
  centerpoint.distances.noise.ghosts[j] <- centerpoint.dist.noise.ghosts
  centerpoint.dist.notLinked.noise.ghosts <- (distVincentyEllipsoid(c(market.long,market.lat), c(median(cases.withNoise.notLinked.withGhosts.long),median(cases.withNoise.notLinked.withGhosts.lat))))/1000
  centerpoint.distances.notLinked.noise.ghosts[j] <- centerpoint.dist.notLinked.noise.ghosts
  centerpoint.dist.linked.noise.ghosts <- (distVincentyEllipsoid(c(market.long,market.lat), c(median(cases.withNoise.linked.withGhosts.long),median(cases.withNoise.linked.withGhosts.lat))))/1000
  centerpoint.distances.linked.noise.ghosts[j] <- centerpoint.dist.linked.noise.ghosts
}


#### PLOTTING ######
a <- data.frame(group="median", value=median.distances.noise, cases = "all")
b <- data.frame(group="median", value=median.distances.notLinked.noise, cases = "not linked")
c <- data.frame(group="median", value=median.distances.linked.noise, cases = "linked")
noise.median <- rbind(a,b,c)
d <- data.frame(group="center-point", value=centerpoint.distances.noise, cases = "all")
e <- data.frame(group="center-point", value=centerpoint.distances.notLinked.noise, cases = "not linked")
f <- data.frame(group="center-point", value=centerpoint.distances.linked.noise, cases = "linked")
noise.centerpoint <- rbind(d,e,f)

noise  <- rbind(noise.median,noise.centerpoint)

stat_order <- c("median", "center-point")  
pdf("noise.pdf",
    width = 10, height = 6)
p1<-ggplot(noise, aes(x=factor(group, stat_order), y=value, fill=cases)) +
  geom_boxplot(outlier.colour="black", outlier.shape=16,
               outlier.size=2, notch=FALSE) +
  geom_segment(aes(x=0.64, xend=0.85, y=quantile(as.numeric(median.distances.all.worldpop.age$Median_distance_huanan), probs = c(0.05)), yend=quantile(as.numeric(median.distances.all.worldpop.age$Median_distance_huanan), probs = c(0.05))), colour="black", linetype="dotted", size=0.85) +
  geom_segment(aes(x=0.9, xend=1.1, y=quantile(as.numeric(median.distances.linked.worldpop.age$Median_distance_huanan), probs = c(0.05)), yend=quantile(as.numeric(median.distances.linked.worldpop.age$Median_distance_huanan), probs = c(0.05))), colour="black", linetype="dotted", size=0.85) +
  geom_segment(aes(x=1.14, xend=1.35, y=quantile(as.numeric(median.distances.notLinked.worldpop.age$Median_distance_huanan), probs = c(0.05)), yend=quantile(as.numeric(median.distances.notLinked.worldpop.age$Median_distance_huanan), probs = c(0.05))), colour="black", linetype="dotted", size=0.85) +

  geom_segment(aes(x=0.64, xend=0.85, y=quantile(as.numeric(median.distances.all.worldpop.age$Median_distance_huanan), probs = c(0.1)), yend=quantile(as.numeric(median.distances.all.worldpop.age$Median_distance_huanan), probs = c(0.1))), colour="black", size=0.85) +
  geom_segment(aes(x=0.9, xend=1.1, y=quantile(as.numeric(median.distances.linked.worldpop.age$Median_distance_huanan), probs = c(0.1)), yend=quantile(as.numeric(median.distances.linked.worldpop.age$Median_distance_huanan), probs = c(0.1))), colour="black", size=0.85) +
  geom_segment(aes(x=1.14, xend=1.35, y=quantile(as.numeric(median.distances.notLinked.worldpop.age$Median_distance_huanan), probs = c(0.1)), yend=quantile(as.numeric(median.distances.notLinked.worldpop.age$Median_distance_huanan), probs = c(0.1))), colour="black", size=0.85) +
  
  geom_segment(aes(x=0.64, xend=0.85, y=quantile(as.numeric(median.distances.all.weibo$Median_distance_huanan), probs = c(0.05)), yend=quantile(as.numeric(median.distances.all.weibo$Median_distance_huanan), probs = c(0.05))), colour="grey", linetype="dotted", size=0.85) +
  geom_segment(aes(x=0.9, xend=1.1, y=quantile(as.numeric(median.distances.linked.weibo$Median_distance_huanan), probs = c(0.05)), yend=quantile(as.numeric(median.distances.linked.weibo$Median_distance_huanan), probs = c(0.05))), colour="grey", linetype="dotted", size=0.85) +
  geom_segment(aes(x=1.14, xend=1.35, y=quantile(as.numeric(median.distances.notLinked.weibo$Median_distance_huanan), probs = c(0.05)), yend=quantile(as.numeric(median.distances.notLinked.weibo$Median_distance_huanan), probs = c(0.05))), colour="grey", linetype="dotted", size=0.85) +
  
  geom_segment(aes(x=0.64, xend=0.85, y=quantile(as.numeric(median.distances.all.weibo$Median_distance_huanan), probs = c(0.1)), yend=quantile(as.numeric(median.distances.all.weibo$Median_distance_huanan), probs = c(0.1))), colour="grey", size=0.85) +
  geom_segment(aes(x=0.9, xend=1.1, y=quantile(as.numeric(median.distances.linked.weibo$Median_distance_huanan), probs = c(0.1)), yend=quantile(as.numeric(median.distances.linked.weibo$Median_distance_huanan), probs = c(0.1))), colour="grey", size=0.85) +
  geom_segment(aes(x=1.14, xend=1.35, y=quantile(as.numeric(median.distances.notLinked.weibo$Median_distance_huanan), probs = c(0.1)), yend=quantile(as.numeric(median.distances.notLinked.weibo$Median_distance_huanan), probs = c(0.1))), colour="grey", size=0.85) +

  geom_segment(aes(x=1.5, xend=2.5, y=quantile(as.numeric(centerpoint.distances.weibo[,1]), probs = c(0.05)), yend=quantile(as.numeric(centerpoint.distances.weibo[,1]), probs = c(0.05))), colour="grey",linetype="dotted", size=0.85) +
  geom_segment(aes(x=1.5, xend=2.5, y=quantile(as.numeric(centerpoint.distances.weibo[,1]), probs = c(0.1)), yend=quantile(as.numeric(centerpoint.distances.weibo[,1]), probs = c(0.1))), colour="grey", size=0.85) +
  geom_segment(aes(x=1.5, xend=2.5, y=quantile(as.numeric(centerpoint.distances.worldpop.age$distance_huanan), probs = c(0.05)), yend=quantile(as.numeric(centerpoint.distances.worldpop.age$distance_huanan), probs = c(0.05))), colour="black", linetype="dotted", size=0.85) +
  geom_segment(aes(x=1.5, xend=2.5, y=quantile(as.numeric(centerpoint.distances.worldpop.age$distance_huanan), probs = c(0.1)), yend=quantile(as.numeric(centerpoint.distances.worldpop.age$distance_huanan), probs = c(0.1))), colour="black", size=0.85) +
  coord_cartesian(ylim=c(0,8)) +
  scale_fill_manual(values=c('#e66101','#fdb863','#b2abd2')) +
  xlab("") + ylab("distance (km)") +
  ggtitle("noise (radius =  1 km)") +
  theme(axis.text = element_text(size=12)) +
  theme(axis.title.y = element_text(size = 12)) +
  theme(plot.title = element_text(size = 14, hjust = 0.5))
p1
dev.off() 


a <- data.frame(group="median", value=median.distances.noise.ghosts, cases = "all")
b <- data.frame(group="median", value=median.distances.notLinked.noise.ghosts, cases = "not linked")
c <- data.frame(group="median", value=median.distances.linked.noise.ghosts, cases = "linked")
noise.ghosts.median <- rbind(a,b,c)
d <- data.frame(group="center-point", value=centerpoint.distances.noise.ghosts, cases = "all")
e <- data.frame(group="center-point", value=centerpoint.distances.notLinked.noise.ghosts, cases = "not linked")
f <- data.frame(group="center-point", value=centerpoint.distances.linked.noise.ghosts, cases = "linked")
noise.ghosts.centerpoint <- rbind(d,e,f)
noise.ghosts  <- rbind(noise.ghosts.median,noise.ghosts.centerpoint)

stat_order <- c("median", "center-point")  
pdf("noise_missing.pdf",
    width = 10, height = 6)
p2<-ggplot(noise.ghosts, aes(x=factor(group, stat_order), y=value, fill=cases)) +
  geom_boxplot(outlier.colour="black", outlier.shape=16,
               outlier.size=2, notch=FALSE) +
  geom_segment(aes(x=0.64, xend=0.85, y=quantile(as.numeric(median.distances.all.ghosts.worldpop.age$Median_distance_huanan), probs = c(0.05)), yend=quantile(as.numeric(median.distances.all.ghosts.worldpop.age$Median_distance_huanan), probs = c(0.05))), colour="black", linetype="dotted", size=0.85) +
  geom_segment(aes(x=0.9, xend=1.1, y=quantile(as.numeric(median.distances.linked.ghosts.worldpop.age$Median_distance_huanan), probs = c(0.05)), yend=quantile(as.numeric(median.distances.linked.ghosts.worldpop.age$Median_distance_huanan), probs = c(0.05))), colour="black", linetype="dotted", size=0.85) +
  geom_segment(aes(x=1.14, xend=1.35, y=quantile(as.numeric(median.distances.notLinked.ghosts.worldpop.age$Median_distance_huanan), probs = c(0.05)), yend=quantile(as.numeric(median.distances.notLinked.ghosts.worldpop.age$Median_distance_huanan), probs = c(0.05))), colour="black", linetype="dotted", size=0.85) +
  
  geom_segment(aes(x=0.64, xend=0.85, y=quantile(as.numeric(median.distances.all.ghosts.worldpop.age$Median_distance_huanan), probs = c(0.1)), yend=quantile(as.numeric(median.distances.all.ghosts.worldpop.age$Median_distance_huanan), probs = c(0.1))), colour="black", size=0.85) +
  geom_segment(aes(x=0.9, xend=1.1, y=quantile(as.numeric(median.distances.linked.ghosts.worldpop.age$Median_distance_huanan), probs = c(0.1)), yend=quantile(as.numeric(median.distances.linked.ghosts.worldpop.age$Median_distance_huanan), probs = c(0.1))), colour="black", size=0.85) +
  geom_segment(aes(x=1.14, xend=1.35, y=quantile(as.numeric(median.distances.notLinked.ghosts.worldpop.age$Median_distance_huanan), probs = c(0.1)), yend=quantile(as.numeric(median.distances.notLinked.ghosts.worldpop.age$Median_distance_huanan), probs = c(0.1))), colour="black", size=0.85) +
  
  geom_segment(aes(x=0.64, xend=0.85, y=quantile(as.numeric(median.distances.all.ghosts.weibo$Median_distance_huanan), probs = c(0.05)), yend=quantile(as.numeric(median.distances.all.ghosts.weibo$Median_distance_huanan), probs = c(0.05))), colour="grey", linetype="dotted", size=0.85) +
  geom_segment(aes(x=0.9, xend=1.1, y=quantile(as.numeric(median.distances.linked.ghosts.weibo$Median_distance_huanan), probs = c(0.05)), yend=quantile(as.numeric(median.distances.linked.ghosts.weibo$Median_distance_huanan), probs = c(0.05))), colour="grey", linetype="dotted", size=0.85) +
  geom_segment(aes(x=1.14, xend=1.35, y=quantile(as.numeric(median.distances.notLinked.ghosts.weibo$Median_distance_huanan), probs = c(0.05)), yend=quantile(as.numeric(median.distances.notLinked.ghosts.weibo$Median_distance_huanan), probs = c(0.05))), colour="grey", linetype="dotted", size=0.85) +
  
  geom_segment(aes(x=0.64, xend=0.85, y=quantile(as.numeric(median.distances.all.ghosts.weibo$Median_distance_huanan), probs = c(0.1)), yend=quantile(as.numeric(median.distances.all.ghosts.weibo$Median_distance_huanan), probs = c(0.1))), colour="grey", size=0.85) +
  geom_segment(aes(x=0.9, xend=1.1, y=quantile(as.numeric(median.distances.linked.ghosts.weibo$Median_distance_huanan), probs = c(0.1)), yend=quantile(as.numeric(median.distances.linked.ghosts.weibo$Median_distance_huanan), probs = c(0.1))), colour="grey", size=0.85) +
  geom_segment(aes(x=1.14, xend=1.35, y=quantile(as.numeric(median.distances.notLinked.ghosts.weibo$Median_distance_huanan), probs = c(0.1)), yend=quantile(as.numeric(median.distances.notLinked.ghosts.weibo$Median_distance_huanan), probs = c(0.1))), colour="grey", size=0.85) +
  
  geom_segment(aes(x=1.5, xend=2.5, y=quantile(as.numeric(centerpoint.distances.weibo[,1]), probs = c(0.05)), yend=quantile(as.numeric(centerpoint.distances.weibo[,1]), probs = c(0.05))), colour="grey",linetype="dotted", size=0.85) +
  geom_segment(aes(x=1.5, xend=2.5, y=quantile(as.numeric(centerpoint.distances.weibo[,1]), probs = c(0.1)), yend=quantile(as.numeric(centerpoint.distances.weibo[,1]), probs = c(0.1))), colour="grey", size=0.85) +
  geom_segment(aes(x=1.5, xend=2.5, y=quantile(as.numeric(centerpoint.distances.worldpop.age$distance_huanan), probs = c(0.05)), yend=quantile(as.numeric(centerpoint.distances.worldpop.age$distance_huanan), probs = c(0.05))), colour="black", linetype="dotted", size=0.85) +
  geom_segment(aes(x=1.5, xend=2.5, y=quantile(as.numeric(centerpoint.distances.worldpop.age$distance_huanan), probs = c(0.1)), yend=quantile(as.numeric(centerpoint.distances.worldpop.age$distance_huanan), probs = c(0.1))), colour="black", size=0.85) +
  coord_cartesian(ylim=c(0,8)) +
  scale_fill_manual(values=c('#e66101','#fdb863','#b2abd2')) +
  xlab("") + ylab("distance (km)") +
  ggtitle("noise (radius =  1 km), missing (n=9)") +
  theme(axis.text = element_text(size=12)) +
  theme(axis.title.y = element_text(size = 12)) +
  theme(plot.title = element_text(size = 14, hjust = 0.5))
p2
dev.off() 

a <- data.frame(group="median", value=median.distances.notLinked.subsampled, cases = "not linked")
b <- data.frame(group="center-point", value=centerpoint.distances.notLinked.subsampled, cases = "not linked")
notLinked.unsampled<-rbind(a,b)

stat_order <- c("median", "center-point")  
pdf("notLinked_subsampled.pdf", 
    width = 4, height = 6)
p3<-ggplot(notLinked.unsampled, aes(x=factor(group, stat_order), y=value, fill=cases)) +
  geom_boxplot(outlier.colour="black", outlier.shape=16,
               outlier.size=2, notch=FALSE) +
  geom_segment(aes(x=0.5, xend=1.5, y=quantile(as.numeric(median.distances.notLinked.sub.worldpop.age$Median_distance_huanan), probs = c(0.05)), yend=quantile(as.numeric(median.distances.notLinked.sub.worldpop.age$Median_distance_huanan), probs = c(0.05))), colour="black", linetype="dotted", size=0.85) +
  geom_segment(aes(x=0.5, xend=1.5, y=quantile(as.numeric(median.distances.notLinked.sub.worldpop.age$Median_distance_huanan), probs = c(0.1)), yend=quantile(as.numeric(median.distances.notLinked.sub.worldpop.age$Median_distance_huanan), probs = c(0.1))), colour="black", size=0.85) +
  geom_segment(aes(x=0.5, xend=1.5, y=quantile(as.numeric(median.distances.notLinked.sub.weibo$Median_distance_huanan), probs = c(0.05)), yend=quantile(as.numeric(median.distances.notLinked.sub.weibo$Median_distance_huanan), probs = c(0.05))), colour="grey", linetype="dotted", size=0.85) +
  geom_segment(aes(x=0.5, xend=1.5, y=quantile(as.numeric(median.distances.notLinked.sub.weibo$Median_distance_huanan), probs = c(0.1)), yend=quantile(as.numeric(median.distances.notLinked.sub.weibo$Median_distance_huanan), probs = c(0.1))), colour="grey", size=0.85) +
  
  geom_segment(aes(x=1.5, xend=2.5, y=quantile(as.numeric(centerpoint.distances.weibo[,1]), probs = c(0.05)), yend=quantile(as.numeric(centerpoint.distances.weibo[,1]), probs = c(0.05))), colour="grey",linetype="dotted", size=0.85) +
  geom_segment(aes(x=1.5, xend=2.5, y=quantile(as.numeric(centerpoint.distances.weibo[,1]), probs = c(0.1)), yend=quantile(as.numeric(centerpoint.distances.weibo[,1]), probs = c(0.1))), colour="grey", size=0.85) +
  geom_segment(aes(x=1.5, xend=2.5, y=quantile(as.numeric(centerpoint.distances.worldpop.age$distance_huanan), probs = c(0.05)), yend=quantile(as.numeric(centerpoint.distances.worldpop.age$distance_huanan), probs = c(0.05))), colour="black", linetype="dotted", size=0.85) +
  geom_segment(aes(x=1.5, xend=2.5, y=quantile(as.numeric(centerpoint.distances.worldpop.age$distance_huanan), probs = c(0.1)), yend=quantile(as.numeric(centerpoint.distances.worldpop.age$distance_huanan), probs = c(0.1))), colour="black", size=0.85) +
  coord_cartesian(ylim=c(0,5)) +
  scale_fill_manual(values=c('#b2abd2')) +
  xlab("") + ylab("distance (km)") +
  ggtitle("not linked, subsampled") +
  theme(axis.text = element_text(size=12)) +
  theme(axis.title.y = element_text(size = 12)) +
  theme(plot.title = element_text(size = 14, hjust = 0.5)) +
  theme(legend.position = "none")
p3
dev.off()

