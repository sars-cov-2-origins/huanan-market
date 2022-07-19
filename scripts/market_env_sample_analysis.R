# https://doi.org/10.5281/zenodo.6786454


# load recommended packages for sparr relative risk analysis
library(readxl)
library(janitor)
library(sf)
library(ggplot2)
library(dplyr)
library(spatstat)
library(sparr)
#color palette tools
library(RColorBrewer)
library(pals)
#sets working directory to file location (if using rstudio)
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
# otherwise use this one
# setwd(getSrcDirectory()[1])

#load geojson
gdf = st_read("../maps/geojson/huanan-market-internal.geojson")
market_boundary = gdf[gdf$group=='Boundary',]

# for boundary-respecting relative risk analysis, we project
# from EPSG: 4326 (lat/lon) into EPSG: 3857
# Also known as the pseudo-mercator projection (used in Google Maps, ArcGIS, etc.)
bound_t = st_transform(market_boundary, 3857)
window <- as.owin(bound_t)

## read in positive sample locations, project to 3857
pos_sites = gdf[gdf$group=='Env-Pos',]
pos_sites$title = lapply(pos_sites$title, FUN=function (x) strsplit(x," (",fixed=T)[[1]][1]) #strip info about sample type
pos_sites = st_transform(pos_sites, 3857)

## read in positive sample locations, project to 3857
neg_sites = gdf[gdf$group=='Env-Neg',]
neg_sites = st_transform(neg_sites, 3857)


### assemble a list of all sampled stalls
# first grab all stalls with positive samples
# Group samples in stalls with multiple positive samples
pos_sites2 = pos_sites
same_stall=list(list('A14','A15','A90','A88','A87'),list('Q61','Q64','Q69','Q70','Q68'),list('A2','A18','A20'),list('F98','F100'),list('Q93','A61'))
stall_center=list(c(114.25658,30.61937),c(114.256469, 30.61947),c(114.256724, 30.619613),c(114.25704,30.61949),c(114.25663,30.61966))
for (i in seq(1, to=length(same_stall), by=1)){
  pos_sites2 = filter(pos_sites2, !(title %in% same_stall[[i]]))
  pos_sites2 = rbind(pos_sites2,st_sf(title = paste("combo",i),group='Env-Pos',label=NA,marker.color='#cc1b15',stroke=NA,fill=NA, geometry=st_transform(st_sfc(st_point(stall_center[[i]]),crs=4326),3857)))
}

#pos_sites2 only contains the center of each positive stall (max 1 per business)
#concatenate with negative stalls
all_sites = rbind(pos_sites2,neg_sites)
all_coords <- matrix(unlist(all_sites$geometry), ncol = 2, byrow = T)

# convert list of sampled stalls to r(spatstat) ppp object format
control_ppp = ppp(x = all_coords[,1], y = all_coords[,2],
                  window = window, check = T)

#Build spatstat ppp object, list of coordinates of all positive samples
pos_coords <- matrix(unlist(pos_sites$geometry), ncol = 2, byrow = T)
#LSCV can be unstable when presented with multiple points with the exact same location
# (due to multiple positives in a single stall)
# so we perturb points by a small, random epsilon,
# orders of magnitude smaller than the variation in the data 
#Note: results are robust to the choice of eps (e.g. results consistent for eps = 1E-5,...,1E-8)
eps=1E-7
pos_ppp <- ppp(x = pos_coords[,1]+rnorm(dim(pos_coords)[1])*eps,
               y = pos_coords[,2]+rnorm(dim(pos_coords)[1])*eps,
               window = window, check = T)

#Build spatstat ppp object, list of coordinates of all negative stalls
neg_coords <- matrix(unlist(neg_sites$geometry), ncol = 2, byrow = T)
neg_ppp = ppp(x = neg_coords[,1], y = neg_coords[,2],
                  window = window, check = T)

#Build spatstat ppp object, list of coordinates of all positive stalls
pos2_coords <- matrix(unlist(pos_sites2$geometry), ncol = 2, byrow = T)
pos2_ppp = ppp(x = pos2_coords[,1], y = pos2_coords[,2],
                  window = window, check = T)

# perform relative risk analysis with sparr
# testing for increased risk, relative to null distribution of stalls sampled
cols <- brewer.pal(5, "Purples")
pal <- colorRampPalette(cols)
# bandwidth selection with least-squares cross-validation using fixed kernel. 
# h0 = LSCV.risk(pos_ppp,control_ppp,type="fixed",seqres=100)
h0 = 13.80062 # bandwidth selected by cross validation

# make relative risk plot
pdf('relative_risk_pos_vs_all.pdf')
## run relative risk analysis, calculate the p-value/tolerance surface
p1 <- risk(pos_ppp,control_ppp,resolution=1024,tolerate=T,log=F,adapt=F,h0=h0)
plot(p1,tol.type = c("upper"), tol.args = list(levels=c(0.05,0.01),lty=2:1, drawlabels = TRUE),col=pal(15))#,zlim=c(0,4.662257))
dev.off()
my.max <- function(x) ifelse( !all(is.na(x)), max(x, na.rm=T), NA)

# plot empirical null distribution, all stalls sampled in the market. 
# use bandwidth selected by LSCV
cols <- c("#FFFFFF", "#5f637d")
pal <- colorRampPalette(cols)
pdf('kde_all_stalls.pdf')
p1 <- bivariate.density(control_ppp,resolution=1024,h0=h0)
plot(p1,col=pal(15),add.pts=F)
plot(control_ppp,cex=0.5,add=T,col='black')
dev.off()

#assemble custom White-Yellow-Orange-Red Colormap
# plot distribution of positive environmental samples
# use bandwidth selected by LSCV
cols <- brewer.pal(5, "YlOrRd")
cols = c("#FFFFFF",cols)
pal <- colorRampPalette(cols)
pdf('kde_positives_sample_level.pdf')
p1 <- bivariate.density(pos_ppp,resolution=1024,h0=h0)
plot(p1,col=pal(15),add.pts=F)
plot(pos_ppp,cex=0.5,add=T,col='black')
dev.off()


