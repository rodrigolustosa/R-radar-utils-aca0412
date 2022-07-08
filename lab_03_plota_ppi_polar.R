# ----------------------------------------------------------------------- #
# Information: 
# Created by: Rodrigo Lustosa
# Creation date: 15 Jun 2022 16:16 (GMT -03)
# ----------------------------------------------------------------------- #

# packages
library(tidyverse)
library(fields)

# directory and file names
file_functions <- "functions.R"
dir_data <- "data"
file_nexrad <- "Nexrad/2010-03-13_2037.dat"


# functions ---------------------------------------------------------------

source(file_functions)


# radar information -------------------------------------------------------

# Nexrad
rmax      <- 100.     # distancia maxima
delta_azi <- 0.5      # amostragem em azimute
gatesize  <- 0.125    # resolucao do gate
x0        <- 0.       # posicao do radar em coordenadas cartesiana
y0        <- 0.       

lat_radar <- -(23 + 33/60 +  42.2/3600)
lon_radar <- -(46 + 44/60 +  06.4/3600)


#################
nbins = round(rmax/gatesize) #+ 1
nazim = round(360./delta_azi)#+ 1

# Miniradar
#rmax      = 21.6     # distancia maxima
#delta_azi = 1.0      # amostragem em azimute
#gatesize  = 0.090    # resolucao do gate
#x0        = 0.       # posicao do radar em coordenadas cartesiana
#y0        = 0.
#lat_radar = -(23 + 33/60 +  42.2/3600)
#lon_radar = -(46 + 44/60 +  06.4/3600)
# 
# #################
# nbins = round(rmax/gatesize) 
# nazim = round(360./delta_azi) #+ 1


# reading and computations ------------------------------------------------

# ########### Calculando as coordenadas polares
xy <- polar_2_cartesiana(nbins,nazim,delta_azi,gatesize,x0,y0)
x <- xy$x
y <- xy$y

# Calculando as coordenadas geograficas
xypolar <- wsr_polar(lon_radar,lat_radar,gatesize,nbins,nazim,delta_azi)
xpolar <- xypolar$xpolar
ypolar <- xypolar$ypolar

# Nexrad
file_radar <- file.path(dir_data,file_nexrad)
chuva <- le_nexrad_ascii(file_radar,nazim,nbins,delta_azi,gatesize)
az    <- seq(0,360-delta_azi,delta_azi)
dist  <- seq(gatesize/2,rmax,gatesize)

# # interpolate (nearest neighborhood) to regular grid (cartesian)
# chuva_cart <- interp_polar_2_cartesiana(chuva,delta_azi,gatesize,rmax)
# # chuva_cart <- interp_polar_2_cartesiana(ifelse(is.na(chuva),-99,chuva),
# #                                         delta_azi,gatesize,rmax)

# # interpolate (nearest neighborhood) to regular grid (lat/lon)
# chuva_geo <- interp_polar_2_geograph(chuva,delta_azi,gatesize,rmax,
#                                      c(lon_radar,lat_radar))
# # chuva_geo <- interp_polar_2_geograph(ifelse(is.na(chuva),-99,chuva),delta_azi,
# #                                      gatesize,rmax,c(lon_radar,lat_radar))



# plots -------------------------------------------------------------------


# plot data (cartesian)
imagePlot(x,y,chuva,zlim=c(0,60),#ifelse(is.na(chuva),-99,chuva),zlim=c(-99,60),
          main="Cartesianas",xlab="Oeste-Leste (km)",ylab="Sul-Norte (km)",
          legend.lab="Z(dBZ)")
abline(h=0,lty=2);abline(v=0,lty=2)

# plot data (lat/lon)
imagePlot(xpolar,ypolar,chuva,zlim=c(0,60),#zlim=c(-99,60),
          # xlim=c(-46.8,-46.7),ylim=c(-23.6,-23.5),
          main="Geográfica (lat/lon)",xlab="Longitude (°)",ylab="Latitude (°)",
          legend.lab="Z(dBZ)")
abline(h=lat_radar,lty=2);abline(v=lon_radar,lty=2)

# # plot interpolated data (cartesian)
# imagePlot(chuva_cart$x,chuva_cart$y,chuva_cart$z,zlim=c(0,60),#zlim=c(-99,60),
#            main="Cartesianas",xlab="Oeste-Leste (km)",ylab="Sul-Norte (km)",
#            legend.lab="Z(dBZ)")
# abline(h=0,lty=2);abline(v=0,lty=2)

# # plot interpolated data (lat/lon)
# imagePlot(chuva_geo$lon,chuva_geo$lat,chuva_geo$z,zlim=c(0,60),#zlim=c(-99,60),
#            # xlim=c(-46.8,-46.7),ylim=c(-23.6,-23.5),
#            main="Geográfica (lat/lon)",xlab="Longitude (°)",ylab="Latitude (°)",
#            legend.lab="Z(dBZ)")
# abline(h=lat_radar,lty=2);abline(v=lon_radar,lty=2)



# tests -------------------------------------------------------------------


# # combinations of (base) contour and filled contour
# contour(chuva)
# filled.contour(chuva)
# # filled.contour(x,y,chuva)
# contour(chuva_cart$z)
# filled.contour(chuva_cart$z)
# contour(chuva_cart$z,add=TRUE)
# filled.contour(chuva_cart$z, plot.axes = {
#   axis(1)
#   axis(2)
#   contour(chuva_cart$z, add = TRUE, lwd = 2)
# }
# )
# contour(x=chuva_cart$x,y=chuva_cart$y,z=chuva_cart$z)


# # transform data to a raster file without interpolation
# library(raster)
# data <- data.frame(chuva=as.vector(chuva),
#                    x=as.vector(x),
#                    y=as.vector(y))
# e <- extent(data[,2:3])
# n <- round(sqrt(length(chuva)))
# r <- raster(e, ncol = round(n/2), nrow = round(n), crs = "+proj=longlat +datum=WGS84")
# # https://www.youtube.com/watch?v=LwCEe9o0vac
# r_new <- rasterize(data[,2:3],r,data[,1], fun=mean,overwrite = TRUE,touches=TRUE,
#                    filename="~/projects/banco_de_dados/temp_2/teste")
# imagePlot(r_new,zlim=c(0,50)#,xlim=c(-47,-46.5),ylim=c(-23.8,-23.3)
#            )


# # plot data in the distance/azimute coordinate
# imagePlot(x=az,y=dist,z=chuva)


# # using ggplot and coord_polar
# az_2   <- rep(az,  times=length(dist))
# dist_2 <- rep(dist,each =length(az))
# data.frame(chuva=as.vector(chuva),az=az_2,dist=dist_2) %>%
#   ggplot(aes(az,dist,z=chuva)) +
#   geom_contour_filled(breaks = seq(0,80,10)) +
#   coord_polar() +
#   scale_x_continuous(breaks=seq(0,270,90))+
#   scale_y_continuous(breaks=seq(0,100,20))+
#   # scale_fill_manual(values = heat.colors(16),drop=FALSE) +
#   scale_fill_viridis_d(drop = FALSE) +
#   # guides(fill = guide_colorsteps(direction = "horizontal",
#   #                                barwidth = unit(par("pin")[1], "in"))) +
#   # theme(legend.position = "bottom")
#   guides(fill = guide_colorsteps(barheight = unit(par("pin")[2], "in")))
