# ----------------------------------------------------------------------- #
# Information: 
# Created by: Rodrigo Lustosa
# Creation date:  8 Jul 2022 12:20 (GMT -03)
# ----------------------------------------------------------------------- #

# packages
library(tidyverse)
library(fields)

# directory and file names
dir_data <- "data"
file_functions <- "functions.R"
file_radar_pattern <- "Nexrad/*.dat"


# functions ---------------------------------------------------------------

source(file_functions)


# radar information -------------------------------------------------------

# ##########################################################
# INSIRA as Posicoes nazimute (iazi) e ngate (jgate) que vc calculou
iazi = 50   + 1 # mais um pois as coordenadas de uma matriz no R vao de 1 a n,
jgate = 536 + 1 # enquanto que no IDL vao de 0 a n-1
# #######################################################

# Nexrad
rmax      <- 100.     # distancia maxima
delta_azi <- 0.5      # amostragem em azimute
gatesize  <- 0.125    # resolucao do gate
x0        <- 0.       # posicao do radar em coordenadas cartesiana
y0        <- 0.       

lat_radar <- -(23 + 33/60 +  42.2/3600)
lon_radar <- -(46 + 44/60 +  06.4/3600)

;   ################# Dimensoes dos campos de chuva e matriz de navegacao
nbins = round(rmax/gatesize) #+ 1
nazim = round(360./delta_azi)#+ 1
# nbins = fix(rmax/gatesize) 
# nazim = fix(360./delta_azi) + 1

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
path <- file.path(dir_data,file_radar_pattern)
file_radar <- Sys.glob(path)
pfiles <- length(files_nexrad)
for(ff in 1:pfiles){
  chuva <- le_nexrad_ascii(file_radar[ff],nazim,nbins,delta_azi,gatesize)
  print(file_radar[ff])
  print(chuva[iazi,jgate])
  imagePlot(xpolar,ypolar,chuva,zlim=c(0,60),#zlim=c(-99,60),
            # xlim=c(-46.8,-46.7),ylim=c(-23.6,-23.5),
            main="Geográfica (lat/lon)",xlab="Longitude (°)",ylab="Latitude (°)",
            legend.lab="Z(dBZ)")
  points(xpolar[iazi,jgate],ypolar[iazi,jgate],pch=8,col="white")
  
}
