# ----------------------------------------------------------------------- #
# Information: 
# Created by: Rodrigo Lustosa
# Creation date:  9 Jul 2022 12:03 (GMT -03)
# ----------------------------------------------------------------------- #

# packages
library(tidyverse)
library(fields)
library(lubridate)
# library(maps)
library(geobr)
library(shiny)


# directory and file names
dir_data <- "data/volscan"
dir_output <- "output"
file_functions <- "functions.R"
file_wsr_polar  <- "Navegacao_3D_CTH.dat"
file_xpol_polar <- "Navegacao_3D_XPOL_Vale.dat"
file_feixe_altura_cth  <- "feixe_altura_range_cth.dat"
file_feixe_altura_xpol <- "feixe_altura_range_xpol.dat"
file_topography <- "topografia_feixe_cth_polar_240km.dat"
file_volscan_cth   <- "cscan_201202101837.dat.gz"
file_volscan_xpol1 <- "level0_2012021018360400dBZ.volscan.bin"
file_volscan_xpol2 <- "level1_2012021018360400dBZ.volscan.bin"

# functions ---------------------------------------------------------------

source(file_functions)

# readings and radar information ------------------------------------------

# Carregando coordenadas polares
lon_cth <- -45.9722
lat_cth <- -23.60
path <- file.path(dir_data,file_wsr_polar)
rlonlat <- load_wsr_polar(path)
rlon <- rlonlat$rlon
rlat <- rlonlat$rlat

lon_xpol = -45.952207
lat_xpol = -23.208683
path <- file.path(dir_data,file_xpol_polar)
rlonlat_xpol <- load_xpol_polar(path)
rlon_xpol <- rlonlat_xpol$rlon
rlat_xpol <- rlonlat_xpol$rlat

# Calculando a distancia e altura do feixe do radar do CTH
# com a configuracao atual
path <- file.path(dir_data,file_feixe_altura_cth)
feixe_alt_cth  <- load_feixe_altura_cth(path)
path <- file.path(dir_data,file_feixe_altura_xpol)
feixe_alt_xpol <- load_feixe_altura_xpol(path)

# Carrega a Topografia
path <- file.path(dir_data,file_topography)
topo <- readBin(path, numeric(), n = 360*240,size = 4, endian = "little")
topo <- array(topo,dim=c(360,240))
topo <- topo/1000. # meters to kilometers

# Lendo o Volscan que ja foi descompacto do formato McGuill
path <- file.path(dir_data,file_volscan_cth)
cth <- le_volscan_cth(path)
cth_clean <- limpa_ppi(feixe_alt_cth$distancia,topo,feixe_alt_cth$altura,cth$scan)

path <- file.path(dir_data,file_volscan_xpol1)
xpol <- le_volscanbin_xpol(path)

path <- file.path(dir_data,file_volscan_xpol2)
xpolc <- le_volscanbin_xpol(path)


# plots information -------------------------------------------------------

lon0 <- min(rlon)+1
lon1 <- max(rlon)-1.
lat0 <- min(rlat)+1.
lat1 <- max(rlat)-1.

lixo  <- array(dim=c(360,240))
xlixo <- rlon[,,1]
ylixo <- rlat[,,1]

# get map contours
# states <- geobr::read_state(simplified = FALSE)
brazil <- geobr::read_country()



# user input --------------------------------------------------------------

# get user inputs from click on the map
xlon0 = rlon_xpol[1,1,1]
ylat0 = rlat_xpol[1,1,1]
xlon1 <- NULL
ylat1 <- NULL
shinyApp(ui, server)

# get cross section from input
cross_section <- 
  cross_section_volscan_cth_xpol(rlon,rlat,feixe_alt_cth$altura,
                                 xlon0,ylat0,xlon1,ylat1,cth_clean$scan_clean,
                                 rlon_xpol,rlat_xpol,feixe_alt_xpol$altura,
                                 xpol$volscan,xpolc$volscan)



# plots -------------------------------------------------------------------

xx1=c(xlon0,xlon1)
yy1=c(ylat0,ylat1)

# CTH
file_output <- str_c('ppi_cth_', format(cth$date,"%y-%m-%d_%H-%M"), '.png')
path <- file.path(dir_output,file_output)
png(path,width = 16,height = 16,units="cm",res=100)
imagePlot(xlixo,ylixo,cth_clean$ppi_final,
          legend.lab = "Z (dBZ)",xlab="Longitude",ylab="Latitude",
          xlim=c(lon0,lon1),ylim=c(lat0,lat1),zlim=c(0.1,71))
# map("world", add = TRUE,resolution = 0)
# plot(states$geom,add=TRUE)
points(lon_xpol,lat_xpol,pch=8)
points(lon_cth,lat_cth,pch=8)
points(xlon1,ylat1,pch=8)
lines(xx1,yy1,col="magenta",lwd=2)
plot(brazil,add=TRUE)
title(paste("CTH -",cth$yymmdd))
box()
dev.off()

# XPOL
file_output <- str_c('ppi_xpol_', format(cth$date,"%y-%m-%d_%H-%M"), '.png')
path <- file.path(dir_output,file_output)
png(path,width = 16,height = 16,units="cm",res=100)
imagePlot(rlon_xpol[,,2],rlat_xpol[,,2],xpol$volscan[2,,],
          legend.lab = "Z (dBZ)",xlab="Longitude",ylab="Latitude",
          xlim=c(lon0,lon1),ylim=c(lat0,lat1),zlim=c(0.1,71))
points(lon_xpol,lat_xpol,pch=8)
points(lon_cth,lat_cth,pch=8)
points(xlon1,ylat1,pch=8)
lines(xx1,yy1,col="magenta",lwd=2)
plot(brazil,add=TRUE)
title(paste("XPOL -",cth$yymmdd))
box()
dev.off()


# CTH
file_output <- str_c('cross_cth_', format(cth$date,"%y-%m-%d_%H-%M"), '.png')
path <- file.path(dir_output,file_output)
png(path,width = 16,height = 16,units="cm",res=100)
imagePlot(cross_section$xdist,cross_section$zdist,cross_section$cross,
          legend.lab = "Z (dBZ)",xlab="Longitude",ylab="Height (km)",
          ylim=c(0,15),zlim=c(.1,71))
lines(cross_section$xdist[,20],cross_section$zdist[,20])
lines(cross_section$xdist[,1],cross_section$zdist[,1])
title(paste("CTH -",cth$yymmdd))
box()
dev.off()

# XPOL
file_output <- str_c('cross_xpol_', format(cth$date,"%y-%m-%d_%H-%M"), '.png')
path <- file.path(dir_output,file_output)
png(path,width = 16,height = 16,units="cm",res=100)
imagePlot(cross_section$xdist_xpol,cross_section$zdist_xpol,
          cross_section$cross_xpol,ylim=c(0,15),zlim=c(.1,71),
          legend.lab = "Z (dBZ)",xlab="Longitude",ylab="Height (km)",)
lines(cross_section$xdist_xpol[,13],cross_section$zdist_xpol[,13])
lines(cross_section$xdist_xpol[,1],cross_section$zdist_xpol[,1])
title(paste("XPOL -",cth$yymmdd))
box()
dev.off()

# XPOLC
file_output <- str_c('cross_xpolc_', format(cth$date,"%y-%m-%d_%H-%M"), '.png')
path <- file.path(dir_output,file_output)
png(path,width = 16,height = 16,units="cm",res=100)
imagePlot(cross_section$xdist_xpol,cross_section$zdist_xpol,
          cross_section$cross_corrected,ylim=c(0,15),zlim=c(.1,71),
          legend.lab = "Z (dBZ)",xlab="Longitude",ylab="Height (km)")
lines(cross_section$xdist_xpol[,13],cross_section$zdist_xpol[,13])
lines(cross_section$xdist_xpol[,1],cross_section$zdist_xpol[,1])
title(paste("XPOLC -",cth$yymmdd))
box()
dev.off()


# save data ---------------------------------------------------------------

# save cross section from xpol
path <- file.path(dir_output,'cross_section_xpol.dat')
con <- file(path,"w")
for (z in 1:13){
  writeLines(paste("Elevacao =",formatC(z,width=2,flag = 0)),con)
  for (i in 1:cross_section$npts_xpol){
    writeLines(paste(sprintf("%5.1f",(i-1)*0.5),
                     sprintf("%5.2f",cross_section$zdist_xpol[i,z]),
                     sprintf("%6.2f",cross_section$cross_xpol[i,z]),
                     sprintf("%6.2f",cross_section$cross_corrected[i,z])
    ),con)
  }
}
close(con)

# save cross section from cth
path <- file.path(dir_output,'cross_section_cth.dat')
con <- file(path,"w")
for (z in 1:20){
  writeLines(paste("Elevacao =",formatC(z,width=2,flag = 0)),con)
  for (i in 1:cross_section$npts){
    writeLines(paste(sprintf("%5.1f",(i-1)*0.5),
                     sprintf("%5.2f",cross_section$zdist[i,z]),
                     sprintf("%6.2f",cross_section$cross[i,z])
    ),con)
  }
}
close(con)



# backups -----------------------------------------------------------------

# imagePlot(cth$scan[2,,],zlim=c(0,60))
# 
# imagePlot(cth_clean$scan_clean[2,,],zlim=c(0,60))
# imagePlot(cth_clean$ppi_final[,],zlim=c(0,60))
# 
# imagePlot(xpolc$volscan[2,,],zlim=c(0,60))
# imagePlot(xpol$volscan[2,,],zlim=c(0,60))
# max(xpol$volscan)
# 
# xpol$volscan[1,,]
# rlon_xpol[,,1]
# rlat_xpol[,,1]
