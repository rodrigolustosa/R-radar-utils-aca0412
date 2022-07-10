# ----------------------------------------------------------------------- #
# Information: Functions
# Created by: Rodrigo Lustosa
# Creation date:  9 Jul 2022 15:18 (GMT -03)
# ----------------------------------------------------------------------- #


# read files --------------------------------------------------------------

# le nexrad ascii
le_nexrad_ascii <- function(filein,nazim,nbins,delta_azi,gatesize){
  chuva <- matrix(0,nrow=nazim,ncol=nbins)
  az    <- vector(length=nazim)*NA
  dist  <- vector(length=nbins)*NA
  dados <- read.table(filein)
  n <- nrow(dados)
  for(i in 1:n){ # faz ate o fim do arquivo
    azimute <- dados[[i,1]]
    range   <- dados[[i,2]]
    zeff    <- dados[[i,3]]
    iazi = round(azimute/delta_azi) + 1
    jdis = floor(range/gatesize) + 1
    if (iazi > 0 & iazi <= nazim & jdis > 0 & jdis <= nbins){
      chuva[iazi,jdis] = zeff
      # az[iazi]   <- azimute
      # dist[jdis] <- range
    }
    
  }
  return(chuva)
}

le_volscanbin_xpol <- function(filein){
  # open file
  con <- file(filein, "rb")
  # read nelev, dz, dx, dy, rangestep
  bins1 <- readBin(con, integer(), n = 4, endian = "little")
  bins2 <- readBin(con, numeric(), n = 1,size = 4, endian = "little")
  nelev     <- bins1[1]
  dz        <- bins1[2]
  dx        <- bins1[3]
  dy        <- bins1[4]
  rangestep <- bins2
  # read angle_vol and volscan
  bins1 <- readBin(con,numeric(),n = nelev      ,size = 4,endian = "little")
  bins2 <- readBin(con,numeric(),n = nelev*dx*dy,size = 4,endian = "little")
  angle_vol <- bins1
  volscan <- array(bins2,dim = c(nelev,dx,dy))
  # close file
  close(con)
  # return
  return(list(volscan   = volscan,
              rangestep = rangestep,
              angle_vol = angle_vol))
}

le_volscan_cth <- function(filein){
  # open compressed file
  con <- gzcon(file(path,"rb"))
  # read date and elev
  dat  <- readBin(con, numeric(), n = 5,size = 4, endian = "little")
  elev <- readBin(con, numeric(), n = 20,size = 4, endian = "little")
  # read scans
  scan <- array(dim=c(20,360,240))
  for(i in 1:20){
    temp <- as.numeric(readBin(con, raw(), n = 240*360, endian = "little"))
    temp <- t(array(temp,dim = c(240,360)))
    scan[i,,] <- temp
  }
  scan <- 17 + scan*0.22 #convertendo para dBZ
  # close file
  close(con)
  # coerce date
  dat <- ymd_hm(str_c(dat[c(3:5,1:2)],collapse = "-"))
  yymmdd <- format(dat,"%d-%m-%y às %H:%M GMT")
  # return
  return(list(date = dat,yymmdd=yymmdd,scan = scan))
}

load_polar <- function(filein,nx,ny,nz){
  # open file
  con <- file(filein, "rb")
  # read rlon and rlat
  bins1 <- readBin(con, numeric(), n = nx*ny*nz,size = 4, endian = "little")
  bins2 <- readBin(con, numeric(), n = nx*ny*nz,size = 4, endian = "little")
  rlon <- array(bins1,dim = c(nx,ny,nz))
  rlat <- array(bins2,dim = c(nx,ny,nz))
  # close file
  close(con)
  # return
  return(list(rlon = rlon,rlat = rlat))
}
load_wsr_polar <- function(filein){
  nx <- 360; ny <- 240; nz <- 20
  return(load_polar(filein,nx,ny,nz))
}
load_xpol_polar <- function(filein){
  nx <- 361; ny <- 666; nz <- 13
  return(load_polar(filein,nx,ny,nz))
}

load_feixe_altura <- function(filein,dim1,dim2){
  # open file
  con <- file(filein, "rb")
  # read altura, distancia and range
  bins1 <- readBin(con, numeric(), n = dim1*dim2,size = 8, endian = "little")
  bins2 <- readBin(con, numeric(), n = dim1*dim2,size = 8, endian = "little")
  bins3 <- readBin(con, numeric(), n =      dim2,size = 8, endian = "little")
  altura    <- array(bins1,dim = c(dim1,dim2))
  distancia <- array(bins2,dim = c(dim1,dim2))
  range     <- bins3
  # close file
  close(con)
  # return
  return(list(altura = altura,distancia = distancia,range = range))
}
load_feixe_altura_cth <- function(filein){
  return(load_feixe_altura(filein,20,240))
}
load_feixe_altura_xpol <- function(filein){
  return(load_feixe_altura(filein,13,666))
}



# corrections -------------------------------------------------------------

limpa_ppi <- function(distancia,topo,altura,scan){
  scan_clean <- array(0,dim=c(20,360,240))
  scan_temp <- scan
  pos <- which(scan <= 17)
  scan_temp[pos] <- 0.
  #print('Limpando')
  for (teta in 0:358){
    for (r in 0:239){
      dist <- distancia[1,r+1]
      k1 <- -1
      if (dist <  60) k1 <- as.integer(dist/0.5)
      if (dist >= 60 & dist < 120) k1 <- as.integer(dist - 60) + 120
      if (dist > 120) k1 <- as.integer( (dist-120)/2 + 180)
      pos_alt <- which(altura[,r] >= topo[teta+1,k1+1]+.5)
      palt <- length(pos_alt)
      if (palt > 0){
        h0 <- pos_alt[1]
        h1 <- pos_alt[palt]
        if (h0 == 0)
          scan_temp[1,teta+1,r+1] <- 0
        else
          scan_temp[1:h0,teta+1,r+1] <- 0
        pos_chuva <- which(scan_temp[,teta+1,r+1] > 17)
        prain <- length(pos_chuva)
        if (prain > 0){
          flag <- 0
          if (prain > 1){
            for (zz in 1:(prain-1)){
              zz0 <- pos_chuva[zz]
              zz1 <- pos_chuva[zz+1]
              if ((zz1-zz0) > 1) flag <- flag+1
            } 
          }
          zz0 <- pos_chuva[1]
          zz1 <- pos_chuva[prain]
          DH <- altura[zz1,r+1]-altura[zz0,r+1]
          if (flag <= 1 & DH >= 4.)
            scan_clean[h0:h1,teta+1,r+1] <- scan_temp[h0:h1,teta+1,r+1]
        }
      }
    }
  }
  ppi_final <- array(0,dim=c(360,240))
  for (r in 0:239){
    pos_alt <- which(altura[,r+1] <= 3.0) # media de 3km
    palt <- length(pos_alt)
    if (palt > 0){
      for (teta in 0:359){
        lixo <- scan_clean[pos_alt,teta+1,r+1]
        posz <- which(lixo > 17)
        pz <- length(posz)
        if (pz > 0){
          zeff <- 10^(lixo[posz]/10)
          dbzm <- sum(zeff)/palt
          ppi_final[teta+1,r+1] <- 10*log10(dbzm)
        }
      }
    } else
      ppi_final[,r+1] = -99.  # nao tem medida neste ponto
  }
  
  return(list(scan_clean = scan_clean,ppi_final =  ppi_final))
}


# reprojections -----------------------------------------------------------

# polar to Cartesian
polar_2_cartesiana <- function(ngates,nazimu,res_azi,res_gate,x0,y0){
  x <- matrix(nrow=nazimu,ncol=ngates)
  y <- matrix(nrow=nazimu,ncol=ngates)
  for (i in 1:ngates){
    raio <- (i-1)*res_gate + res_gate/2.
    for (j in 1:nazimu){
      teta <- (j-1)*res_azi + res_azi/2.
      gama <- 450 - teta
      gama <- ifelse(gama > 360,gama - 360,gama)
      gama <- gama*pi/180
      x[j,i] = x0 + raio*cos(gama)
      y[j,i] = y0 + raio*sin(gama)
    }
  }
  return(list(x=x,y=y))
}

# The LL_ARC_DISTANCE function returns a two-element vector containing the 
# longitude and latitude [lon, lat] of a point given arc distance 
# (-pi ≤ Arc_Dist ≤ pi), and azimuth (Az), from a specified location Lon_lat0. 
# Values are in radians unless the keyword DEGREES is set.
ll_arc_distance <- function(Lon_lat0,Arc_Dist,Az,degree=FALSE){
  # algorithm from:
  # https://www.igismap.com/formula-to-find-bearing-or-heading-angle-between-two-points-latitude-longitude/
  if(degree){ # Arc_Dist is always between -pi and pi
    Lon_lat0 <- Lon_lat0*pi/180
    # Arc_Dist <- Arc_Dist*pi/180
    Az       <-       Az*pi/180
  }
  lon0 <- Lon_lat0[1]
  lat0 <- Lon_lat0[2]
  lat1 <-  asin(sin(lat0) * cos(Arc_Dist)  + cos (lat0) * sin (Arc_Dist) * cos(Az))
  lon1 <- lon0 + atan2(sin(Az)*sin(Arc_Dist)*cos(lat0), cos(Arc_Dist) - sin(lat0)*sin(lat1))
  if(degree){
    lat1 <- lat1*180/pi
    lon1 <- lon1*180/pi
  }
  return(c(lon1,lat1))
}

ll_2_dist_az <- function(Lon_lat0,Lon_lat1,degree=FALSE){
  # algorithm from:
  # https://www.igismap.com/formula-to-find-bearing-or-heading-angle-between-two-points-latitude-longitude/
  # https://www.igismap.com/haversine-formula-calculate-geographic-distance-earth/
  if(degree){
    Lon_lat0 <- Lon_lat0*pi/180
    Lon_lat1 <- Lon_lat1*pi/180
  }
  R <- 6371
  lon0 <- Lon_lat0[1]
  lat0 <- Lon_lat0[2]
  lon1 <- Lon_lat1[1]
  lat1 <- Lon_lat1[2]
  delta_lat <- lat1 - lat0
  delta_lon <- lon1 - lon0
  a <- sin(delta_lat/2)^2 + cos(lat0)*cos(lat1)*sin(delta_lon/2)^2
  Dist <- R*2*atan2(sqrt(a), sqrt(1-a))
  X <- cos(lat1)*sin(delta_lon)
  Y <- cos(lat0)*sin(lat1)-sin(lat0)*cos(lat1)*cos(delta_lon)
  Az <- atan2(X,Y)
  # Az <- acos(sin(lat1) - sin(lat0) * cos(Dist))/(cos(lat0)*sin(Dist))
  if(degree){
    Az <- Az*180/pi
  }
  return(c(Dist,Az))
}

# azimute/radius to latitude/longitude coordinates
wsr_polar <- function(lon_radar,lat_radar,res,nbins,nazim,delta_azim){
  xpolar <- matrix(nrow=nazim,ncol=nbins)
  ypolar <- matrix(nrow=nazim,ncol=nbins)
  
  res1 <- res
  lonlat <- c(lon_radar,lat_radar)
  raio <- (0:(nbins-1))*res1 + res1/2
  raio <- raio/6371#*180/pi
  
  for(teta in 1:nazim){
    alfa <- (teta-1)*delta_azim
    for (rr in 1:nbins) {
      resul <- ll_arc_distance(lonlat,raio[rr],alfa,degree=TRUE)
      xpolar[teta,rr] <- resul[1]
      ypolar[teta,rr] <- resul[2]
    }
  }
  return(list(xpolar=xpolar,ypolar=ypolar))
}



# cross sections ----------------------------------------------------------

cross_section_volscan_cth_xpol <- function(xlon,ylat,zalt,x1,y1,x2,y2,volscan,
                                           xlon_xpol,ylat_xpol,zalt_xpol,volscan_xpol,
                                           volscan_corrected){
  # xlon,ylat = 360 x 240 x 20 # azimute x range x elevacoes
  # zalt = 20 elevacoes x 240 range
  # volscan = 20 x 360 x 240
  
  slope <- (y2-y1)/(x2-x1)
  distgama <- dis1(x1,y1,x2,y2)
  dist1 <- distgama$r1
  gama1 <- distgama$gama
  
  lonlat <- c(x1,y1)
  
  npts <- as.integer(dist1) + 1
  
  xdist <- array(dim=c(npts,20))
  ydist <- array(dim=c(npts,20))
  zdist <- array(dim=c(npts,20))
  cross <- array(dim=c(npts,20))
  
  
  distancia <- array(dim=c(npts,20))
  
  resolucao_xpol <- 0.5
  
  npts_xpol <- as.integer(dist1/resolucao_xpol) + 1
  
  xdist_xpol <- array(dim=c(npts_xpol,13))
  ydist_xpol <- array(dim=c(npts_xpol,13))
  zdist_xpol <- array(dim=c(npts_xpol,13))
  cross_xpol <- array(dim=c(npts_xpol,13))
  cross_corrected <- array(dim=c(npts_xpol,13))
  distancia_xpol <- array(dim=c(npts_xpol,13))
  
  
  # CTH
  for (i in 1:npts){
    raio <- (i-1)/6371.
    resul <- ll_arc_distance(lonlat,raio,gama1,degree=TRUE)
    xpos <- resul[1]
    ypos <- y1 + slope*(xpos-x1)
    #plots,xpos,ypos,psym=1,color=61
    distancia[i,] <- i-1
    
    for (z in 1:20){
      dist1 <- sqrt((xlon[,,z]-xpos)^2 + (ylat[,,z]-ypos)^2)*111.195
      
      p1 <- which.min(dist1)
      lixo <- dist1[p1] # min(dist1)
      ipos <- (p1-1) %% 360 + 1L
      jpos <- as.integer((p1-1)/360) + 1L
      
      xdist[i,z] <- xpos
      ydist[i,z] <- ypos
      zdist[i,z] <- zalt[z,jpos]
      cross[i,z] <- volscan[z,ipos,jpos]
    }
  }
  
  
  # XPOL
  for (i in 1:npts_xpol){
    raio <- (i-1)*resolucao_xpol/6371.
    resul <- ll_arc_distance(lonlat,raio,gama1,degree=TRUE)
    xpos <- resul[1]
    ypos <- y1 + slope*(xpos-x1)
    
    distancia_xpol[i,] <- i-1
    
    for (z in 1:13){
      dist1 <- sqrt((xlon_xpol[,,z]-xpos)^2 + (ylat_xpol[,,z]-ypos)^2)*111.195
      
      p1 <- which.min(dist1)
      lixo <- dist1[p1] # min(dist1)
      ipos <- (p1-1) %% 361 + 1L
      jpos <- as.integer((p1-1)/361) + 1L
      
      xdist_xpol[i,z] <- xpos
      ydist_xpol[i,z] <- ypos
      zdist_xpol[i,z] <- zalt_xpol[z,jpos]
      cross_xpol[i,z] <- volscan_xpol[z,ipos,jpos]
      cross_corrected[i,z] <- volscan_corrected[z,ipos,jpos]
    }
  }
  return(list(cross=cross,distancia=distancia,xdist=xdist,ydist=ydist,zdist=zdist,
              cross_xpol=cross_xpol,distancia_xpol=distancia_xpol,
              xdist_xpol=xdist_xpol,ydist_xpol=ydist_xpol,zdist_xpol=zdist_xpol,
              volscan_corrected=volscan_corrected,cross_corrected=cross_corrected,
              npts=npts,npts_xpol=npts_xpol))
}



# interpolations ----------------------------------------------------------

interp_polar_2_cartesiana <- function(z_polar,res_azi,res_gate,rmax,x0=0,y0=0){
  azimu <- seq(0,360-res_azi,res_azi)
  dist  <- seq(res_gate/2,rmax,res_gate)
  ngates <- length(dist)
  nazimu <- length(azimu)
  rmax_squared <- (rmax - res_gate/2)^2
  # consider frontiers between azimute = 0 and = 360
  z_polar_new <- rbind( z_polar[nazimu,], z_polar,  z_polar[1,])
  azimu_new   <-     c(azimu[nazimu]-360,   azimu, azimu[1]+360)
  # consider frontiers in radius = 0. points in opposite azimuth
  azimu_new_op <- azimu_new - 180
  azimu_new_op <- ifelse(azimu_new_op < 0,azimu_new_op + 360,azimu_new_op)
  z_polar_new <- cbind(approx(azimu_new, z_polar_new[,1], xout = azimu_new_op, 
                              method = "linear")$y,
                       z_polar_new)
  dist_new <- c(-dist[1],dist)
  ngates_new <- length(dist_new)
  nazimu_new <- length(azimu_new)
  # chuva[nazimu,ngates]
  # new interpolated data
  res <- res_gate
  new_r <- seq(res/2,rmax,res)
  x_delta <- c(-rev(new_r),new_r)
  y_delta <- c(-rev(new_r),new_r)
  n <- length(x_delta)
  z <- matrix(nrow=n,ncol=n)
  
  for (i in 1:n){
    # print(i)
    for (j in 1:n){
      r_squared <- x_delta[i]^2 + y_delta[j]^2
      if (r_squared <= rmax_squared){
        r <- sqrt(r_squared)
        # compute azimuth
        s_x <- sign(x_delta[i])
        s_y <- sign(y_delta[j])
        ang <- (atan(y_delta[j]/x_delta[i]) +
                  s_x*(s_x-1)/2*pi)
        az <- 90 - ang*180/pi
        s_az <- sign(az)
        az   <- az + 180*(s_az-1)*s_az
        # az <- 90 - atan2(y_delta[j],x_delta[i])*180/pi
        # az <- ifelse(az < 0,az + 360,az) # more pc consuming
        
        # print(paste(x_delta[i],y_delta[j],r,az))
        
        # closest azimuths
        iazi = floor(az/delta_azi + 0.5) + 1 + 1 # considering new column in z_polar_new
        # d_az <- (azimu_new - az)
        # iazi <- which.min(abs(d_az)) # more pc consuming
        
        # closest radius
        jdis  = floor((r + gatesize/2)/gatesize + 0.5) + 1 # considering new row in z_polar_new
        # d_r <- (dist_new - r)
        # jdis <- which.min(abs(d_r)) # more pc consuming
        
        z[i,j] <- z_polar_new[iazi,jdis]
        
        
        # # closest azimuths
        # d_az <- (azimu_new - az)
        # iazi_min <- which(d_az == max(d_az[d_az <  0]))[1]
        # iazi_max <- which(d_az == min(d_az[d_az >= 0]))[1]
        # # closest radius
        # d_r <- (new_r - r)
        # jdis_min <- which(d_r == max(d_r[d_r <  0]))[1]
        # jdis_max <- which(d_r == min(d_r[d_r >= 0]))[1]
        # # interpolation with fixed radius (r min)
        # z_polar_new_r_min_az <- z_polar_new[iazi_min,jdis_min] +
        #   (z_polar_new[iazi_max,jdis_min] - z_polar_new[iazi_min,jdis_min])/
        #   (azimu_new[iazi_max] - azimu_new[iazi_min])*(az - azimu_new[iazi_min])
        # # interpolation with fixed radius (r max)
        # z_polar_new_r_max_az <- z_polar_new[iazi_min,jdis_max] +
        #   (z_polar_new[iazi_max,jdis_max] - z_polar_new[iazi_min,jdis_max])/
        #   (azimu_new[iazi_max] - azimu_new[iazi_min])*(az - azimu_new[iazi_min])
        # # interpolation with fixed azimuth (r)
        # z[i,j] <- z_polar_new_r_min_az +
        #   (z_polar_new_r_max_az - z_polar_new_r_min_az)/
        #   (dist_new[jdis_max] - dist_new[jdis_min])*(r - dist_new[jdis_min])
        
        # print(iazi_min)
        # print(jdis_min)
        # print(dim(z_polar_new))
        
      }
    }
  }
  return(list(x=x0+x_delta,y=y0+y_delta,z=z))
}

interp_polar_2_geograph <- function(z_polar,res_azi,res_gate,rmax,Lon_lat0){
  azimu <- seq(0,360-res_azi,res_azi)
  dist  <- seq(res_gate/2,rmax,res_gate)
  ngates <- length(dist)
  nazimu <- length(azimu)
  # consider frontiers between azimute = 0 and = 360
  z_polar_new <- rbind( z_polar[nazimu,], z_polar,  z_polar[1,])
  azimu_new   <-     c(azimu[nazimu]-360,   azimu, azimu[1]+360)
  # consider frontiers in radius = 0. points in opposite azimuth
  azimu_new_op <- azimu_new - 180
  azimu_new_op <- ifelse(azimu_new_op < 0,azimu_new_op + 360,azimu_new_op)
  z_polar_new <- cbind(approx(azimu_new, z_polar_new[,1], xout = azimu_new_op, 
                              method = "linear")$y,
                       z_polar_new)
  dist_new <- c(-dist[1],dist)
  ngates_new <- length(dist_new)
  nazimu_new <- length(azimu_new)
  # chuva[nazimu,ngates]
  # new interpolated data
  nbins  <- round(rmax/gatesize)
  n <- nbins*2
  lat_max <- ll_arc_distance(Lon_lat0,(rmax - res_gate/2)/6371, 0,degree=TRUE)[2]
  lon_max <- ll_arc_distance(Lon_lat0,(rmax - res_gate/2)/6371,90,degree=TRUE)[1]
  lat_min <- ll_arc_distance(Lon_lat0,(rmax - res_gate/2)/6371,180,degree=TRUE)[2]
  lon_min <- ll_arc_distance(Lon_lat0,(rmax - res_gate/2)/6371,270,degree=TRUE)[1]
  lats <- seq(lat_min,lat_max,length.out=n)
  lons <- seq(lon_min,lon_max,length.out=n)
  z <- matrix(nrow=n,ncol=n)
  for (i in 1:n){
    print(i)
    for (j in 1:n){
      # compute radius and azimuth
      r_az <- ll_2_dist_az(Lon_lat0,c(lons[i],lats[j]),degree=TRUE)
      r  <- r_az[1]
      if (r <= rmax){
        # correct sign of azimuth
        s_az <- sign(r_az[2])
        az <- r_az[2] + 180*(s_az-1)*s_az
        
        # print(paste(x_delta[i],y_delta[j],r,az))
        
        # closest azimuths
        iazi = floor(az/delta_azi + 0.5) + 1 + 1 # considering new column in z_polar_new
        
        # closest radius
        jdis  = floor((r + gatesize/2)/gatesize + 0.5) + 1 # considering new row in z_polar_new
        
        z[i,j] <- z_polar_new[iazi,jdis]
        
        
        
      }
    }
  }
  return(list(lon=lons,lat=lats,z=z))
}



# other computations ------------------------------------------------------

dis1 <- function(lon1,lat1,lon2,lat2){
  pi2 <- pi/2.
  rad <- pi/180.
  grd <- 180./pi
  degtokm <- 111.195
  
  # Latitude of North Pole
  latnp <- pi2
  
  
  lon0 <- lon1*rad
  loo1 <- lon2*rad
  
  lat0 <- lat1*rad
  laa1 <- lat2*rad
  
  # Determing which quadrant
  siglat <- lat2-lat1
  siglon <- lon2-lon1
  if (siglat >= 0){
    if (siglon >= 0){
      quad <- 1
    } else {
      quad <- 4
    }
  } else {
    if (siglon >= 0){
      quad <- 2
    } else {
      quad <- 3
    }
  }
  
  
  #print(quad)
  
  c1 <- loo1-lon0
  
  if (c1 != 0.){
    
    c1<-1./tan(c1/2.)
    
    lla1 <- laa1-lat0
    llb1 <- laa1+lat0
    
    #print(lla1);print(llb1)
    
    if (lla1 != 0.){
      yxa1 <- c1*sin(lla1/2.)/cos(llb1/2.)
      yxb1 <- c1*cos(lla1/2.)/sin(llb1/2.)
      
      yxa1 <- atan(yxa1)
      yxb1 <- atan(yxb1)
      
      r1 <- tan(lla1/2.)*sin(yxb1)/sin(yxa1)
      r1 <- atan(r1)*2.
      
      b1 <- r1
      a1 <- latnp-lat0
      c1 <- latnp-laa1
      s <- (a1+b1+c1)/2.
      r <- abs(sin(s-a1)*sin(s-b1)*sin(s-c1)/sin(s))
      
      r <- sqrt(r)
      gama <- 2.*atan(r/sin(s-c1))
      
      gama <- gama*grd
      zeta <- gama
      
      # if quad eq 4 then print,"**** ",gama
      
      if (quad == 2) gama <- 180 + gama
      if (quad == 1) gama <- 180 + gama
      if (quad == 3) gama <- 180.-gama
      if (quad == 4) gama <- 180.-gama
      
    } else {
      # Applying the Right Spherical Triangle rule
      lla1 <- latnp-laa1
      c1 <- loo1-lon0
      r1 <- sin(c1)*sin(lla1)
      r1 <- asin(r1)
      if (siglon >= 0) gama <- 90.
      if (siglon <  0) gama <- 270.
      zeta <- gama
      
    }
    
  } else {
    # In the Latitude Great Circle
    r1 <- abs(laa1-lat0)
    if (siglat >= 0) gama <- 0.
    if (siglat <  0) gama <- 180.
    zeta <- gama
  }
  
  # Converting to Degrees
  r1 <- r1*grd
  
  # Converting to Km (picking just the positive)
  r1 <- abs(r1*degtokm)
  
  return(list(r1=r1,gama=gama))
}



# shiny app ---------------------------------------------------------------


ui <- fluidPage(
  fluidRow(column(6,plotOutput("plot1")),
           column(6,plotOutput("plot2",click = "plot_click"))),
  verbatimTextOutput("info")
)

server <- function(input, output, session) {
  # clicks on the plot
  clicks <- reactiveValues(x=NULL,y=NULL)
  # compute closest possible coordinate from click
  observe({
    if(!is.null(input$plot_click)){
      clicks$x <- input$plot_click$x
      clicks$y <- input$plot_click$y
      dist1 <- sqrt((rlon_xpol[,,1]-clicks$x)^2 + (rlat_xpol[,,1]-clicks$y)^2)*111.195
      pmin <- which.min(dist1)
      rmin <- dist1[pmin]
      azimute <- (pmin-1) %% 361 + 1
      idist   <- (pmin-1)/361 + 1
      xlon1 <<- rlon_xpol[azimute,idist,1]
      ylat1 <<- rlat_xpol[azimute,idist,1]
    }
  })
  # plot CTH
  output$plot1 <- renderPlot({
    imagePlot(xlixo,ylixo,cth_clean$ppi_final,legend.lab = "Z (dBZ)",
              xlim=c(lon0,lon1),ylim=c(lat0,lat1),zlim=c(0.1,60))
    # map("world", add = TRUE,resolution = 0)
    # plot(states$geom,add=TRUE)
    plot(brazil,add=TRUE)
    points(lon_xpol,lat_xpol,pch=8)
    points(lon_cth,lat_cth,pch=8)
    title(paste("CTH -",cth$yymmdd))
    box()
    # plot click
    if(!is.null(clicks$x)){
      points(xlon1,ylat1,pch=8)
      xx1=c(xlon0,xlon1)
      yy1=c(ylat0,ylat1)
      lines(xx1,yy1,col="magenta",lwd=2)
    }
  })
  # plot XPOL
  output$plot2 <- renderPlot({
    imagePlot(rlon_xpol[,,2],rlat_xpol[,,2],xpol$volscan[2,,],
              legend.lab = "Z (dBZ)",
              xlim=c(lon0,lon1),ylim=c(lat0,lat1),zlim=c(0.1,60))
    plot(brazil,add=TRUE)
    points(lon_xpol,lat_xpol,pch=8)
    points(lon_cth,lat_cth,pch=8)
    title(paste("XPOL -",cth$yymmdd))
    box()
    # plot click
    if(!is.null(clicks$x)){
      points(xlon1,ylat1,pch=8)
      xx1=c(xlon0,xlon1)
      yy1=c(ylat0,ylat1)
      lines(xx1,yy1,col="magenta",lwd=2)
    }
  })
  # explanation
  output$info <- renderText({
    "Click no Azimute do XPOL que vc quer fazer a Secao transversal. Feche esta pagina apos selecionar."
  })
}
