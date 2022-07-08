

# polar to cartesiana
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
  if(degree){
    Lon_lat0 <- Lon_lat0*pi/180
    Arc_Dist <- Arc_Dist*pi/180
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


wsr_polar <- function(lon_radar,lat_radar,res,nbins,nazim,delta_azim){
  xpolar <- matrix(nrow=nazim,ncol=nbins)
  ypolar <- matrix(nrow=nazim,ncol=nbins)
  
  res1 <- res
  lonlat <- c(lon_radar,lat_radar)
  raio <- (0:(nbins-1))*res1 + res1/2
  raio <- raio/6371*180/pi
  
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
  lat_max <- ll_arc_distance(Lon_lat0,(rmax - res_gate/2)/6371/pi*180, 0,degree=TRUE)[2]
  lon_max <- ll_arc_distance(Lon_lat0,(rmax - res_gate/2)/6371/pi*180,90,degree=TRUE)[1]
  lat_min <- ll_arc_distance(Lon_lat0,(rmax - res_gate/2)/6371/pi*180,180,degree=TRUE)[2]
  lon_min <- ll_arc_distance(Lon_lat0,(rmax - res_gate/2)/6371/pi*180,270,degree=TRUE)[1]
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
