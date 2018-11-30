## helper_functions_popstar.R
## Anne-Elise Nieblas (adapted from helper_functions.R by Sylvain Bonhommeau, and then from helper_functions_ss3_standard.R by AE Nieblas)
## 16/11/2018
## DESCRIPTION: This script is read by the basic_popstar_shny.Rmd for the POPSTAR basic tag diagnostic plots shiny. It loads the standardized
## .nc tag files, melts the arrays into data frames, adds/calculates new variables and assigns colnames.
## INPUTS: model, run, dir_save (where to find the .nc files)
## OUTPUTS: CPUE data frame (Fleet, Year, Exp, Obs,Dev,StDev);
##          LEN  data frame (Fleet, Year, Bin, Exp, Obs)
##          


#################################### POINTS DATA ##########################################################
for_ncPoints <- function(x, dir_save){
  nc          <-nc_open(paste0(dir_save,x,'.nc'))
  DEPTH_val   <-ncvar_get(nc,'depth')
  
  # dims        <-list(nc$dim$TIME$vals)
  tt          <- nc$var$latitude$dim[[1]]$vals
  
  dd          <- as.POSIXct(tt,origin='1970-01-01',tz='GMT')
  dims        <-list(dd)
  
  ## checking for singleton dimensions
  dr          <-NULL
  for(d in 1:length(dims)){if(length(dims[[d]])==1){dr=c(dr,d)}}
  dims[dr]    <-NULL
  dimnames(DEPTH_val)<-dims
  
  depth_val      <-melt(DEPTH_val)
  time_posix    <-as.POSIXct(depth_val$Var1,origin='1970-01-01',tz='GMT')
  depth         <-depth_val$value*(-1)
  
  # recruitment deviation lci (recdev$Value - 1.96 * recdev$Parm_StDev)
  LAT_val      <-ncvar_get(nc,'latitude')
  lat_val      <-melt(LAT_val)
  ## recruitment deviation uci (recdev$Value + 1.96 * recdev$Parm_StDev)
  LON_val      <-ncvar_get(nc,'longitude')
  lon_val      <-melt(LON_val)
  
  TEMP_val    <-ncvar_get(nc,'sea_water_temperature')
  temp_val    <-melt(TEMP_val)
  
  SST_val     <-ncvar_get(nc,'sea_surface_temperature')
  sst_val     <-melt(SST_val)
  
  # CHL_val     <-ncvar_get(nc,'chlorophyll_a')
  # chl_val     <-melt(CHL_val)
  
  
  points      <-cbind(time_posix,lon_val,lat_val$value,depth,temp_val$value,sst_val$value)
  colnames(points)<-c('time','numtime','longitude','latitude','depth','in_situ_temp','sst')
  
  
  Points_data <- points
  return(Points_data)
}

read_POINTS_data <- function(model,run, dir_save){
  comb_run     <- expand.grid(model=model)
  data_runs    <- cbind(File=paste(comb_run$model,sep=""), comb_run)
  run_name     <- paste(run)
  Model_run    <- paste(model,run)
  POINTS_data    <- ddply(data_runs, .(File), function(x) data.frame(Model_run=paste(x$model, x$run),ldply(lapply(x$File, for_ncPoints, dir_save))))[,-1]
  return(POINTS_data)
}


########################################## POLY1 DATA ######################################################
for_ncPoly1 <- function(x, dir_save){
  nc          <-nc_open(paste0(dir_save,x,'.nc'))
  # nc          <-nc_open(paste(dir_save,x, sep=""))
  # R0<-ncvar_get(nc,'Recruit_0')
  LON_1<-ncvar_get(nc,'longitude_uncertainty_1')
  
  dims        <-list(nc$dim$n_values_1$vals,nc$dim$time_uncertainty$vals)
  coldims     <-c('nvalues','time_UncertaintyArea')
  
  # check for singleton dimensions
  dr          <-NULL
  for(d in 1:length(dims)){if(length(dims[[d]])==1){dr=c(dr,d)}}
  dims[dr]    <-NULL
  if(length(dr)>0){coldims     <-coldims[-dr]}
  dimnames(LON_1)<-dims
  
  lon_1          <-melt(LON_1)
  
  
  LAT_1      <-ncvar_get(nc,'latitude_uncertainty_1')
  lat_1      <-melt(LAT_1)
  
  SST_1_mean  <-ncvar_get(nc,'sea_surface_temperature_uncertainty_1_mean')
  sst_1_mean <-melt(SST_1_mean)
  
  SST_1_min  <-ncvar_get(nc,'sea_surface_temperature_uncertainty_1_min')
  sst_1_min <-melt(SST_1_min)
  
  SST_1_max  <-ncvar_get(nc,'sea_surface_temperature_uncertainty_1_max')
  sst_1_max <-melt(SST_1_max)
  
  SST_1_q25  <-ncvar_get(nc,'sea_surface_temperature_uncertainty_1_q25')
  sst_1_q25 <-melt(SST_1_q25)
  
  SST_1_q75  <-ncvar_get(nc,'sea_surface_temperature_uncertainty_1_q75')
  sst_1_q75 <-melt(SST_1_q75)
  
  SST_1_sd  <-ncvar_get(nc,'sea_surface_temperature_uncertainty_1_sd')
  sst_1_sd <-melt(SST_1_sd)
  
  # CHL_1_mean  <-ncvar_get(nc,'chlorophyll_a_uncertainty_1_mean')
  # chl_1_mean  <-melt(CHL_1_mean)
  
  poly1      <-cbind(lon_1,lat_1$value,sst_1_mean$value,sst_1_min$value,sst_1_max$value,sst_1_sd$value,sst_1_q25$value,sst_1_q75$value)
  colnames(poly1)<-c('nval','time','longitude','latitude','mean_sst','min_sst','max_sst','sd_sst','q25_sst','q75_sst')
  
  Poly1_data     <-poly1
  return(Poly1_data)
}

read_POLY1_data <- function(model,run, dir_save){
  comb_run     <- expand.grid(model=model)
  data_runs    <- cbind(File=paste(comb_run$model, sep=""), comb_run)
  run_name     <- paste(run)
  Model_run    <- paste(model,run)
  POLY1_data    <- ddply(data_runs, .(File), function(x) data.frame(Model_run=paste(x$model, x$run),ldply(lapply(x$File, for_ncPoly1, dir_save))))[,-1]
  return(POLY1_data)
}

########################################## POLY2 DATA ######################################################
for_ncPoly2 <- function(x, dir_save){
  nc          <-nc_open(paste0(dir_save,x,'.nc'))
  # nc          <-nc_open(paste(dir_save,x, sep=""))
  # R0<-ncvar_get(nc,'Recruit_0')
  LON_2<-ncvar_get(nc,'longitude_uncertainty_2')
  
  dims        <-list(nc$dim$n_values_2$vals,nc$dim$time_uncertainty$vals)
  coldims     <-c('nvalues','time_UncertaintyArea')
  
  # check for singleton dimensions
  dr          <-NULL
  for(d in 1:length(dims)){if(length(dims[[d]])==1){dr=c(dr,d)}}
  dims[dr]    <-NULL
  if(length(dr)>0){coldims     <-coldims[-dr]}
  dimnames(LON_2)<-dims
  
  lon_2          <-melt(LON_2)
  
  
  LAT_2      <-ncvar_get(nc,'latitude_uncertainty_2')
  lat_2      <-melt(LAT_2)
  
  SST_2_mean  <-ncvar_get(nc,'sea_surface_temperature_uncertainty_2_mean')
  sst_2_mean <-melt(SST_2_mean)
  
  SST_2_min  <-ncvar_get(nc,'sea_surface_temperature_uncertainty_2_min')
  sst_2_min <-melt(SST_2_min)
  
  SST_2_max  <-ncvar_get(nc,'sea_surface_temperature_uncertainty_2_max')
  sst_2_max <-melt(SST_2_max)
  
  SST_2_q25  <-ncvar_get(nc,'sea_surface_temperature_uncertainty_2_q25')
  sst_2_q25 <-melt(SST_2_q25)
  
  SST_2_q75  <-ncvar_get(nc,'sea_surface_temperature_uncertainty_2_q75')
  sst_2_q75 <-melt(SST_2_q75)
  
  SST_2_sd  <-ncvar_get(nc,'sea_surface_temperature_uncertainty_2_sd')
  sst_2_sd <-melt(SST_2_sd)
  
  # CHL_2_mean  <-ncvar_get(nc,'chlorophyll_a_uncertainty_2_mean')
  # chl_2_mean  <-melt(CHL_2_mean)
  
  poly2      <-cbind(lon_2,lat_2$value,sst_2_mean$value,sst_2_min$value,sst_2_max$value,sst_2_sd$value,sst_2_q25$value,sst_2_q75$value)
  colnames(poly2)<-c('nval','time','longitude','latitude','mean_sst','min_sst','max_sst','sd_sst','q25_sst','q75_sst')
  
  Poly2_data     <-poly2
  return(Poly2_data)
}

read_POLY2_data <- function(model,run, dir_save){
  comb_run     <- expand.grid(model=model)
  data_runs    <- cbind(File=paste(comb_run$model, sep=""), comb_run)
  run_name     <- paste(run)
  Model_run    <- paste(model,run)
  POLY2_data    <- ddply(data_runs, .(File), function(x) data.frame(Model_run=paste(x$model, x$run),ldply(lapply(x$File, for_ncPoly2, dir_save))))[,-1]
  return(POLY2_data)
}


########################################## POLY3 DATA ######################################################
for_ncPoly3 <- function(x, dir_save){
  nc          <-nc_open(paste0(dir_save,x,'.nc'))
  # nc          <-nc_open(paste(dir_save,x, sep=""))
  # R0<-ncvar_get(nc,'Recruit_0')
  LON_3<-ncvar_get(nc,'longitude_uncertainty_3')
  
  dims        <-list(nc$dim$n_values_3$vals,nc$dim$time_uncertainty$vals)
  coldims     <-c('nvalues','time_UncertaintyArea')
  
  # check for singleton dimensions
  dr          <-NULL
  for(d in 1:length(dims)){if(length(dims[[d]])==1){dr=c(dr,d)}}
  dims[dr]    <-NULL
  if(length(dr)>0){coldims     <-coldims[-dr]}
  dimnames(LON_3)<-dims
  
  lon_3          <-melt(LON_3)
  
  
  LAT_3      <-ncvar_get(nc,'latitude_uncertainty_3')
  lat_3      <-melt(LAT_3)
  
  SST_3_mean  <-ncvar_get(nc,'sea_surface_temperature_uncertainty_3_mean')
  sst_3_mean <-melt(SST_3_mean)
  
  SST_3_min  <-ncvar_get(nc,'sea_surface_temperature_uncertainty_3_min')
  sst_3_min <-melt(SST_3_min)
  
  SST_3_max  <-ncvar_get(nc,'sea_surface_temperature_uncertainty_3_max')
  sst_3_max <-melt(SST_3_max)
  
  SST_3_q25  <-ncvar_get(nc,'sea_surface_temperature_uncertainty_3_q25')
  sst_3_q25 <-melt(SST_3_q25)
  
  SST_3_q75  <-ncvar_get(nc,'sea_surface_temperature_uncertainty_3_q75')
  sst_3_q75 <-melt(SST_3_q75)
  
  SST_3_sd  <-ncvar_get(nc,'sea_surface_temperature_uncertainty_3_sd')
  sst_3_sd <-melt(SST_3_sd)
  
  # CHL_3_mean  <-ncvar_get(nc,'chlorophyll_a_uncertainty_3_mean')
  # chl_3_mean  <-melt(CHL_3_mean)
  
  poly3      <-cbind(lon_3,lat_3$value,sst_3_mean$value,sst_3_min$value,sst_3_max$value,sst_3_sd$value,sst_3_q25$value,sst_3_q75$value)
  colnames(poly3)<-c('nval','time','longitude','latitude','mean_sst','min_sst','max_sst','sd_sst','q25_sst','q75_sst')
  
  Poly3_data     <-poly3
  return(Poly3_data)
}

read_POLY3_data <- function(model,run, dir_save){
  comb_run     <- expand.grid(model=model)
  data_runs    <- cbind(File=paste(comb_run$model, sep=""), comb_run)
  run_name     <- paste(run)
  Model_run    <- paste(model,run)
  POLY3_data   <- ddply(data_runs, .(File), function(x) data.frame(Model_run=paste(x$model, x$run),ldply(lapply(x$File, for_ncPoly3, dir_save))))[,-1]
  return(POLY3_data)
}

