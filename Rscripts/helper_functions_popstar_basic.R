## helper_functions_popstar.R
## Anne-Elise Nieblas (adapted from helper_functions.R by Sylvain Bonhommeau, and then from helper_functions_ss3_standard.R by AE Nieblas)
## 16/11/2018
## DESCRIPTION: This script is read by the basic_popstar_shny.Rmd for the POPSTAR basic tag diagnostic plots shiny. It loads the standardized
## .nc tag files, melts the arrays into data frames, adds/calculates new variables and assigns colnames.
## INPUTS: model, run, dir_save (where to find the .nc files)
## OUTPUTS: CPUE data frame (Fleet, Year, Exp, Obs,Dev,StDev);
##          LEN  data frame (Fleet, Year, Bin, Exp, Obs)
##          

# ########################################## CPUE DATA ######################################################
# for_ncCPUE <- function(x, dir_save){
#   nc            <-nc_open(paste(dir_save,x, sep=""))
#   exp           <-ncvar_get(nc,'Exp_cpue')
#   
#   ## replace fleet ids by fleet names for dimnames
#   globatt      <-ncatt_get(nc,0)
#   if(length(nc$dim$fleet$vals)==length(unlist(strsplit(globatt$abundance_index_id,', ')))){
#     fleetnames <- unlist(strsplit(globatt$abundance_index_names,', '))
#   }else{
#     fleetnames <- c(unlist(strsplit(globatt$abundance_index_names,', ')),'NA')
#   }
#   
#   dims          <-list(fleetnames,nc$dim$year$vals,nc$dim$season$vals)
#   coldims       <-c('Fleet','Year','Season')
#   
#   ## check for singleton dimensions
#   dr            <-NULL
#   for(d in 1:length(dims)){if(length(dims[[d]])==1){dr=c(dr,d)}}
#   dims[dr]      <-NULL
#   if(length(dr)>0){coldims       <-coldims[-dr]}
#   dimnames(exp) <-dims
#   
#   Exp          <-melt(exp)
#   obs           <-ncvar_get(nc,'Obs_cpue')
#   Obs           <-melt(obs)
#   Dev           <-Obs$value-Exp$value
#   stDev         <-Dev/sd(Dev,na.rm=T)
#   CPUE          <-cbind(Exp,Obs$value,Dev,stDev)
#   
#   colnames(CPUE)<-c(coldims,'Exp','Obs','Dev','stDev')
#   CPUE_data     <- CPUE
#   
#   return(CPUE_data)
# }
# 
# read_CPUE_data <- function(model,run, dir_save){
#   comb_run     <- expand.grid(model=model,run=run)
#   data_runs    <- cbind(File=paste("data_ss324_", comb_run$model,"_", comb_run$run, ".nc", sep=""), comb_run)
#   run_name     <- paste(run)
#   Model_run    <- paste(model,run)
#   CPUE_data    <- ddply(data_runs, .(File), function(x) data.frame(Model_run=paste(x$model, x$run),ldply(lapply(x$File, for_ncCPUE, dir_save))))[,-1]
#   return(CPUE_data)
# }
# 
# ########################################## LENGTH DATA ######################################################
# for_ncLEN <- function(x, dir_save){
#   nc           <-nc_open(paste(dir_save,x, sep=""))
#   Lobs         <-ncvar_get(nc,'Obs_len')
#   globatt      <-ncatt_get(nc,0)
#   
#   # replacing fleet ids by fleetnames for dimnames of resulting data frame.
#   if(length(nc$dim$fleet$vals)==length(unlist(strsplit(globatt$abundance_index_id,', ')))){
#     fleetnames <- unlist(strsplit(globatt$abundance_index_names,', '))
#   }else{
#     fleetnames <- c(unlist(strsplit(globatt$abundance_index_names,', ')),'NA')
#   }
#   dims         <-list(fleetnames,nc$dim$year$vals,nc$dim$season$vals,nc$dim$part$vals,nc$dim$lenbin$vals)
#   coldims      <-c('Fleet','Year','Season','Part','Bin')
#   
#   # checking for singleton dimensions
#   dr           <-NULL
#   for(d in 1:length(dims)){if(length(dims[[d]])==1){dr=c(dr,d)}}
#   dims[dr]    <-NULL
#   if(length(dr)>0){coldims     <-coldims[-dr]}
#   dimnames(Lobs)<-dims
#   
#   lobs        <-melt(Lobs)
#   Lexp        <-ncvar_get(nc,'Exp_len')
#   lexp        <-melt(Lexp)
#   Len         <-cbind(lobs,lexp$value)
#   colnames(Len)<-c(coldims,'Obs','Exp')
#   
#   LEN_data    <-Len
#   return(LEN_data)
# }
# 
# read_LEN_data <- function(model,run, dir_save){
#   comb_run    <- expand.grid(model=model,run=run)
#   data_runs   <- cbind(File=paste("data_ss324_", comb_run$model,"_", comb_run$run, ".nc", sep=""), comb_run)
#   run_name    <- paste(run)
#   Model_run   <- paste(model,run)
#   LEN_data    <- ddply(data_runs, .(File), function(x) data.frame(Model_run=paste(x$model, x$run),ldply(lapply(x$File, for_ncLEN, dir_save))))[,-1]
#   return(LEN_data)
# }

#################################### DEPTH PROFILE DATA ##########################################################
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


########################################## RECRUIT DATA ######################################################
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
  
  # CHL_1_mean  <-ncvar_get(nc,'chlorophyll_a_uncertainty_1_mean')
  # chl_1_mean  <-melt(CHL_1_mean)
  
  poly1      <-cbind(lon_1,lat_1$value,sst_1_mean$value)
  colnames(poly1)<-c('nval','time','longitude','latitude','mean_sst')
  
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


# ##################################### RECRUITMENT DEVIATION DATA #################################################
# for_ncRecDev  <- function(x, dir_save){
#   nc          <-nc_open(paste(dir_save,x, sep=""))
#   RD_val      <-ncvar_get(nc,'RecDev_val')
#   
#   dims        <-list(nc$dim$year$vals)
#   
#   ## checking for singleton dimensions
#   dr          <-NULL
#   for(d in 1:length(dims)){if(length(dims[[d]])==1){dr=c(dr,d)}}
#   dims[dr]    <-NULL
#   dimnames(RD_val)<-dims
#   
#   rd_val      <-melt(RD_val)
#   # recruitment deviation lci (recdev$Value - 1.96 * recdev$Parm_StDev)
#   RD_lci      <-ncvar_get(nc,'RecDev_lci')
#   rd_lci      <-melt(RD_lci)
#   ## recruitment deviation uci (recdev$Value + 1.96 * recdev$Parm_StDev)
#   RD_uci      <-ncvar_get(nc,'RecDev_uci')
#   rd_uci      <-melt(RD_uci)
#   
#   recdev      <-cbind(rd_val,rd_lci$value,rd_uci$value)
#   colnames(recdev)<-c('Year','value','lci','uci')
#   
#   RecDev_data <- recdev
#   return(RecDev_data)
# }
# 
# read_RecDev_data <- function(model,run, dir_save){
#   comb_run    <- expand.grid(model=model,run=run)
#   data_runs   <- cbind(File=paste("data_ss324_", comb_run$model,"_", comb_run$run, ".nc", sep=""), comb_run)
#   run_name    <- paste(run)
#   Model_run   <- paste(model,run)
#   RecDev_data <- ddply(data_runs, .(File), function(x) data.frame(Model_run=paste(x$model, x$run),ldply(lapply(x$File, for_ncRecDev, dir_save))))[,-1]
#   return(RecDev_data)
# }
# ########################################## SSB DATA ######################################################
# for_ncSSB <- function(x, dir_save){
#   nc          <- nc_open(paste(dir_save,x, sep=""))
#   SSB         <- ncvar_get(nc,'SpawnBio')
#   
#   xx=nc$dim$era$units
#   xx.val=strsplit(xx,' ')
#   xx.vals=strsplit(xx.val[[1:length(xx.val)]],':')
#   era.vals<-NULL
#   for(x in 1:length(xx.vals)){era.vals=c(era.vals,xx.vals[[x]][2])}
#   
#   dims        <- list(nc$dim$area$vals,nc$dim$year$vals,nc$dim$season$vals,era.vals)
#   coldims     <- c('Area','Year','Season','Era')
#   
#   # check for singleton dimensions
#   dr          <- NULL
#   for(d in 1:length(dims)){if(length(dims[[d]])==1){dr=c(dr,d)}}
#   dims[dr]    <- NULL
#   if(length(dr)>0){coldims     <- coldims[-dr]}
#   dimnames(SSB)<-dims
#   
#   ssb         <- melt(SSB)
#   colnames(ssb)<-c(coldims,'spbio')
#   ssb$spbio   <- ssb$spbio/1000
#   SSB_data    <- ssb
#   return(SSB_data)
# }
# 
# read_SSB_data <- function(model,run, dir_save){
#   comb_run    <- expand.grid(model=model,run=run)
#   data_runs   <- cbind(File=paste("data_ss324_", comb_run$model,"_", comb_run$run, ".nc", sep=""), comb_run)
#   run_name    <- paste(run)
#   Model_run   <- paste(model,run)
#   SSB_data    <- ddply(data_runs, .(File), function(x) data.frame(Model_run=paste(x$model, x$run),ldply(lapply(x$File, for_ncSSB, dir_save))))[,-1]
#   return(SSB_data)
# }
# 
# ########################################## BIOMASS DATA ######################################################
# for_ncBio <- function(x, dir_save){
#   nc          <- nc_open(paste(dir_save,x, sep=""))
#   Bioall      <- ncvar_get(nc,'Bio_all_ts')
#   
#   ## MUST PARSE UNITS STRING INTO VALUES - ERA
#   xx=nc$dim$era$units
#   xx.val=strsplit(xx,' ')
#   xx.vals=strsplit(xx.val[[1:length(xx.val)]],':')
#   era.vals<-NULL
#   for(x in 1:length(xx.vals)){era.vals=c(era.vals,xx.vals[[x]][2])}
#   dims        <- list(nc$dim$area$vals,nc$dim$year$vals,nc$dim$season$vals,era.vals)
#   coldims     <- c('Area','Year','Season','Era')
#   
#   # check for singleton dimensions
#   dr          <- NULL
#   for(d in 1:length(dims)){if(length(dims[[d]])==1){dr=c(dr,d)}}
#   dims[dr]    <- NULL
#   if(length(dr)>0){coldims     <- coldims[-dr]}
#   dimnames(Bioall)<-dims
#   
#   bioall      <- melt(Bioall)
#   colnames(bioall)<-c(coldims,'bio')
#   Bio_data    <- bioall
#   return(Bio_data)
# }
# 
# read_Bio_data <- function(model,run, dir_save){
#   comb_run    <- expand.grid(model=model,run=run)
#   data_runs   <- cbind(File=paste("data_ss324_", comb_run$model,"_", comb_run$run, ".nc", sep=""), comb_run)
#   run_name    <- paste(run)
#   Model_run   <- paste(model,run)
#   Bio_data    <- ddply(data_runs, .(File), function(x) data.frame(Model_run=paste(x$model, x$run),ldply(lapply(x$File, for_ncBio, dir_save))))[,-1]
#   return(Bio_data)
# }
# 
# 
# 
# ########################################## FFMSY DATA ######################################################
# for_ncffmsy <- function(x, dir_save){
#   nc          <- nc_open(paste(dir_save,x, sep=""))
#   FFmsy       <- ncvar_get(nc,'FFmsy')
#   
#   dims        <- list(nc$dim$year$vals)
#   
#   # check for singleton dimensions
#   dr          <- NULL
#   for(d in 1:length(dims)){if(length(dims[[d]])==1){dr=c(dr,d)}}
#   dims[dr]    <- NULL
#   dimnames(FFmsy)<-dims
#   
#   ffmsy       <- melt(FFmsy)
#   
#   colnames(ffmsy)<-c('Year','FFmsy')
#   ffmsy_data  <- ffmsy
#   return(ffmsy_data)
# }
# 
# read_ffmsy_data <- function(model,run, dir_save){
#   comb_run    <- expand.grid(model=model,run=run)
#   data_runs   <- cbind(File=paste("data_ss324_", comb_run$model,"_", comb_run$run, ".nc", sep=""), comb_run)
#   run_name    <- paste(run)
#   Model_run   <- paste(model,run)
#   ffmsy_data  <- ddply(data_runs, .(File), function(x) data.frame(Model_run=paste(x$model, x$run),ldply(lapply(x$File, for_ncffmsy, dir_save))))[,-1]
#   return(ffmsy_data)
# }
# 
# #################################### SIZE SELECTIVITY DATA ################################################
# for_ncSsel <- function(x, dir_save){
#   nc          <- nc_open(paste(dir_save,x, sep=""))
#   Ssel        <- ncvar_get(nc,'sizeselex')
# 
#   globatt  <- ncatt_get(nc,0)
#   end_yr   <- globatt$time_coverage_end
#   
#   ## MUST PARSE UNITS STRING INTO VALUES
#   xx=nc$dim$size_factor$units
#   xx.val=strsplit(xx,' ')
#   xx.vals=strsplit(xx.val[[1:length(xx.val)]],':')
#   size_vals<-NULL
#   for(x in 1:length(xx.vals)){size_vals=c(size_vals,xx.vals[[x]][2])}
#   
#   # replacing fleet ids by fleet names for dimnames of resulting data frame.
#   if(length(nc$dim$fleet$vals)==length(unlist(strsplit(globatt$abundance_index_id,', ')))){
#     fleetnames <- unlist(strsplit(globatt$abundance_index_names,', '))
#   }else{
#     fleetnames <- c(unlist(strsplit(globatt$abundance_index_names,', ')),'NA')
#   }
#   dims        <- list(fleetnames,nc$dim$year$vals,nc$dim$gender$vals,size_vals,nc$dim$sizeselex_size$vals)
#   coldims     <- c('Fleet','Year','Gender','Factor','Size')
#   dr          <- NULL
#   for(d in 1:length(dims)){if(length(dims[[d]])==1){dr=c(dr,d)}}
#   dims[dr]    <- NULL
#   if(length(dr)>0){coldims     <- coldims[-dr]}
#   dimnames(Ssel)<-dims
#   
#   ssel        <- melt(Ssel)
#   ssel$endyr  <- end_yr
#   colnames(ssel)<-c(coldims,'sel','endyr')
#   Ssel_data   <- ssel
#   return(Ssel_data)
# }
# 
# read_Ssel_data <- function(model,run, dir_save){
#   comb_run    <- expand.grid(model=model,run=run)
#   data_runs   <- cbind(File=paste("data_ss324_", comb_run$model,"_", comb_run$run, ".nc", sep=""), comb_run)
#   run_name    <- paste(run)
#   Model_run   <- paste(model,run)
#   Ssel_data   <- ddply(data_runs, .(File), function(x) data.frame(Model_run=paste(x$model, x$run),ldply(lapply(x$File, for_ncSsel, dir_save))))[,-1]
#   return(Ssel_data)
# }
# 
# ### AGE SELECTIVITY
# for_ncasel <- function(x, dir_save){
#   nc          <- nc_open(paste(dir_save,x, sep=""))
#   asel        <- ncvar_get(nc,'ageselex')
#   
#   globatt  <- ncatt_get(nc,0)
#   end_yr   <- globatt$time_coverage_end
#   
#   ## MUST PARSE UNITS STRING INTO VALUES
#   xx=nc$dim$age_factor$units
#   xx.val=strsplit(xx,' ')
#   xx.vals=strsplit(xx.val[[1:length(xx.val)]],':')
#   age_vals<-NULL
#   for(x in 1:length(xx.vals)){age_vals=c(age_vals,xx.vals[[x]][2])}
#   
#   # replacing fleet ids by fleet names for dimnames of resulting data frame.
#   if(length(nc$dim$fleet$vals)==length(unlist(strsplit(globatt$abundance_index_id,', ')))){
#     fleetnames <- unlist(strsplit(globatt$abundance_index_names,', '))
#   }else{
#     fleetnames <- c(unlist(strsplit(globatt$abundance_index_names,', ')),'NA')
#   }
#   dims        <- list(fleetnames,nc$dim$year$vals,nc$dim$season$vals,nc$dim$gender$vals,nc$dim$morph$vals,age_vals,nc$dim$age$vals)
#   coldims     <- c('Fleet','Year','Season','Gender','Morph','Factor','Age')
#   dr          <- NULL
#   for(d in 1:length(dims)){if(length(dims[[d]])==1){dr=c(dr,d)}}
#   dims[dr]    <- NULL
#   if(length(dr)>0){coldims     <- coldims[-dr]}
#   dimnames(asel)<-dims
#   
#   agesel        <- melt(asel)
#   agesel$endyr  <- end_yr
#   colnames(agesel)<-c(coldims,'sel','endyr')
#   Asel_data   <- agesel
#   return(Asel_data)
# }
# 
# read_Asel_data <- function(model,run, dir_save){
#   comb_run    <- expand.grid(model=model,run=run)
#   data_runs   <- cbind(File=paste("data_ss324_", comb_run$model,"_", comb_run$run, ".nc", sep=""), comb_run)
#   run_name    <- paste(run)
#   Model_run   <- paste(model,run)
#   Asel_data   <- ddply(data_runs, .(File), function(x) data.frame(Model_run=paste(x$model, x$run),ldply(lapply(x$File, for_ncasel, dir_save))))[,-1]
#   return(Asel_data)
# }
# 
# # ##### RELOADS UPDATED ss3.24.R to INFRASTRUCTURE ####
# # overwrite<-T #SET TO "TRUE" IF THE FILES ALREADY ON THE WORKSPACE SHOULD BE OVERWRITTEN
# # outputs_WS <- paste("/Home",username,"Workspace/IOTC_SS3_Shiny/Rscripts",sep="/")
# # listWS(outputs_WS) #GET THE LIST OF FILES AND FOLDERS IN ONE SUB-FOLDER
# # helper_SS3=paste('/home/anne.elise.nieblas/SS3/IOTC_SS3_Shiny/Rscripts/helper_functions_SS3.R',sep='') # FILE WITH THE FUNCTION TO WRITE OGC 19115 metadata
# # uploadWS(outputs_WS,helper_SS3,overwrite)
# # #####################################################
