## Setup really simple exmaple for testing
library(INLA)
library(VAST)
dir.create(here::here("Testing"))
example <- load_example( data_set="EBS_pollock" )
dat <- subset(example$sampling_data, Year==2013)
settings <- make_settings(n_x=100, Region=example$Region,
                          purpose="index2", bias.correct=FALSE )
settings$FieldConfig[1:2, 1:2] <- 0

fit = fit_model(settings=settings,
                working_dir=paste(here::here("Testing"), "/", sep = ""),
                Lat_i=dat$Lat, Lon_i=dat$Lon, t_i=dat$Year,
                b_i=dat$Catch_KG, a_i=dat$AreaSwept_km2)

## Also works to use a single .dll in parallel
fit.fn <- function(){
  library(VAST)
  fit_model(settings=settings,
             working_dir=paste(here::here("Testing"), "/", sep = ""),
             Lat_i=dat$Lat, Lon_i=dat$Lon, t_i=dat$Year,
             b_i=dat$Catch_KG, a_i=dat$AreaSwept_km2)
}
sfInit(parallel = TRUE, cpus = 2)
sfExportAll()
test<- sfLapply(1:2,  function(i) fit.fn())

## Adding some more complexity...
dat1<- data.frame(dat, "Group" = rep("A", nrow(dat)))
dat2<- data.frame(dat, "Group" = rep("B", nrow(dat)))
dat_all<- dat1 %>%
  rbind(., dat2) %>%
  group_by(Group) %>%
  nest()

cores_avail<- detectCores()
library(doFuture)
registerDoFuture()
plan(multisession, workers = 2)

out_folder<- paste(here::here("Testing"), "/", sep = "")
vast_files<- c(paste(here::here("Testing"), "VAST_v12_0_0.cpp", sep = "/"), paste(here::here("Testing"), "VAST_v12_0_0.so", sep = "/"), paste(here::here("Testing"), "VAST_v12_0_0.o", sep = "/"))
catch_limit<- 150
covs<- "Catch_KG"
vast_scale_func<- function(x, type){
  if(type == "JT"){
    x.out<- x/100
    return(x.out)
  } else {
    x.out<- as.numeric(scale(abs(x)))
    return(x.out)
  }
}

foreach(i = 1:nrow(dat_all)) %dopar% {
  library(VAST)
  
  # Data stuff
  dat_use<- dat_all$data[[i]] %>% 
    filter(., Catch_KG <= catch_limit)
  
  dat_use<- dat_use %>%
    mutate_at(., {{covs}}, vast_scale_func, type = "AJA")
  
  group_use<- dat_all$Group[[i]]
  
  
  # Create output folder
  outfolder<- paste(out_folder, paste("VAST", group_use, sep = "_"), sep = "")
  if(!file.exists(outfolder)){
    dir.create(outfolder)
  }
  
  # Copy over the VAST files
  file.copy(vast_files, outfolder)
  
  # Create sub-folder
  outfile<- paste(outfolder, "/", paste(group_use, sep = "_"), sep = "")
  if(!file.exists(outfile)){
    dir.create(outfile)
  }
  
  # Text file to print progress
  progress_out<- "Starting"
  write(progress_out, file = paste(outfile, "/", "progress.txt", sep = ""), append = FALSE)
  
  out<- fit_model(settings=settings,
                  "working_dir" = paste(outfile, "/", sep = ""), "CompileDir" = outfolder,
                  Lat_i=dat_use$Lat, Lon_i=dat_use$Lon, t_i=dat_use$Year,
                  b_i=dat_use$Catch_KG, a_i=dat_use$AreaSwept_km2)
  file_out<- paste("Model", group_use, sep = "_")
  saveRDS(out, file = paste(outfile, "modfit.rds", sep = "/"))
  
  # Text file to print progress
  progress_out<- max(dat_use$Catch_KG)
  write(progress_out, file = paste(outfile, "/", "progress.txt", sep = ""), append = TRUE)
}

