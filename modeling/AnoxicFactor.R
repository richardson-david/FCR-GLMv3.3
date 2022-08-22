#code from Robert Ladwig, 18 Aug 2022
#modified by CCC for FCR

pacman::p_load(tidyverse, lubridate, ncdf4, GLMr, glmtools,ncdf4,pracma)

sim_folder<-"/Users/cayelan/Dropbox/ComputerFiles/SCC/FCR-GLMv3.3"
setwd(sim_folder)
nc_file <- file.path(sim_folder, 'output/output.nc') #defines the output.nc file 
mod_nc <- nc_file
var="OXY_oxy"
int_step = 0.1 #m of interpolation (0.1 = 10 cm)
nml_file <- paste0(sim_folder,"/glm3.nml")
nml <- read_nml(nml_file) 


#get oxygen data using Robert's modified GLM-AED output
get_wattemp <- function(mod_nc, reference = "surface", var, int_step){
  ncin <- nc_open(mod_nc)
  watdep <- ncvar_get(ncin, "z")
  wattemp <- ncvar_get(ncin, var)
  
  time <- ncvar_get(ncin, "time")
  time.units <- ncatt_get(ncin, "time", "units")
  #sub("^.*\\s2","",time.units$value)
  time.start <- as.POSIXct(strptime(sub("hours since ","",time.units$value), 
                                    format = "%Y-%m-%d %H:%M:%S"))
  datetime <- time.start + time*3600
  
  layer <- ncvar_get(ncin, "NS")
  nc_close(ncin)
  watdep[which(watdep == max(watdep))] <- NaN
  wattemp[which(wattemp == max(wattemp))] <- NaN
  
  sim.watdep <- 0.0*watdep - 999
  for (i in 1:length(datetime)){
    max_depth <- watdep[layer[i],i]
    sim.watdep[1,i] <- max_depth - watdep[1,i]/2 
    for (j in 2:layer[i]){
      sim.watdep[j, i] <- max_depth - (watdep[j,i] + watdep[j-1, i ])/2 
    }
  }
  
  int.dep <- rev(seq(0.25,round(max(watdep,na.rm = TRUE)),int_step))
  sim.wattemp <- matrix(0, nrow = length(int.dep), ncol= length(datetime))
  for (i in 1:length(datetime)){
    sim.approx <- approx(na.omit(sim.watdep[,i]), na.omit(wattemp[,i]), int.dep)
    sim.wattemp[1:length(int.dep),i] <- sim.approx$y
  }
  
  return(list('sim' = sim.wattemp, 'time' = datetime, 'depth' = int.dep))
}

#now run the function above for oxygen 
oxy_data <- get_wattemp(mod_nc, reference = 'surface', var ='OXY_oxy',int_step=0.1)

# get oxycline depth (1 mg/L)
crit <- 1.0 * 1000/32 #could be 2 mg/L or other!
obs.crit <- crit * 32/1000

#get hypsometry
#nml$morphometry$H
H <- abs(nml$morphometry$H - max(nml$morphometry$H))
A <- nml$morphometry$A

#setting empty vectors to be filled
min.crit <- c()
oxy.dep <- c()
oxy.area <- c()
oxy.value <- c()
oxy_data$doy = yday(oxy_data$time)

for (ii in 1:ncol(oxy_data$sim)){
  if (oxy_data$doy[ii] > 140 && oxy_data$doy[ii] <300){ #set days of year to focus on days 140-300
    if (min(oxy_data$sim[,ii], na.rm = TRUE) <= crit ){#(min(abs(oxy_data$sim[,ii] - crit), na.rm = TRUE) <= crit ){#(min(oxy_data$sim[,ii], na.rm = TRUE) < crit ){#
      min.crit <- (which(oxy_data$sim[,ii] <= crit))
      if (max(diff(min.crit)) >1){
        badpt <- which(diff(min.crit) > 1)
        min.crit <- min.crit[-c((min(badpt)+1):length(min.crit))]
      }
      oxy.dep <- append(oxy.dep, oxy_data$dep[max(min.crit)])
      #print(oxy.dep)
      oxy.area <- append(oxy.area, pracma::interp1(H, A, oxy.dep[ii], method = "spline"))
      oxy.value = append(oxy.value, oxy_data$sim[which(abs(oxy_data$sim[,ii] - crit) == min(abs(oxy_data$sim[,ii] - crit), na.rm = TRUE)),ii])
    } else {
      min.crit <- append(min.crit, NA)
      oxy.dep <- append(oxy.dep, NA)
      oxy.area <- append(oxy.area, NA)
      oxy.value <- append(oxy.value, NA)
    }
  }
  else {
    min.crit <- append(min.crit, NA)
    oxy.dep <- append(oxy.dep, NA)
    oxy.area <- append(oxy.area, NA)
    oxy.value <- append(oxy.value, NA)
  }
}

af_data <- data.frame('time' = oxy_data$time, 'depth' = oxy.dep, 'area' = oxy.area)

filter.af_data <- af_data %>% ##  AF FACTOR
  mutate(year = lubridate::year(time)) %>%
  mutate(month = lubridate::month(time)) %>%
  dplyr::filter(month > 4 & month < 10) %>%
  drop_na()
annual.filter.af_data <- filter.af_data %>%
  group_by(year) %>%
  dplyr::summarise(sum_area = sum(area)) %>%
  mutate(norm_sum_area = sum_area / max(nml$morphometry$A))
