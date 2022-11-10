#written by CCC on 16 July 2018 to develop SSS inflow file for FCR GLM
#updated 1 June 2020 to become "tidy" and update inflow nutrient fractions for 2013-2019

setwd("../inputs")

#load packages
library(zoo)
library(tidyverse)
library(lubridate)
library(magrittr)

#creating new dataframe with list of all dates
datelist<-seq.Date(as.Date("2015/07/07"),as.Date("2021/12/31"), "days") #changed from May 15, 2013 because of NA in flow
datelist<-as.data.frame(datelist)
colnames(datelist)=c("time")
datelist$time<-as.POSIXct(strptime(datelist$time, "%Y-%m-%d", tz="EST"))

#now let's merge with chemistry
#first pull in FCR chem data from 2013-2021 from EDI
#inUrl1  <- "https://pasta.lternet.edu/package/data/eml/edi/199/10/aa2ccc23688fc908f9d61cb217210a3d"
#infile1 <- paste0(getwd(),"/chemistry_2013_2021.csv")
#download.file(inUrl1,infile1,method="curl",extra='-k')

FCRchem <- read.csv("chemistry_2013_2021.csv", header=T) %>%
  select(Reservoir:DIC_mgL) %>%
  filter(Reservoir=="FCR") %>%
  filter(Site==50) %>%
  filter(Depth_m==8) %>%
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d", tz="EST"))) %>%
  rename(time = DateTime) %>%
  select(time:DIC_mgL)

alldata<-merge(datelist, FCRchem, by="time", all.x=TRUE)

#pull in SSS operations file 
inflowoxy <- read.csv("Calc_HOX_flow_DO_20211102.csv", header=T) %>%
  #within this file, the eductor nozzles increased flow rate by a factor of 4, so 227 LPM = 1135 LPM
  select(time,SSS_m3.day,mmol.O2.m3.day) %>%
  mutate(SSS_m3.day = SSS_m3.day * (1/86400))  %>%
  #convert m3/day to m3/second as needed by GLM
  mutate(SALT = rep(0,length(inflowoxy$time))) %>%
  #add salt column as needed by GLM
  mutate(time = as.POSIXct(strptime(time, "%Y-%m-%d", tz="EST"))) %>%
  rename(FLOW = SSS_m3.day, OXY_oxy=mmol.O2.m3.day)

#some diagnostic plots
plot(inflowoxy$time, inflowoxy$FLOW, type = "o")
plot(inflowoxy$time, inflowoxy$OXY_oxy, type = "l", col = "red", ylab="mmol O2/m3/d added per day",
     xlab="Time")

alldata<-merge(alldata, inflowoxy, by="time", all.x=TRUE)

#read in lab dataset of dissolved silica, measured in summer 2014 only
#inUrl1  <- "https://pasta.lternet.edu/package/data/eml/edi/542/1/791ec9ca0f1cb9361fa6a03fae8dfc95" 
#infile1 <- paste0(getwd(),"/silica_master_df.csv")
#download.file(inUrl1,infile1,method="curl", extra='-k')

silica <- read.csv("silica_master_df.csv", header=T) %>%
  dplyr::filter(Reservoir == "FCR") %>% 
  dplyr::filter(Site == 50) %>% #100 = weir inflow site
  dplyr::filter(Depth_m == 9) %>% #because no Si was measured at 8m
  select(DateTime, DRSI_mgL) %>%
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d", tz="EST"))) %>%
  rename(time = DateTime)

#diagnostic plot of silica
plot(silica$time, silica$DRSI_mgL)
hist(silica$DRSI_mgL)
median(silica$DRSI_mgL) #this median concentration is going to be used to set as 
#the constant Si inflow conc in both wetland & weir inflows

#read in dataset of CH4 from EDI
#inUrl1  <- "https://pasta.lternet.edu/package/data/eml/edi/551/6/38d72673295864956cccd6bbba99a1a3"
#infile1 <- paste0(getwd(),"/final_GHG_2015-2021.csv")
#download.file(inUrl1,infile1,method="curl")

ghg <- read.csv("final_GHG_2015-2021.csv", header=T) %>%
  dplyr::filter(Reservoir == "FCR") %>%
  dplyr::filter(Site == 50) %>% #weir inflow
  dplyr::filter(Depth_m == 8) %>% 
  select(DateTime, ch4_umolL) %>%
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d", tz="EST"))) %>%
  rename(time = DateTime, CAR_ch4 = ch4_umolL) %>%
  group_by(time) %>%
  drop_na %>% 
  summarise(CAR_ch4 = mean(CAR_ch4)) %>%
  dplyr::filter(CAR_ch4<0.19) #remove outliers
plot(ghg$time, ghg$CAR_ch4)

datelist2<-seq.Date(as.Date(first(ghg$time)),as.Date(last(ghg$time)), "days")
datelist2<-as.data.frame(datelist2)
colnames(datelist2)=c("time")
datelist2$time<-as.POSIXct(strptime(datelist2$time, "%Y-%m-%d", tz="EST"))

ghg1 <- merge(datelist2, ghg, by="time", all.x=TRUE) 
ghg1$CAR_ch4 <- na.fill(na.approx(ghg1$CAR_ch4), "extend")
plot(ghg1$time, ghg1$CAR_ch4) 

alldata<-merge(alldata, ghg1, by="time", all.x=TRUE) %>% 
  group_by(time) %>% 
  summarise_all(mean, na.RM=TRUE)

#read in lab dataset of pH at inflow
#inUrl1  <- "https://pasta.lternet.edu/package/data/eml/edi/198/10/b3bd353312f9e37ca392e2a5315cc9da" 
#infile1 <- paste0(getwd(),"/YSI_PAR_profiles_2013-2021.csv")
#download.file(inUrl1,infile1,method="curl", extra='-k')

pH <- read.csv("YSI_PAR_profiles_2013-2021.csv", header=T) %>%
  dplyr::filter(Reservoir == "FCR") %>% 
  dplyr::filter(Site == 50) %>% #50 = deephole
  dplyr::filter(Depth_m==8) %>% 
  select(DateTime, Depth_m, pH) %>%
  drop_na() %>% 
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d", tz="EST"))) %>%
  rename(time = DateTime,
         CAR_pH = pH)

#diagnostic plot of pH
plot(pH$time, pH$CAR_pH)
hist(pH$CAR_pH)
median(pH$CAR_pH) #this median concentration is going to be used to set up the SSS file

#need to merge the inflow OXY file with 8m chemistry to create SSS inflows

#some cool long-term plots of 8m chemistry
plot(alldata$time, alldata$SRP_ugL)
plot(alldata$time, alldata$DOC_mgL)
plot(alldata$time, alldata$NO3NO2_ugL)
plot(alldata$time, alldata$NH4_ugL)
plot(alldata$time, alldata$TN_ugL)
plot(alldata$time, alldata$TP_ugL)
plot(alldata$time, alldata$DIC_mgL)

alldata$time[which(duplicated(alldata$time))]

lastrow <- length(alldata$time) #need for extend function below
#now need to interpolate missing values in chem; setting 1st and last value in time series as medians
#then linearly interpolating the middle missing values
alldata$TN_ugL[1]<-median(na.exclude(alldata$TN_ugL))
alldata$TN_ugL[lastrow]<-median(na.exclude(alldata$TN_ugL)) 
alldata$TN_ugL<-na.fill(na.approx(alldata$TN_ugL),"extend")

alldata$TP_ugL[1]<-median(na.exclude(alldata$TP_ugL))
alldata$TP_ugL[lastrow]<-median(na.exclude(alldata$TP_ugL))
alldata$TP_ugL<-na.fill(na.approx(alldata$TP_ugL),"extend")

alldata$NH4_ugL[1]<-median(na.exclude(alldata$NH4_ugL))
alldata$NH4_ugL[lastrow]<-median(na.exclude(alldata$NH4_ugL))
alldata$NH4_ugL<-na.fill(na.approx(alldata$NH4_ugL),"extend")

alldata$NO3NO2_ugL[1]<-median(na.exclude(alldata$NO3NO2_ugL))
alldata$NO3NO2_ugL[lastrow]<-median(na.exclude(alldata$NO3NO2_ugL))
alldata$NO3NO2_ugL<-na.fill(na.approx(alldata$NO3NO2_ugL),"extend")

alldata$SRP_ugL[1]<-median(na.exclude(alldata$SRP_ugL))
alldata$SRP_ugL[lastrow]<-median(na.exclude(alldata$SRP_ugL))
alldata$SRP_ugL<-na.fill(na.approx(alldata$SRP_ugL),"extend")

alldata$DOC_mgL[1]<-median(na.exclude(alldata$DOC_mgL))
alldata$DOC_mgL[lastrow]<-median(na.exclude(alldata$DOC_mgL))
alldata$DOC_mgL<-na.fill(na.approx(alldata$DOC_mgL),"extend")

alldata$DIC_mgL[1]<-median(na.exclude(alldata$DIC_mgL))
alldata$DIC_mgL[lastrow]<-median(na.exclude(alldata$DIC_mgL))
alldata$DIC_mgL<-na.fill(na.approx(alldata$DIC_mgL),"extend")

alldata <- alldata[(!duplicated(alldata$time)),]#remove duplicated dates

#need to convert mass observed data into mmol/m3 units for two DOC pools
SSS_inflow <- alldata %>% 
  mutate(NIT_amm = NH4_ugL*1000*0.001*(1/18.04)) %>% 
  mutate(NIT_nit = NO3NO2_ugL*1000*0.001*(1/62.00)) %>% #as all NO2 is converted to NO3
  mutate(PHS_frp = SRP_ugL*1000*0.001*(1/94.9714)) %>% 
  mutate(OGM_doc = DOC_mgL*1000*(1/12.01)* 0.10) %>% #assuming 10% of total DOC is in labile DOC pool (Wetzel page 753)
  mutate(OGM_docr = 1.5*DOC_mgL*1000*(1/12.01)* 0.90) %>% #assuming 90% of total DOC is in recalcitrant DOC pool
  mutate(TN_ugL = TN_ugL*1000*0.001*(1/14)) %>% 
  mutate(TP_ugL = TP_ugL*1000*0.001*(1/30.97)) %>% 
  mutate(OGM_poc = 0.1*(OGM_doc+OGM_docr)) %>% #assuming that 10% of DOC is POC (Wetzel page 755)
  mutate(OGM_don = (5/6)*(TN_ugL-(NIT_amm+NIT_nit))*0.10) %>% #DON is ~5x greater than PON (Wetzel page 220)
  mutate(OGM_donr = (5/6)*(TN_ugL-(NIT_amm+NIT_nit))*0.90) %>% #to keep mass balance with DOC, DONr is 90% of total DON
  mutate(OGM_pon = (1/6)*(TN_ugL-(NIT_amm+NIT_nit))) %>% #detemined by subtraction
  mutate(OGM_dop = 0.3*(TP_ugL-PHS_frp)*0.10) %>% #Wetzel page 241, 70% of total organic P = particulate organic; 30% = dissolved organic P
  mutate(OGM_dopr = 0.3*(TP_ugL-PHS_frp)*0.90) %>% #to keep mass balance with DOC & DON, DOPr is 90% of total DOP
  mutate(OGM_pop = TP_ugL-(OGM_dop+OGM_dopr+PHS_frp)) %>% # #In lieu of having the adsorbed P pool activated in the model, need to have higher complexed P
  mutate(CAR_dic = DIC_mgL*1000*(1/52.515)) %>% #Long-term median pH of FCR is 6.5, at which point CO2/HCO3 is about 50-50
#given this disparity, using a 50-50 weighted molecular weight (44.01 g/mol and 61.02 g/mol, respectively)
  mutate(CAR_pH = median(pH$CAR_pH)) %>% 
  mutate(SIL_rsi = median(pH$CAR_pH)*1000*(1/60.08)) %>%  #setting the Silica concentration to the median 2014 9m concentration for consistency
  select(time, FLOW, SALT, OXY_oxy, NIT_amm:SIL_rsi, CAR_ch4)

#reality check of mass balance: these histograms should be at zero minus rounding errors
#hist(SSS_inflow$TP_ugL - (SSS_inflow$PHS_frp + SSS_inflow$OGM_dop + SSS_inflow$OGM_pop))
#hist(SSS_inflow$TN_ugL - (SSS_inflow$NIT_amm + SSS_inflow$NIT_nit + SSS_inflow$OGM_don + SSS_inflow$OGM_pon))

SSS_inflow1 <- SSS_inflow %>%
  mutate(TRC_tr1 = rep(1,length(SSS_inflow$time)),
         TRC_age = rep(0,length(SSS_inflow$time)),
         NCS_ss1 = rep(1,length(SSS_inflow$time)),
         NCS_ss2 = rep(1,length(SSS_inflow$time)),
         CAR_ch4_bub = rep(0,length(SSS_inflow$time)),
         PHS_frp_ads = PHS_frp*0.2, #from Matt pers comm 11 April 2022
         OGM_cpom = OGM_poc, #guessing for now
         PHY_cyano = rep(0,length(SSS_inflow$time)),
         PHY_cyano_IN = rep(0,length(SSS_inflow$time)),
         PHY_cyano_IP = rep(0,length(SSS_inflow$time)),
         PHY_green = rep(0,length(SSS_inflow$time)),
         PHY_green_IN = rep(0,length(SSS_inflow$time)),
         PHY_green_IP = rep(0,length(SSS_inflow$time)),
         PHY_diatom = rep(0,length(SSS_inflow$time)),
         PHY_diatom_IN = rep(0,length(SSS_inflow$time)),
         PHY_diatom_IP = rep(0,length(SSS_inflow$time)),
         BIV_filtfrac = rep(0,length(SSS_inflow$time))) %>% 
  select(time, FLOW, SALT, 
         TRC_tr1, TRC_age, 
         NCS_ss1, 
         NCS_ss2, 
         OXY_oxy, CAR_dic, 
         CAR_pH, CAR_ch4, 
         CAR_ch4_bub, 
         SIL_rsi, NIT_amm, NIT_nit, PHS_frp, 
         PHS_frp_ads,
         OGM_doc, OGM_poc, OGM_don, OGM_pon, OGM_dop, OGM_pop, OGM_docr, OGM_donr, 
         OGM_dopr, 
         OGM_cpom, 
         PHY_cyano, PHY_cyano_IN, PHY_cyano_IP, PHY_green, 
         PHY_green_IN, PHY_green_IP, PHY_diatom, PHY_diatom_IN, PHY_diatom_IP,
         BIV_filtfrac) %>% 
  mutate_if(is.numeric, round, 4) #round to 4 digits 

#now, need to set water temperature of this file to CTD observations at 8 m, the depth
# the HOx injects water into the hypolimnion
#inUrl1  <- "https://pasta.lternet.edu/package/data/eml/edi/200/12/0a62d1946e8d9a511bc1404e69e59b8c" 
#infile1 <- paste0(getwd(),"/CTD_2013_2021.csv")
#download.file(inUrl1,infile1,method="curl")

CTD<-read.csv("CTD_2013_2021.csv", header=TRUE) #now need to get temp at 8m for inflow
CTD8 <- CTD %>%
  select(Reservoir:Temp_C) %>%
  dplyr::filter(Reservoir=="FCR") %>%
  dplyr::filter(Site==50) %>%
  rename(time=Date, TEMP=Temp_C) %>%
  mutate(time = as.POSIXct(strptime(time, "%Y-%m-%d", tz="EST"))) %>%
  mutate(Depth_m = round(Depth_m, digits=0)) %>%
  group_by(time) %>% 
  dplyr::filter(Depth_m==8) %>%
  summarise(TEMP=mean(TEMP)) %>%
  dplyr::filter(TEMP<17) #remove outlier from 2014

#diagnostic plot to check for 8m water temp
plot(CTD8$time, CTD8$TEMP, type = "o")

#make final SSS inflow file, setting 12/31 of each year to 4oC in lieu of CTD data for interpolation
SSS_inflowALL<-merge(SSS_inflow1,CTD8, by="time",all.x=TRUE)
SSS_inflowALL$TEMP[1]<-12.5
SSS_inflowALL$TEMP[178]<-4
SSS_inflowALL$TEMP[544]<-4
SSS_inflowALL$TEMP[909]<-4
SSS_inflowALL$TEMP[1274]<-4
SSS_inflowALL$TEMP[1639]<-4
SSS_inflowALL$TEMP[2005]<-4
SSS_inflowALL$TEMP[2370]<-4 #set last row as 4oC in prep for freezing
SSS_inflowALL$TEMP<-na.fill(na.approx(SSS_inflowALL$TEMP),"extend")
#SSS_inflowALL$CAR_ch4 <-na.fill(na.approx(SSS_inflowALL$CAR_ch4), "extend")
plot(SSS_inflowALL$time, SSS_inflowALL$TEMP, type = "o")

SSS_inflowALL<-SSS_inflowALL %>%
  select(time, FLOW, TEMP, SALT:BIV_filtfrac)  %>% #get all of the columns in order
  mutate_if(is.numeric, round, 4) #round to 4 digits 

SSS_inflowALL[which(duplicated(SSS_inflowALL$time)),] #identify if there are repeated dates

SSS_inflowALL <- SSS_inflowALL[(!duplicated(SSS_inflowALL$time)),] #remove repeated dates

#et voila! the final inflow file for the SSS for 2 pools of DOC
write.csv(SSS_inflowALL, "FCR_SSS_inflow_2015_2022_20221109_allfractions_2DOCpools.csv", row.names = FALSE)

# #NOW NEED TO MAKE file for 1 pool of DOC
# #need to convert mass observed data into mmol/m3 units for two DOC pools
# SSS_inflow <- alldata %>% 
#   select(time, FLOW, SALT, TN_ugL:CAR_ch4, OXY_oxy) %>%
#   mutate(NIT_amm = NH4_ugL*1000*0.001*(1/18.04)) %>% 
#   mutate(NIT_nit = NO3NO2_ugL*1000*0.001*(1/62.00)) %>% #as all NO2 is converted to NO3
#   mutate(PHS_frp = SRP_ugL*1000*0.001*(1/94.9714)) %>% 
#   mutate(OGM_doc = DOC_mgL*1000*(1/12.01)) %>% 
#   mutate(TN_ugL = TN_ugL*1000*0.001*(1/14)) %>% 
#   mutate(TP_ugL = TP_ugL*1000*0.001*(1/30.97)) %>% 
#   mutate(DRSI_mgL = rep(median(silica$DRSI_mgL),length(SSS_inflow$time))) %>% #setting inflow to median 9m concentration from 2014
#   mutate(OGM_poc = 0.1*(OGM_doc)) %>% #assuming that 10% of DOC is POC (Wetzel page 755)
#   mutate(OGM_don = (5/6)*(TN_ugL-(NIT_amm+NIT_nit))) %>% #DON is ~5x greater than PON (Wetzel page 220)
#   mutate(OGM_pon = (1/6)*(TN_ugL-(NIT_amm+NIT_nit))) %>%
#   mutate(OGM_dop = 0.3*(TP_ugL-PHS_frp)) %>% #Wetzel page 241, 70% of total organic P = particulate organic; 30% = dissolved organic P
#   mutate(OGM_pop = 0.7*(TP_ugL-PHS_frp)) %>% 
#   mutate(PHS_frp_ads = PHS_frp) %>% #Following Farrell et al. 2020 EcolMod
#   mutate(CAR_dic = DIC_mgL*1000*(1/52.515)) %>% #Long-term avg pH of FCR is 6.5, at which point CO2/HCO3 is about 50-50
#   #given this disparity, using a 50-50 weighted molecular weight (44.01 g/mol and 61.02 g/mol, respectively)
#   mutate(SIL_rsi = DRSI_mgL*1000*(1/60.08))  #setting the Silica concentration to the median 2014 inflow concentration for consistency
# 
# #reality check of mass balance: these histograms should be at zero minus rounding errors
# hist(SSS_inflow$TP_ugL - (SSS_inflow$PHS_frp + SSS_inflow$OGM_dop + SSS_inflow$OGM_pop))
# hist(SSS_inflow$TN_ugL - (SSS_inflow$NIT_amm + SSS_inflow$NIT_nit + SSS_inflow$OGM_don + SSS_inflow$OGM_pon))
# 
# SSS_inflow <- SSS_inflow %>%
#   select(time:SALT, OXY_oxy:CAR_dic, CAR_ch4, SIL_rsi)
# 
# #now need to get water temperature at 8m to set as the inflow SSS temp
# #need to pull in from EDI
# inUrl1  <- "https://pasta.lternet.edu/package/data/eml/edi/200/10/2461524a7da8f1906bfc3806d594f94c" 
# infile1 <- paste0(getwd(),"/CTD_final_2013_2019.csv")
# download.file(inUrl1,infile1,method="curl")
# 
# CTD<-read.csv("CTD_final_2013_2019.csv", header=TRUE) #now need to get temp at 8m for inflow
# CTD8 <- CTD %>%
#   select(Reservoir:Temp_C) %>%
#   filter(Reservoir=="FCR") %>%
#   filter(Site==50) %>%
#   rename(time=Date, TEMP=Temp_C) %>%
#   mutate(time = as.POSIXct(strptime(time, "%Y-%m-%d", tz="EST"))) %>%
#   mutate(Depth_m = round(Depth_m, digits=0)) %>%
#   group_by(time) %>% 
#   filter(Depth_m==8) %>%
#   summarise(TEMP=mean(TEMP)) %>%
#   filter(TEMP<17) #remove crazy outlier from 2014
# 
# #diagnostic plot to check for 8m water temp
# plot(CTD8$time, CTD8$TEMP, type = "o")
# 
# #make final SSS inflow file, setting 12/31 of each year to 4oC in lieu of CTD data for interpolation
# SSS_inflowALL<-merge(SSS_inflow,CTD8, by="time",all.x=TRUE)
# SSS_inflowALL$TEMP[231]<-4
# SSS_inflowALL$TEMP[596]<-4
# SSS_inflowALL$TEMP[961]<-4
# SSS_inflowALL$TEMP[1329]<-4
# SSS_inflowALL$TEMP[1694]<-4
# SSS_inflowALL$TEMP[2059]<-4
# SSS_inflowALL$TEMP[2424]<-4 #set last row as 4oC in prep for freezing
# SSS_inflowALL$TEMP<-na.fill(na.approx(SSS_inflowALL$TEMP),"extend")
# #SSS_inflowALL$CAR_ch4 <-na.fill(na.approx(SSS_inflowALL$CAR_ch4), "extend")
# plot(SSS_inflowALL$time, SSS_inflowALL$TEMP, type = "o")
# 
# #need to fill in missing data for earlier part of 2013-2015 CH4 data
# #in lieu of having data, setting concentrations to mean observations during 2015-2019
# #first thing first, need to figure out what mean concentration is on May 15 (day 135)
# missingdata <- ghg1 %>% 
#   mutate(DOY = yday(time)) %>%
#   group_by(DOY) %>%
#   summarise(CAR_ch4=mean(CAR_ch4))
# plot(missingdata$DOY, missingdata$CAR_ch4)
# 
# #setting missing 2013-2014 data as mean values observed at 8m during 2015-2019; 
# #I recognize this isn't perfect but AR1 & mechanistic models don't do a good job here either :(
# for(i in 1:length(SSS_inflowALL$time)){
#   if(is.na(SSS_inflowALL$CAR_ch4[i])){
#     SSS_inflowALL$CAR_ch4[i] = missingdata$CAR_ch4[which(missingdata$DOY==yday(SSS_inflowALL$time[i]))]
#   }
# }  
# 
# plot(SSS_inflowALL$time, SSS_inflowALL$CAR_ch4) #not great, but not horrible and keeps pattern
# 
# SSS_inflowALL<-SSS_inflowALL %>%
#   select(time, FLOW, TEMP, SALT:OGM_doc, OGM_poc:SIL_rsi)  %>% #get all of the columns in order
#   mutate_if(is.numeric, round, 4) #round to 4 digits 
# 
# SSS_inflowALL[which(duplicated(SSS_inflowALL$time)),] #identify if there are repeated dates
# 
# SSS_inflowALL <- SSS_inflowALL[(!duplicated(SSS_inflowALL$time)),] #remove repeated dates
# 
# #et voila! the final inflow file for 2 DOC pools
# write.csv(SSS_inflowALL, "FCR_SSS_inflow_2013_2019_20200607_allfractions_1DOCpool.csv", row.names = FALSE)
# 

##################################################################
#now make SSS outflow file
SSSoutflow<-SSS_inflowALL[,c(1,2)]
write.csv(SSSoutflow, "FCR_SSS_outflow_2013_2019_20200601.csv", row.names = FALSE)
