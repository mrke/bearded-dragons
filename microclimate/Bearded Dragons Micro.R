
# R Implementation of an integration of the microclimate model of Warren Porter's Niche Mapper system 
# Michael Kearney November 2013

# This version uses the Australia Water Availability Project (AWAP) daily 5km climate
# layers for Australia for air temperature, relative humidity, rainfall and cloud cover
# and uses monthly soil moisture estimates (splined) and the Australia Soils database to
# obtain soil properties, including their change through time due to soil moisture.
# Humidity is only from 1971 onward. Cloud cover is only from 1990 onward (and is based
# on daily solar layers relative to clear sky calculations from NicheMapR).
# It also uses a global monthly soil moisture estimate from NOAA CPC Soil Moisture http://140.172.38.100/psd/thredds/catalog/Datasets/cpcsoil/catalog.html
# Aerosol attenuation can also be computed based on the Global Aerosol Data Set (GADS)
# Koepke, P., M. Hess, I. Schult, and E. P. Shettle. 1997. Global Aerosol Data Set. Max-Planck-Institut for Meteorologie, Hamburg
# by choosing the option 'rungads<-1' 

# required R packages
# raster
# sp
# ncdf
# XML
# dismo
# chron
# rgdal
# zoo
# RODBC

######################### model modes ###########################################################
mac<-0 # choose mac (1) or pc (0) 
writecsv<-0 # make Fortran code write output as csv files
write_input<-0 # write csv files of final input to working directory? 1=yes, 0=no.
runshade<-1 # run the model twice, once for each shade level (1) or just for the first shade level (0)?
runmoist<-1 # run soil moisture model (0=no, 1=yes)?
snowmodel<-0 # run the snow model (0=no, 1=yes)? - note that this runs slower
basic<-0 # for use with a simplified demo script 
shore<-0 # include tide effects (if 0, an empty matrix of tide effects is created)
rungads<-1 # use the Global Aerosol Database?
#########################################################################################################


############## location and climatic data  ###################################
spatial<-"c:/Australian Environment/" # place where climate input files are kept
sitemethod <- 1 # 0=specified single site long/lat, 1=place name search using geodis (needs internet)
longlat<-c(147.829971,-42.557127) #central plateau c(146.5666667,-41.85) Orford c(147.829971,-42.557127) 
loc <- "Walpeup, Victoria" # type in a location here, used if option 1 is chosen above
timezone<-0 # if timezone=1 (needs internet), uses GNtimezone function in package geonames to correct to local time zone (excluding daylight saving correction)
dailywind<-1 # use daily windspeed database?
terrain<-0 # include terrain (slope, aspect, horizon angles) (1) or not (0)?
soildata<-1 # include soil data for Australia (1) or not (0)?
ystart<-2013
yfinish<-2013
nyears<-yfinish-ystart+1# integer, number of years for which to run the microclimate model, only for AWAP data (!!max 10 years!!)

DEP <- c(0., 2.5, 5.,  10., 15., 20., 30., 50.,  100., 200.) # Soil nodes (cm) - keep spacing close near the surface, last value is where it is assumed that the soil temperature is at the annual mean air temperature


soilpro<-read.csv("c:/git/Tiliqua_rugosa/microclimate/soilprops.txt",header=FALSE)
colnames(soilpro)<-c('i','site','long','lat','desc','blkdens','clay','silt','sand')
soilpro<-subset(soilpro,site==1)
soilpro<-soilpro[,5:9]
soilpro[,2]<-1.4
soilpro[,3]<-90
soilpro[,4]<-0.0
soilpro[,5]<-10
#    
soil_depths<-c(2.5,7.5,22.5,45,80,150)
# plot(soilpro$clay~soil_depths,ylim=c(0,100),col='red',type='l')
# points(soilpro$sand~soil_depths,ylim=c(0,100),col='orange',type='l')
# points(soilpro$silt~soil_depths,ylim=c(0,100),col='grey',type='l')
# title(main=loc)
# legend("topleft", inset=.05,
#        legend=round(soilpro[1,3:5],1),bty="n", 
#        horiz=TRUE, bg=NULL, cex=0.8)

DEP2<-rep(0,18)
j<-1
for(i in 1:length(DEP2)){
  if(i%%2==0){
    DEP2[i]<-DEP2[i-1]+(DEP[j]-DEP2[i-1])/2
  }else{
    DEP2[i]<-DEP[j]
    j<-j+1
  }
}
DEP2<-as.data.frame(floor(DEP2))
colnames(DEP2)<-"DEPTH"
 
par(mfrow=c(2,2))


CampNormTbl9_1<-read.csv('../micro_australia/CampNormTbl9_1.csv')
Richter<-read.csv('../micro_australia/Richter_Table1_SI.csv')
dclay<-0.001 #mm
dsilt<-0.026 #mm
dsand<-1.05 #mm
a<-(soilpro$clay/100)*log(dclay) + (soilpro$sand/100)*log(dsand) + (soilpro$silt/100)*log(dsilt)
b.1<-(((soilpro$clay/100)*log(dclay)^2+(soilpro$sand/100)*log(dsand)^2+(soilpro$silt/100)*log(dsilt)^2)-a^2)^(1/2)
dg<-exp(a)
sigma_g<-exp(b.1)
PES<-(0.5*dg^(-1/2))*-1
b<--2*PES+0.2*sigma_g
PE<-PES*(soilpro$blkdens/1.3)^(0.67*b)
KS<-0.004*(1.3/soilpro$blkdens)^(1.3*b)*exp(-6.9*soilpro$clay/100-3.7*soilpro$silt/100)
BD<-soilpro$blkdens
   

# plot(KS~soil_depths,xlim=c(-1,200),ylim=c(0.000017,0.0058))
KS_spline <-spline(soil_depths,KS,n=201,xmin=0,xmax=200,method='natural')
#points(KS_spline$y,col='red',type='l')
KS_spline<-as.data.frame(cbind(KS_spline$x,KS_spline$y))
colnames(KS_spline)<-c('DEPTH','VALUE')
KS<-merge(DEP2,KS_spline)
KS<-c(KS[1,2],KS[,2])
KS[KS<0.000017]<-0.000017
   
#plot(PE~soil_depths,xlim=c(-1,200),ylim=c(-15,0))
PE_spline <-spline(soil_depths,PE,n=201,xmin=0,xmax=200,method='natural')
#points(PE_spline$y,col='red',type='l')
PE_spline<-as.data.frame(cbind(PE_spline$x,PE_spline$y))
colnames(PE_spline)<-c('DEPTH','VALUE')
PE<-merge(DEP2,PE_spline)
PE<-c(-1*PE[1,2],-1*PE[,2])
   
#plot(b~soil_depths,xlim=c(-1,200),ylim=c(2,24))
b_spline <-spline(soil_depths,b,n=201,xmin=0,xmax=200,method='natural')
#points(b_spline$y,col='red',type='l')
b_spline<-as.data.frame(cbind(b_spline$x,b_spline$y))
colnames(b_spline)<-c('DEPTH','VALUE')
b<-merge(DEP2,b_spline)
BB<-c(b[1,2],b[,2])
   
#plot(BD~soil_depths,xlim=c(-1,200),ylim=c(1,1.6))
BD_spline <-spline(soil_depths,BD,n=201,xmin=0,xmax=200,method='natural')
#points(BD_spline$y,col='red',type='l')
BD_spline<-as.data.frame(cbind(BD_spline$x,BD_spline$y))
colnames(BD_spline)<-c('DEPTH','VALUE')
BD<-merge(DEP2,BD_spline)
BD<-c(BD[1,2],BD[,2])

par(mfrow=c(1,1))

############# microclimate model parameters ################################
EC <- 0.0167238 # Eccenricity of the earth's orbit (current value 0.0167238, ranges between 0.0034 to 0.058)
RUF <- 0.004 # Roughness height (m), , e.g. sand is 0.05, grass may be 2.0, current allowed range: 0.001 (snow) - 2.0 cm.
# Next for parameters are segmented velocity profiles due to bushes, rocks etc. on the surface, IF NO EXPERIMENTAL WIND PROFILE DATA SET ALL THESE TO ZERO!
Z01 <- 0. # Top (1st) segment roughness height(m)
Z02 <- 0. # 2nd segment roughness height(m)
ZH1 <- 0. # Top of (1st) segment, height above surface(m)
ZH2 <- 0. # 2nd segment, height above surface(m)
SLE <- 0.96 # Substrate longwave IR emissivity (decimal %), typically close to 1
ERR <- 2.5 # Integrator error for soil temperature calculations
Thcond <- rep(2.5,10) # soil minerals thermal conductivity (W/mC)
Density <- rep(2650,10) # soil minerals density (kg/m3)
SpecHeat <- rep(870,10) # soil minerals specific heat (J/kg-K)
BulkDensity <- rep(1300,10) # soil bulk density (kg/m3)
cap<-1 # organic cap present on soil surface? (cap has lower conductivity - 0.2 W/mC - and higher specific heat 1920 J/kg-K)
SatWater <- rep(0.26,10) # volumetric water content at saturation (0.1 bar matric potential) (m3/m3)
Clay <- rep(20,10)#CLAY # clay content for matric potential calculations (%)
SoilMoist <- 0 # fractional soil moisture (decimal %)
rainmult<-1 # rain multiplier for surface soil moisture (use to induce runoff), proportion
SoilMoist_Init<-rep(0.2,10) # initial soil water content, m3/m3
evenrain<-1 # spread daily rainfall evenly across 24hrs (1) or one event at midnight (2)
maxpool<-500 # max depth for water pooling on the surface, mm (to account for runoff)
soiltype<-10
CampNormTbl9_1<-read.csv('../micro_australia/CampNormTbl9_1.csv')
fieldcap<-CampNormTbl9_1[soiltype,7] # field capacity, mm
wilting<-CampNormTbl9_1[soiltype,8]  # use value from digital atlas of Australian soils # wilting point, mm
# PE<-rep(CampNormTbl9_1[soiltype,4],19) #air entry potential J/kg 
# KS<-rep(CampNormTbl9_1[soiltype,6],19) #saturated conductivity, kg s/m3
# BB<-rep(CampNormTbl9_1[soiltype,5],19) #soil 'b' parameterPE<-rep(CampNormTbl9_1[soiltype,4],19) #air entry potential J/kg 
PE[1:9]<-CampNormTbl9_1[3,4] #air entry potential J/kg 
KS[1:9]<-CampNormTbl9_1[3,6] #saturated conductivity, kg s/m3
BB[1:9]<-CampNormTbl9_1[3,5] #soil 'b' parameter
PE[10:13]<-CampNormTbl9_1[4,4] #air entry potential J/kg 
KS[10:13]<-CampNormTbl9_1[4,6] #saturated conductivity, kg s/m3
BB[10:13]<-CampNormTbl9_1[4,5] #soil 'b' parameter
L<-c(0,0,8.18990859,7.991299442,7.796891252,7.420411664,7.059944542,6.385001059,5.768074989,4.816673431,4.0121088,1.833554792,0.946862989,0.635260544,0.804575,0.43525621,0.366052856,0,0)*10000
LAI<-0.0 # leaf area index, used to partition traspiration/evaporation from PET
REFL<-0.10 # soil reflectance (decimal %)
slope<-0. # slope (degrees, range 0-90)
aspect<-180. # aspect (degrees, 0 = North, range 0-360)
hori<-rep(0,24) # enter the horizon angles (degrees) so that they go from 0 degrees azimuth (north) clockwise in 15 degree intervals
PCTWET<-0. # percentage of surface area acting as a free water surface (%)
CMH2O <- 1. # precipitable cm H2O in air column, 0.1 = VERY DRY; 1.0 = MOIST AIR CONDITIONS; 2.0 = HUMID, TROPICAL CONDITIONS (note this is for the whole atmospheric profile, not just near the ground)  
TIMAXS <- c(1.0, 1.0, 0.0, 0.0)   # Time of Maximums for Air Wind RelHum Cloud (h), air & Wind max's relative to solar noon, humidity and cloud cover max's relative to sunrise            											
TIMINS <- c(0.0, 0.0, 1.0, 1.0)   # Time of Minimums for Air Wind RelHum Cloud (h), air & Wind min's relative to sunrise, humidity and cloud cover min's relative to solar noon
minshade<-0. # minimum available shade (%)
maxshade<-80. # maximum available shade (%)
manualshade<-1 # if using soildata, which includes shade, this will override the data from the database and force max shade to be the number specified above
Usrhyt <- 1# local height (cm) at which air temperature, relative humidity and wind speed calculatinos will be made 
rainwet<-1.5 # mm rain that causes soil to become 90% wet
snowtemp<-1.5 # temperature at which precipitation falls as snow (used for snow model)
snowdens<-0.325 # snow density (mg/m3)
densfun<-c(0.001369,0.1095) # slope and intercept of linear model of snow density as a function of day of year - if it is c(0,0) then fixed density used
snowmelt<-0.9 # proportion of calculated snowmelt that doesn't refreeze
undercatch<-1.0 # undercatch multipier for converting rainfall to snow
rainmelt<-0.0125#85 # paramter in equation that melts snow with rainfall as a function of air temp, start with 0.0125
warm<-0 # uniform warming of air temperature input to simulate climate change
loop<-0 # if doing multiple years, this shifts the starting year by the integer value


# run the model
maindir<-getwd()
setwd('/git/micro_australia/')
niche<-list(writecsv=writecsv,densfun=densfun,L=L,LAI=LAI,SoilMoist_Init=SoilMoist_Init,evenrain=evenrain,runmoist=runmoist,maxpool=maxpool,PE=PE,KS=KS,BB=BB,BD=BD,loop=loop,warm=warm,rainwet=rainwet,manualshade=manualshade,dailywind=dailywind,terrain=terrain,soildata=soildata,loc=loc,ystart=ystart,yfinish=yfinish,nyears=nyears,RUF=RUF,SLE=SLE,ERR=ERR,DEP=DEP,Thcond=Thcond,Density=Density,SpecHeat=SpecHeat,BulkDensity=BulkDensity,Clay=Clay,SatWater=SatWater,SoilMoist=SoilMoist,CMH2O=CMH2O,TIMAXS=TIMAXS,TIMINS=TIMINS,minshade=minshade,maxshade=maxshade,Usrhyt=Usrhyt,REFL=REFL,slope=slope,aspect=aspect,hori=hori,rungads=rungads,cap=cap,write_input=write_input,spatial=spatial,snowmodel=snowmodel,snowtemp=snowtemp,snowdens=snowdens,snowmelt=snowmelt,undercatch=undercatch,rainmelt=rainmelt,rainmult=rainmult,runshade=runshade)
source('NicheMapR_Setup_micro.R')
nicheout<-NicheMapR(niche)
setwd(maindir)

microdir<-'microclimate/'


# get output
dim<-nicheout$dim
metout<-as.data.frame(nicheout$metout[1:(dim*24),]) # above ground microclimatic conditions, min shade
shadmet<-as.data.frame(nicheout$shadmet[1:(dim*24),]) # above ground microclimatic conditions, max shade
soil<-as.data.frame(nicheout$soil[1:(dim*24),]) # soil temperatures, minimum shade
shadsoil<-as.data.frame(nicheout$shadsoil[1:(dim*24),]) # soil temperatures, maximum shade
soilmoist<-as.data.frame(nicheout$soilmoist[1:(dim*24),]) # soil water content, minimum shade
shadmoist<-as.data.frame(nicheout$shadmoist[1:(dim*24),]) # soil water content, maximum shade
humid<-as.data.frame(nicheout$humid[1:(dim*24),]) # soil humidity, minimum shade
shadhumid<-as.data.frame(nicheout$shadhumid[1:(dim*24),]) # soil humidity, maximum shade
soilpot<-as.data.frame(nicheout$soilpot[1:(dim*24),]) # soil water potential, minimum shade
shadpot<-as.data.frame(nicheout$shadpot[1:(dim*24),]) # soil water potential, maximum shade
rainfall<-as.data.frame(nicheout$RAINFALL)
MAXSHADES<-as.data.frame(nicheout$MAXSHADES)
elev<-as.numeric(nicheout$ALTT)
REFL<-as.numeric(nicheout$REFL)
longlat<-as.matrix(nicheout$longlat)
ectoin<-rbind(elev,REFL,longlat,0,0,ystart,ystart+nyears-1)


write.csv(metout,paste(microdir,'metout.csv',sep=""))
write.csv(soil,paste(microdir,'soil.csv',sep=""))
write.csv(soilpot,paste(microdir,'soilpot.csv',sep=""))
write.csv(humid,paste(microdir,'humid.csv',sep=""))
write.csv(soilmoist,paste(microdir,'soilmoist.csv',sep=""))
if(runshade==0){
  write.csv(metout,paste(microdir,'shadmet.csv',sep=""))
  write.csv(soil,paste(microdir,'shadsoil.csv',sep=""))
  write.csv(humid,paste(microdir,'shadhumid.csv',sep=""))
  write.csv(soilpot,paste(microdir,'shadpot.csv',sep=""))
  write.csv(soilmoist,paste(microdir,'shadmoist.csv',sep=""))
}else{
  write.csv(shadmet,paste(microdir,'shadmet.csv',sep=""))
  write.csv(shadsoil,paste(microdir,'shadsoil.csv',sep=""))
  write.csv(shadhumid,paste(microdir,'shadhumid.csv',sep=""))
  write.csv(shadpot,paste(microdir,'shadpot.csv',sep=""))
  write.csv(shadmoist,paste(microdir,'shadmoist.csv',sep=""))
}
write.csv(rainfall,paste(microdir,'rainfall.csv',sep=""))
write.csv(ectoin,paste(microdir,'ectoin.csv',sep=""))
write.csv(DEP,paste(microdir,'DEP.csv',sep=""))
write.csv(MAXSHADES,paste(microdir,'MAXSHADES.csv',sep=""))

if(!require(geonames)){
  stop('package "geonames" is required.')
}
tzone<-paste("Etc/GMT-10",sep="") # doing it this way ignores daylight savings!

dates<-seq(ISOdate(ystart,1,1,tz=tzone)-3600*12, ISOdate((ystart+nyears),1,1,tz=tzone)-3600*13, by="hours") 
dates<-subset(dates, format(dates, "%m/%d")!= "02/29") # remove leap years
dates<-subset(dates, !duplicated(as.matrix(dates[2110:2120])))
dates<-unique(dates)
metout<-cbind(dates,metout)
shadmet<-cbind(dates,shadmet)
soil<-cbind(dates,soil)
shadsoil<-cbind(dates,shadsoil)
soilmoist<-cbind(dates,soilmoist)
shadmoist<-cbind(dates,shadmoist)
humid<-cbind(dates,humid)
shadhumid<-cbind(dates,shadhumid)
soilpot<-cbind(dates,soilpot)
shadpot<-cbind(dates,shadpot)

dates2<-seq(ISOdate(ystart,1,1,tz=tzone)-3600*12, ISOdate((ystart+nyears),1,1,tz=tzone)-3600*13, by="days") 
dates2<-subset(dates2, format(dates2, "%m/%d")!= "02/29") # remove leap years
rainfall<-as.data.frame(cbind(dates2,rainfall))
colnames(rainfall)<-c('dates','rainfall')

plot(metout$SNOWDEP~dates,type='l')
points(shadmet$SNOWDEP~dates,type='l',col='blue')
#points(rainfall$rainfall~dates2,type='h',col='light blue')

dstart<-as.POSIXct(as.Date('01/01/2000', "%d/%m/%Y"))-3600*11
dfinish<-as.POSIXct(as.Date('31/12/2013', "%d/%m/%Y"))-3600*10
plotsoilmoist<-subset(soilmoist,  soilmoist$dates > dstart & soilmoist$dates < dfinish )
plotshadmoist<-subset(shadmoist,  shadmoist$dates > dstart & shadmoist$dates < dfinish )
plothumid<-subset(humid,  humid$dates > dstart & humid$dates < dfinish )
plotsoilpot<-subset(soilpot,  soilpot$dates > dstart & soilpot$dates < dfinish )
plotsoil<-subset(soil,  soil$dates > dstart & soil$dates < dfinish )
plotmetout<-subset(metout,  metout$dates > dstart & metout$dates < dfinish )

plot(plotsoilmoist$dates, plotsoilmoist[,4]*100,type='l',col = "red",lty=1,ylim = c(0,50),main=CampNormTbl9_1[soiltype,1],ylab='moisture (% vol)',xlab='date')
plot(plotsoilmoist$dates, plotsoilmoist[,5]*100,type='l',col = 3,lty=1,ylim = c(0,50),main=CampNormTbl9_1[soiltype,1],ylab='moisture (% vol)',xlab='date')
points(plotsoilmoist$dates, plotsoilmoist[,6]*100,type='l',col = 4,lty=1,ylim = c(0,50),main=CampNormTbl9_1[soiltype,1],ylab='moisture (% vol)',xlab='date')
points(plotsoilmoist$dates, plotsoilmoist[,7]*100,type='l',col = 5,lty=1,ylim = c(0,50),main=CampNormTbl9_1[soiltype,1],ylab='moisture (% vol)',xlab='date')
points(plotsoilmoist$dates, plotsoilmoist[,8]*100,type='l',col = 6,lty=1,ylim = c(0,50),main=CampNormTbl9_1[soiltype,1],ylab='moisture (% vol)',xlab='date')
points(plotsoilmoist$dates, plotsoilmoist[,9]*100,type='l',col = 7,lty=1,ylim = c(0,50),main=CampNormTbl9_1[soiltype,1],ylab='moisture (% vol)',xlab='date')
points(plotsoilmoist$dates, plotsoilmoist[,10]*100,type='l',col = 8,lty=1,ylim = c(0,50),main=CampNormTbl9_1[soiltype,1],ylab='moisture (% vol)',xlab='date')
points(plotsoilmoist$dates, plotsoilmoist[,11]*100,type='l',col = 9,lty=1,ylim = c(0,50),main=CampNormTbl9_1[soiltype,1],ylab='moisture (% vol)',xlab='date')
points(plotsoilmoist$dates, plotsoilmoist[,12]*100,type='l',col = 10,lty=1,ylim = c(0,50),main=CampNormTbl9_1[soiltype,1],ylab='moisture (% vol)',xlab='date')
points(plotsoilmoist$dates, plotsoilmoist[,13]*100,type='l',col = 11,lty=1,ylim = c(0,50),main=CampNormTbl9_1[soiltype,1],ylab='moisture (% vol)',xlab='date')
points(rainfall$rainfall~rainfall$dates,type='h',col='dark blue')
#abline(11,0)

plot(plotshadmoist$dates, plotshadmoist[,4]*100,type='l',col = "red",lty=1,ylim = c(0,50),main=CampNormTbl9_1[soiltype,1],ylab='moisture (% vol)',xlab='date')
plot(plotshadmoist$dates, plotshadmoist[,5]*100,type='l',col = 3,lty=1,ylim = c(0,50),main=CampNormTbl9_1[soiltype,1],ylab='moisture (% vol)',xlab='date')
points(plotshadmoist$dates, plotshadmoist[,6]*100,type='l',col = 4,lty=1,ylim = c(0,50),main=CampNormTbl9_1[soiltype,1],ylab='moisture (% vol)',xlab='date')
points(plotshadmoist$dates, plotshadmoist[,7]*100,type='l',col = 5,lty=1,ylim = c(0,50),main=CampNormTbl9_1[soiltype,1],ylab='moisture (% vol)',xlab='date')
points(plotshadmoist$dates, plotshadmoist[,8]*100,type='l',col = 6,lty=1,ylim = c(0,50),main=CampNormTbl9_1[soiltype,1],ylab='moisture (% vol)',xlab='date')
points(plotshadmoist$dates, plotshadmoist[,9]*100,type='l',col = 7,lty=1,ylim = c(0,50),main=CampNormTbl9_1[soiltype,1],ylab='moisture (% vol)',xlab='date')
points(plotshadmoist$dates, plotshadmoist[,10]*100,type='l',col = 8,lty=1,ylim = c(0,50),main=CampNormTbl9_1[soiltype,1],ylab='moisture (% vol)',xlab='date')
points(plotshadmoist$dates, plotshadmoist[,11]*100,type='l',col = 9,lty=1,ylim = c(0,50),main=CampNormTbl9_1[soiltype,1],ylab='moisture (% vol)',xlab='date')
points(plotshadmoist$dates, plotshadmoist[,12]*100,type='l',col = 10,lty=1,ylim = c(0,50),main=CampNormTbl9_1[soiltype,1],ylab='moisture (% vol)',xlab='date')
points(plotshadmoist$dates, plotshadmoist[,13]*100,type='l',col = 11,lty=1,ylim = c(0,50),main=CampNormTbl9_1[soiltype,1],ylab='moisture (% vol)',xlab='date')
points(rainfall$rainfall~rainfall$dates,type='h',col='dark blue')
abline(11,0)

plot(plothumid$dates, plothumid[,4]*100,type='l',col = "red",lty=1,ylim = c(0,100),main=CampNormTbl9_1[soiltype,1],ylab='relative humdity (%)',xlab='date')
points(plothumid$dates, plothumid[,5]*100,type='l',col = 3,lty=1,ylim = c(0,50),main=CampNormTbl9_1[soiltype,1],ylab='relative humdity (%)',xlab='date')
points(plothumid$dates, plothumid[,6]*100,type='l',col = 4,lty=1,ylim = c(0,50),main=CampNormTbl9_1[soiltype,1],ylab='relative humdity (%)',xlab='date')
points(plothumid$dates, plothumid[,7]*100,type='l',col = 5,lty=1,ylim = c(0,50),main=CampNormTbl9_1[soiltype,1],ylab='relative humdity (%)',xlab='date')
points(plothumid$dates, plothumid[,8]*100,type='l',col = 6,lty=1,ylim = c(0,50),main=CampNormTbl9_1[soiltype,1],ylab='relative humdity (%)',xlab='date')
points(plothumid$dates, plothumid[,9]*100,type='l',col = 7,lty=1,ylim = c(0,50),main=CampNormTbl9_1[soiltype,1],ylab='relative humdity (%)',xlab='date')
points(plothumid$dates, plothumid[,10]*100,type='l',col = 8,lty=1,ylim = c(0,50),main=CampNormTbl9_1[soiltype,1],ylab='relative humdity (%)',xlab='date')
points(plothumid$dates, plothumid[,11]*100,type='l',col = 9,lty=1,ylim = c(0,50),main=CampNormTbl9_1[soiltype,1],ylab='relative humdity (%)',xlab='date')
points(plothumid$dates, plothumid[,12]*100,type='l',col = 10,lty=1,ylim = c(0,50),main=CampNormTbl9_1[soiltype,1],ylab='relative humdity (%)',xlab='date')
points(plothumid$dates, plothumid[,13]*100,type='l',col = 11,lty=1,ylim = c(0,50),main=CampNormTbl9_1[soiltype,1],ylab='relative humdity (%)',xlab='date')

plot(plotsoilpot$dates, plotsoilpot[,4],type='l',col = "red",lty=1,ylim = c(-5000,0),main=CampNormTbl9_1[soiltype,1],ylab='water potential (J/kg)',xlab='date')
points(plotsoilpot$dates, plotsoilpot[,5],type='l',col = 3,lty=1,ylim = c(-300,0),main=CampNormTbl9_1[soiltype,1],ylab='water potential (J/kg)',xlab='date')
points(plotsoilpot$dates, plotsoilpot[,6],type='l',col = 4,lty=1,ylim = c(0,50),main=CampNormTbl9_1[soiltype,1],ylab='relative humdity (%)',xlab='date')
points(plotsoilpot$dates, plotsoilpot[,7],type='l',col = 5,lty=1,ylim = c(0,50),main=CampNormTbl9_1[soiltype,1],ylab='relative humdity (%)',xlab='date')
points(plotsoilpot$dates, plotsoilpot[,8],type='l',col = 6,lty=1,ylim = c(0,50),main=CampNormTbl9_1[soiltype,1],ylab='relative humdity (%)',xlab='date')
points(plotsoilpot$dates, plotsoilpot[,9],type='l',col = 7,lty=1,ylim = c(0,50),main=CampNormTbl9_1[soiltype,1],ylab='relative humdity (%)',xlab='date')
points(plotsoilpot$dates, plotsoilpot[,10],type='l',col = 8,lty=1,ylim = c(0,50),main=CampNormTbl9_1[soiltype,1],ylab='relative humdity (%)',xlab='date')
points(plotsoilpot$dates, plotsoilpot[,11],type='l',col = 9,lty=1,ylim = c(0,50),main=CampNormTbl9_1[soiltype,1],ylab='relative humdity (%)',xlab='date')
points(plotsoilpot$dates, plotsoilpot[,12],type='l',col = 10,lty=1,ylim = c(0,50),main=CampNormTbl9_1[soiltype,1],ylab='relative humdity (%)',xlab='date')
points(plotsoilpot$dates, plotsoilpot[,13],type='l',col = 11,lty=1,ylim = c(0,50),main=CampNormTbl9_1[soiltype,1],ylab='relative humdity (%)',xlab='date')

plot(plotsoil$dates, plotsoil[,4],type='l',col = "red",lty=1,ylim = c(-20,80),main=CampNormTbl9_1[soiltype,1],ylab='temperature (C)',xlab='date')
points(plotsoil$dates, plotsoil[,5],type='l',col = 3,lty=1,ylim = c(0,50),main=CampNormTbl9_1[soiltype,1],ylab='relative humdity (%)',xlab='date')
points(plotsoil$dates, plotsoil[,6],type='l',col = 4,lty=1,ylim = c(0,50),main=CampNormTbl9_1[soiltype,1],ylab='relative humdity (%)',xlab='date')
points(plotsoil$dates, plotsoil[,7],type='l',col = 5,lty=1,ylim = c(0,50),main=CampNormTbl9_1[soiltype,1],ylab='relative humdity (%)',xlab='date')
points(plotsoil$dates, plotsoil[,8],type='l',col = 6,lty=1,ylim = c(0,50),main=CampNormTbl9_1[soiltype,1],ylab='relative humdity (%)',xlab='date')
points(plotsoil$dates, plotsoil[,9],type='l',col = 7,lty=1,ylim = c(0,50),main=CampNormTbl9_1[soiltype,1],ylab='relative humdity (%)',xlab='date')
points(plotsoil$dates, plotsoil[,10],type='l',col = 8,lty=1,ylim = c(0,50),main=CampNormTbl9_1[soiltype,1],ylab='relative humdity (%)',xlab='date')
points(plotsoil$dates, plotsoil[,11],type='l',col = 9,lty=1,ylim = c(0,50),main=CampNormTbl9_1[soiltype,1],ylab='relative humdity (%)',xlab='date')
points(plotsoil$dates, plotsoil[,12],type='l',col = 10,lty=1,ylim = c(0,50),main=CampNormTbl9_1[soiltype,1],ylab='relative humdity (%)',xlab='date')
points(plotsoil$dates, plotsoil[,13],type='l',col = 11,lty=1,ylim = c(0,50),main=CampNormTbl9_1[soiltype,1],ylab='relative humdity (%)',xlab='date')



plot(soilmoist$WC5cm~soilmoist$dates,type='l')