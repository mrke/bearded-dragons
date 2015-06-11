basedir<-getwd()
workdir<-"/vlsci/VR0212/shared/NicheMapR_Working/projects/sleepy/"

args <- (commandArgs(TRUE))
simnum<-as.numeric(args[1])
bioregion<-as.numeric(args[2])

barcoo<-paste('/scratch/VR0212/bio',bioregion,'/',sep="")
load(paste(barcoo,'longlat.bin',sep=''))
longlat <- data[simnum,1:2]

numsites<-ceiling(nrow(data)/2/1000)
jstart<-numsites*(simnum-1)+1
jfinish<-numsites*(simnum-1)+numsites

if(jstart<=nrow(data)){
  
  if(jfinish>nrow(data)){
    jfinish<-nrow(data)
  }
  
  
  for(jobnum in jstart:jfinish){
    
    longlat <- data[jobnum,1:2]

#longlat<-c(139.3109, -33.888)
#longlats<-read.csv("/vlsci/VR0212/shared/NicheMapR_Working/bioregion_points_final.csv")
#longlats<-subset(longlats,RASTERVALU!=0)
#lng<-which.min(abs(longlats[,1]-longlat[1]))
#coarselong<-longlats[lng,1]
#longlats<-subset(longlats,longlats$long==longlats[lng,1])
#lt<-which.min(abs(longlats[,2]-longlat[2]))
#coarselat<-longlats[lt,2]
#bioregion<-longlats[lt,3]
#barcoo<-paste('/scratch/VR0212/bio',bioregion,'/',sep="")
#load(paste(barcoo,'longlat.bin',sep=''))
#jobnum<-as.numeric(rownames(subset(data,data$V2==coarselong & data$V3==coarselat)))
#quadrangle<-jobnum


spatial<-"c:/Australian Environment/" # place where climate input files are kept
setwd("/vlsci/VR0212/shared/NicheMapR_Working/")
############## location and climatic data  ###################################
sitemethod <- 0 # 0=specified single site long/lat, 1=place name search using geodis (needs internet)
#longlat<-c(139.3109, -33.888) #c(139.3109, -33.888) #Mt Mary site#c(139.35, -33.93)<- Kerr and Bull 2004 site #
loc <- "Broken Hill, Australia" # type in a location here, used if option 1 is chosen above
timezone<-0 # if timezone=1 (needs internet), uses GNtimezone function in package geonames to correct to local time zone (excluding daylight saving correction)
rungads<-1 # use the Global Aerosol Database?
dailywind<-1 # use daily windspeed database?
terrain<-0 # include terrain (slope, aspect, horizon angles) (1) or not (0)?
soildata<-1 # include soil data for Australia (1) or not (0)?
snowmodel<-0 # run snow version? (slower!)
ystart <- 1990# start year for weather generator calibration dataset or AWAP database
yfinish <- 2009# end year for weather generator calibration dataset
nyears<-yfinish-ystart+1# integer, number of years for which to run the microclimate model, only for AWAP data (!!max 10 years!!)

#longlats<-read.csv("/vlsci/VR0212/shared/NicheMapR_Working/bioregion_points_final.csv")
#longlats<-subset(longlats,RASTERVALU!=0)
#lng<-which.min(abs(longlats[,1]-longlat[1]))
#coarselong<-longlats[lng,1]
#longlats<-subset(longlats,longlats$long==longlats[lng,1])
#lt<-which.min(abs(longlats[,2]-longlat[2]))
#coarselat<-longlats[lt,2]
#bioregion<-longlats[lt,3]
#barcoo<-paste('/scratch/VR0212/bio',bioregion,'/',sep="")
#load(paste(barcoo,'longlat.bin',sep=''))
#jobnum<-as.numeric(rownames(subset(data,data$V2==coarselong & data$V3==coarselat)))
quadrangle<-jobnum


############# microclimate model parameters ################################
EC <- 0.0167238 # Eccenricity of the earth's orbit (current value 0.0167238, ranges between 0.0034 to 0.058)
RUF <- 0.004 # Roughness height (m), , e.g. sand is 0.05, grass may be 2.0, current allowed range: 0.001 (snow) - 2.0 cm.
# Next for parameters are segmented velocity profiles due to bushes, rocks etc. on the surface, IF NO EXPERIMENTAL WIND PROFILE DATA SET ALL THESE TO ZERO!
Z01 <- 0. # Top (1st) segment roughness height(m)
Z02 <- 0. # 2nd segment roughness height(m)
ZH1 <- 0. # Top of (1st) segment, height above surface(m)
ZH2 <- 0. # 2nd segment, height above surface(m)
SLE <- 0.96 # Substrate longwave IR emissivity (decimal %), typically close to 1
ERR <- 2.0 # Integrator error for soil temperature calculations
DEP <- c(0.,1.5,  3.5, 5.,  10,  15,  30.,  60.,  100.,  200.) # Soil nodes (cm) - keep spacing close near the surface, last value is where it is assumed that the soil temperature is at the annual mean air temperature

#if(sitemethod==1){
#longlat <- geocode(loc)[1, 3:4] # assumes first geocode match is correct
#}

#longlat<-c(139.3109, -33.888)
prevdir<-getwd()
setwd('/hsm/VR0212/shared/CSIRO Soil and Landscape Grid')
x<-cbind(longlat[1],longlat[2])
library('raster')
library('rgdal')
files<-c('density/BDW_000_005_EV_N_P_AU_NAT_C_20140801.tif','density/BDW_005_015_EV_N_P_AU_NAT_C_20140801.tif','density/BDW_015_030_EV_N_P_AU_NAT_C_20140801.tif','density/BDW_030_060_EV_N_P_AU_NAT_C_20140801.tif','density/BDW_060_100_EV_N_P_AU_NAT_C_20140801.tif','density/BDW_100_200_EV_N_P_AU_NAT_C_20140801.tif')
bdw<-rep(NA,6)
for(ii in 1:6){
a<-raster(files[ii])
bdw[ii]<-extract(a,x)
rm(a)
gc()
}
files<-c('clay/CLY_000_005_EV_N_P_AU_NAT_C_20140801.tif','clay/CLY_005_015_EV_N_P_AU_NAT_C_20140801.tif','clay/CLY_015_030_EV_N_P_AU_NAT_C_20140801.tif','clay/CLY_030_060_EV_N_P_AU_NAT_C_20140801.tif','clay/CLY_060_100_EV_N_P_AU_NAT_C_20140801.tif','clay/CLY_100_200_EV_N_P_AU_NAT_C_20140801.tif')
cly<-rep(NA,6)
for(ii in 1:6){
a<-raster(files[ii])
cly[ii]<-extract(a,x)
rm(a)
gc()
}
files<-c('silt/SLT_000_005_EV_N_P_AU_NAT_C_20140801.tif','silt/SLT_005_015_EV_N_P_AU_NAT_C_20140801.tif','silt/SLT_015_030_EV_N_P_AU_NAT_C_20140801.tif','silt/SLT_030_060_EV_N_P_AU_NAT_C_20140801.tif','silt/SLT_060_100_EV_N_P_AU_NAT_C_20140801.tif','silt/SLT_100_200_EV_N_P_AU_NAT_C_20140801.tif')
slt<-rep(NA,6)
for(ii in 1:6){
a<-raster(files[ii])
slt[ii]<-extract(a,x)
rm(a)
gc()
}
files<-c('sand/SND_000_005_EV_N_P_AU_NAT_C_20140801.tif','sand/SND_005_015_EV_N_P_AU_NAT_C_20140801.tif','sand/SND_015_030_EV_N_P_AU_NAT_C_20140801.tif','sand/SND_030_060_EV_N_P_AU_NAT_C_20140801.tif','sand/SND_060_100_EV_N_P_AU_NAT_C_20140801.tif','sand/SND_100_200_EV_N_P_AU_NAT_C_20140801.tif')
snd<-rep(NA,6)
for(ii in 1:6){
a<-raster(files[ii])
snd[ii]<-extract(a,x)
rm(a)
gc()
}
setwd(prevdir)
soilpro<-as.data.frame(cbind(bdw,cly,slt,snd))
colnames(soilpro)<-c('blkdens','clay','silt','sand')
# pre-extracted
# soilpro<-read.csv(paste(workdir,"/soilprops.txt",sep=""),header=FALSE)
# colnames(soilpro)<-c('i','site','long','lat','desc','blkdens','clay','silt','sand')
# soilpro<-subset(soilpro,site==1)
# soilpro<-soilpro[,5:9]
   
soil_depths<-c(2.5,7.5,22.5,45,80,150)
#plot(soilpro$clay~soil_depths,ylim=c(0,100),col='red',type='l')
#points(soilpro$sand~soil_depths,ylim=c(0,100),col='orange',type='l')
#points(soilpro$silt~soil_depths,ylim=c(0,100),col='grey',type='l')
#title(main=loc)
#legend("topleft", inset=.05,
#       legend=round(soilpro[1,3:5],1),bty="n", 
#       horiz=TRUE, bg=NULL, cex=0.8)

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
 
#par(mfrow=c(2,2))
#par(mar = c(2,2,1,2) + 0.1) # margin spacing stuff 
#par(mar = c(5,5,5,5) + 0.1) # margin spacing stuff 

CampNormTbl9_1<-read.csv('/hsm/VR0212/shared/NicheMapR_Working/CampNormTbl9_1.csv')
Richter<-read.csv('/hsm/VR0212/shared/NicheMapR_Working/Richter_Table1_SI.csv')
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
   

#plot(KS~soil_depths,xlim=c(-1,200),ylim=c(0.000017,0.0058))
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
  
Thcond <- 2.5 # soil minerals thermal conductivity (W/mC)
Density <- 2560. # soil minerals density (kg/m3)
SpecHeat <- 870. # soil minerals specific heat (J/kg-K)
BulkDensity <- 1300 # soil bulk density (kg/m3)
cap<-1 # organic cap present on soil surface? (cap has lower conductivity - 0.2 W/mC - and higher specific heat 1920 J/kg-K)
SatWater <- 0.26 # volumetric water content at saturation (0.1 bar matric potential) (m3/m3)
Clay <- 20 # clay content for matric potential calculations (%)
SoilMoist <- 0 # fractional soil moisture (decimal %)
rainmult<-1 # rain multiplier for surface soil moisture (use to induce runoff), proportion
runmoist<-1 # run soil moisture model (0=no, 1=yes)?
SoilMoist_Init<-c(0.1,0.12,0.15,0.3,0.4,0.4,0.4,0.4,0.4,0.4) # initial soil water content, m3/m3
evenrain<-0 # spread daily rainfall evenly across 24hrs (1) or one event at midnight (0)
maxpool<-10000#6 # max depth for water pooling on the surface, mm (to account for runoff)
L<-c(0,0,rep(4,9),1.8,0.95,0.85,0.8,0.4,0.366,0,0)*10000
L<-c(0,0,8.18990859,7.991299442,7.796891252,7.420411664,7.059944542,6.385001059,5.768074989,4.816673431,4.0121088,1.833554792,0.946862989,0.635260544,0.804575,0.43525621,0.366052856,0,0)*10000
LAI<-0.1 # leaf area index, used to partition traspiration/evaporation from PET
PE[1:9]<-CampNormTbl9_1[3,4] #air entry potential J/kg 
KS[1:9]<-CampNormTbl9_1[3,6] #saturated conductivity, kg s/m3
BB[1:9]<-CampNormTbl9_1[3,5] #soil 'b' parameter
PE[10:13]<-CampNormTbl9_1[4,4] #air entry potential J/kg 
KS[10:13]<-CampNormTbl9_1[4,6] #saturated conductivity, kg s/m3
BB[10:13]<-CampNormTbl9_1[4,5] #soil 'b' parameter
REFL<-0.2 # soil reflectance (decimal %)
slope<-0. # slope (degrees, range 0-90)
aspect<-180. # aspect (degrees, 0 = North, range 0-360)
hori<-rep(0,24) # enter the horizon angles (degrees) so that they go from 0 degrees azimuth (north) clockwise in 15 degree intervals
PCTWET<-0 # percentage of surface area acting as a free water surface (%)
CMH2O <- 1. # precipitable cm H2O in air column, 0.1 = VERY DRY; 1.0 = MOIST AIR CONDITIONS; 2.0 = HUMID, TROPICAL CONDITIONS (note this is for the whole atmospheric profile, not just near the ground)  
TIMAXS <- c(1.0, 1.0, 0.0, 0.0)   # Time of Maximums for Air Wind RelHum Cloud (h), air & Wind max's relative to solar noon, humidity and cloud cover max's relative to sunrise          												
TIMINS <- c(0.0, 0.0, 1.0, 1.0)   # Time of Minimums for Air Wind RelHum Cloud (h), air & Wind min's relative to sunrise, humidity and cloud cover min's relative to solar noon
minshade<-0. # minimum available shade (%)
maxshade<-90. # maximum available shade (%)
runshade<-1. # run the model twice, once for each shade level (1) or just for the first shade level (0)?
manualshade<-1 # if using soildata, which includes shade, this will override the data from the database and force max shade to be the number specified above
Usrhyt <- 3# local height (cm) at which air temperature, relative humidity and wind speed calculatinos will be made 
rainwet<-1.5 # mm rain that causes soil to become 90% wet
snowtemp<-1.5 # temperature at which precipitation falls as snow (used for snow model)
snowdens<-0.4 # snow density (mg/m3)
snowmelt<-1. # proportion of calculated snowmelt that doesn't refreeze
undercatch<-1. # undercatch multipier for converting rainfall to snow
rainmelt<-0.016 # paramter in equation that melts snow with rainfall as a function of air temp
write_input<-0 # write csv files of final input to working directory? 1=yes, 0=no.
warm<-0 # uniform warming of air temperature input to simulate climate change
loop<-0 # if doing multiple years, this shifts the starting year by the integer value

# run the model
niche<-list(L=L,LAI=LAI,SoilMoist_Init=SoilMoist_Init,evenrain=evenrain,runmoist=runmoist,maxpool=maxpool,PE=PE,KS=KS,BB=BB,BD=BD,loop=loop,warm=warm,rainwet=rainwet,manualshade=manualshade,dailywind=dailywind,terrain=terrain,soildata=soildata,loc=loc,ystart=ystart,yfinish=yfinish,nyears=nyears,RUF=RUF,SLE=SLE,ERR=ERR,DEP=DEP,Thcond=Thcond,Density=Density,SpecHeat=SpecHeat,BulkDensity=BulkDensity,Clay=Clay,SatWater=SatWater,SoilMoist=SoilMoist,CMH2O=CMH2O,TIMAXS=TIMAXS,TIMINS=TIMINS,minshade=minshade,maxshade=maxshade,Usrhyt=Usrhyt,REFL=REFL,slope=slope,aspect=aspect,hori=hori,rungads=rungads,cap=cap,write_input=write_input,spatial=spatial,snowmodel=snowmodel,snowtemp=snowtemp,snowdens=snowdens,snowmelt=snowmelt,undercatch=undercatch,rainmelt=rainmelt,rainmult=rainmult,runshade=runshade)
source('NicheMapR_Setup_micro.R')
nicheout<-NicheMapR(niche)

# get output
metout1<-as.data.frame(nicheout$metout) # above ground microclimatic conditions, min shade
shadmet1<-as.data.frame(nicheout$shadmet) # above ground microclimatic conditions, max shade
soil1<-as.data.frame(nicheout$soil) # soil temperatures, minimum shade
shadsoil1<-as.data.frame(nicheout$shadsoil) # soil temperatures, maximum shade
longlat<-as.matrix(nicheout$longlat)
longlat<-nicheout$longlat
elevation<-as.numeric(nicheout$ALTT)
REFL<-as.numeric(nicheout$REFL)

library(deSolve)


# subset microclimate output files for relevant dates
ystart<-1990
yfinish<-2009
nyears<-yfinish-ystart+1
month<-1
# chose period to simulate
daystart<-paste(substr(ystart,3,4),'/01/01',sep="") # y/m/d
dayfin<-daystart # y/m/d

# key parameters to play with
mass<-319 # grams
vtmax<-38 # voluntary maximum Tb
vtmin<-32 # voluntary minimum Tb
baskthresh<-18 # min temp before animal will move to a basking spot
abs_min<-0.62 # minimum animal solar absorptivity
abs_max<-0.90 # maximum animal solar absorptivity
abs_ref<-0.76 # animal solar absorptivity, no colour change

windfact<-1 # factor to multiply predicted wind by

#metout<-read.csv(paste(microin,'/metout.csv',sep=""))[,-1]
#shadmet<-read.csv(paste(microin,'/shadmet.csv',sep=""))[,-1]
#soil<-read.csv(paste(microin,'/soil.csv',sep=""))[,-1]
#shadsoil<-read.csv(paste(microin,'/shadsoil.csv',sep=""))[,-1]
#rainfall<-read.csv(paste(microin,'/rainfall.csv',sep=""))[,-1]
tzone<-paste("Etc/GMT-",10,sep="") # doing it this way ignores daylight savings!
dates<-seq(ISOdate(ystart,1,1,tz=tzone)-3600*12, ISOdate((ystart+nyears),1,1,tz=tzone)-3600*13, by="hours")
dates<-subset(dates, format(dates, "%m/%d")!= "02/29") # remove leap years
#dates<-dates+3600*1.5
metout1<-cbind(dates,metout1)
shadmet1<-cbind(dates,shadmet1)
shadsoil1<-cbind(dates,shadsoil1)
soil1<-cbind(dates,soil1)

setwd('/hsm/VR0212/shared/NicheMapR_Working/projects/beardies')
source('OneLumpAnalytical.R') # load the analytical one lump model
source('OneLump_varenv_noskin.R') # load source for ode solver version without evaporation and Tskin

for(ystart in 1990:2009){
yfinish<-ystart
nyears<-yfinish-ystart+1
month<-1
# chose period to simulate
daystart<-paste(substr(ystart,3,4),'/01/01',sep="") # y/m/d
dayfin<-daystart # y/m/d
days<-as.numeric(as.POSIXlt(dayfin)-as.POSIXlt(daystart))
#metout<-subset(metout, format(metout$dates, "%Y")== ystart & as.numeric(format(metout$dates, "%m"))==month)
#soil<-subset(soil, format(soil$dates, "%Y")== ystart & as.numeric(format(soil$dates, "%m"))<=month)
#shadmet<-subset(shadmet, format(shadmet$dates, "%Y")== ystart & as.numeric(format(shadmet$dates, "%m"))==month)
#shadsoil<-subset(shadsoil, format(shadsoil$dates, "%Y")== ystart & as.numeric(format(shadsoil$dates, "%m"))==month)
metout<-subset(metout1, format(metout1$dates, "%Y")== ystart)
soil<-subset(soil1, format(soil1$dates, "%Y")== ystart)
shadmet<-subset(shadmet1, format(shadmet1$dates, "%Y")== ystart)
shadsoil<-subset(shadsoil1, format(shadsoil1$dates, "%Y")== ystart)

# combine relevant input fields
micro_sun_all<-cbind(metout[,1:5],metout[,8],soil[,4],metout[,13:15],metout[,6])
colnames(micro_sun_all)<-c('dates','JULDAY','TIME','TALOC','TA1.2m','VLOC','TS','ZEN','SOLR','TSKYC','RHLOC')
micro_shd_all<-cbind(shadmet[,1:5],shadmet[,8],shadsoil[,4],shadmet[,13:15],shadmet[,6])
colnames(micro_shd_all)<-c('dates','JULDAY','TIME','TALOC','TA1.2m','VLOC','TS','ZEN','SOLR','TSKYC','RHLOC')




time<-seq(0,(days+1)*60*24,60) #60 minute intervals from microclimate output
time<-time[-1]
time3<-seq(0,(days+1)*60*24,20)
time3<-time3[-1]
times2<-seq(0,(days+1)*60*24,2) #two minute intervals for prediction
time<-time*60 # minutes to seconds
times2<-times2*60 # minutes to seconds
time3<-time3*60



# constants
cp<-3073 #specific heat of flesh, J/kg-C
emis<-0.95 #emissivity of skin, -
Fo_e<-0.8 #config factor, object to IR environment, -
rho<-932 #animal density, kg/m3
# 'lometry' determines whether standard or custom shapes/surface area/volume relationships are used.
# 0=plate,1=cyl,2=ellips,3=lizard (desert iguana),4=frog (leopard frog),
# 5=custom (cylinder geometry is automatically invoked when container model operates)
lometry<-3 # organism shape (see above)
# 'custallom' below operates if lometry=5, and consists of 4 pairs of values representing 
# the parameters a and b of a relationship AREA=a*mass^b, where AREA is in cm2 and mass is in g.
# The first pair are a and b for total surface area, then a and b for ventral area, then for  
# sillhouette area normal to the sun, then sillhouette area perpendicular to the sun
customallom<-c(10.4713,.688,0.425,0.85,3.798,.683,0.694,.743) # custom allometry coefficients (see above)
shape_a<-1. 
shape_b<-3.16666666667
shape_c<-0.6666666667
FATOSK<-0.4 # configuration factor to sky
FATOSB<-0.4 # configuration factor to substrate
kflesh<-0.5 # thermal conductivity of flesh W/mK
posture<-'b' # pointing normal 'n' or parallel 'p' to the sun's rays, or average 'b'?
press<-101325 #atmospheric pressure, pa
sub_reflect<-REFL # solar reflectance of substrate
pctdif<-0.1 # proportion of solar energy that is diffuse (rather than direct beam)
q<-0 # metabolic rate (W/m3)

#elevation<-read.csv(paste(microin,'/ectoin.csv',sep=""))[1,2] # elevation
pressure<-101325 # air pressure

plotxy<-0

times_sec<-seq(0,3600*23*days,3600) # hours of day in seconds

sumstats<-matrix(data = NA, nrow = nrow(metout)/24, ncol = 8, byrow = FALSE, dimnames = NULL)
contourplot<-matrix(data = NA, nrow = nrow(metout), ncol = 5, byrow = FALSE, dimnames = NULL)

for(simday in 1:(nrow(metout)/24)){
micro_sun<-subset(micro_sun_all, micro_sun_all$JULDAY==(ystart-1990)*365+simday)
micro_shd<-subset(micro_shd_all,micro_shd_all$JULDAY==(ystart-1990)*365+simday)
#micro_shd<-subset(micro_shd_all, as.numeric(format(as.POSIXlt(micro_shd_all$dates), "%d"))==simday)

# use approxfun to create interpolations for the required environmental variables
Qsolf<- approxfun(time, micro_sun[,9], rule = 2)
Tradf<- approxfun(time, rowMeans(cbind(micro_sun[,7],micro_sun[,10]),na.rm=TRUE), rule = 2) 
Tradfshd<- approxfun(time, rowMeans(cbind(micro_shd[,7],micro_shd[,10]),na.rm=TRUE), rule = 2) 
velf<- approxfun(time, micro_sun[,6]*windfact, rule = 2)
Tairf<- approxfun(time, micro_sun[,4], rule = 2)
Tashdf<- approxfun(time, micro_shd[,4], rule = 2)
Zenf<- approxfun(time, micro_sun[,8], rule = 2)
RHf<-approxfun(time, micro_sun[,11], rule = 2)
Tsubf<- approxfun(time, micro_sun[,7], rule = 2)

times<-seq(0,3600*24*(days+1),10) # sequence of seconds for a day
hours<-times/3600



# # plot Tb in open
# if(colourchanger==1){
#   abs<-abs_max # hottest possible
# }else{
#   abs<-abs_ref
# }
# posture<-'n' # hottest possible
# indata<-list(thresh=vtmax,q=q,cp=cp,emis=emis,Fo_e=Fo_e,rho=rho,abs=abs,lometry=lometry,customallom=customallom,shape_a=shape_a,shape_b=shape_b,shape_c=shape_c,posture=posture,FATOSK=FATOSK,FATOSB=FATOSB,mass=mass,sub_reflect=sub_reflect,pctdif=pctdif)
# Tc_init<-micro_sun[1,4]
# Tbs_ode<-as.data.frame(ode(y=Tc_init,times=times,func=onelump_varenv,parms=indata))
# colnames(Tbs_ode)<-c('time','Tb','Tskin','Tcfinal','timethresh')
# Tbs_ode$time<-Tbs_ode$time/3600
 dates5<-seq(ISOdate(paste(substr(ystart,1,2),substr(daystart,1,2),sep=''),substr(daystart,4,5),substr(daystart,7,8),tz=tzone)-3600*12, ISOdate(paste(substr(ystart,1,2),substr(dayfin,1,2),sep=''),substr(dayfin,4,5),substr(dayfin,7,8),tz=tzone)-3600*12+3600*24, 10)
# Tbs_ode$datetime<-dates5


colourchanger<-0
if(colourchanger==1){
  abs<-abs_max # hottest possible
}else{
  abs<-0.85
}

interval<-24*240 # number of intervals across the day to compute
step<-(3600*24)/interval
halfstep<-step/2
shade<-0.9

dayresults<-matrix(data = NA, nrow = interval*(days+1), ncol = 9, byrow = FALSE, dimnames = NULL)
debresults<-matrix(data = NA, nrow = interval*(days+1), ncol = 23, byrow = FALSE, dimnames = NULL)

out<-0
bask<-0
for(i in 0:(interval*(days+1)-1)){
  # now select a starting hour in the shade and run model with time-varying environment
  starthr<-i
  hour<-trunc((i/interval*24))+1
  stp<-3600*24/interval
  Zen<- Zenf(i*stp)
  vel<- velf(i*stp)
  Qsol_sun<- Qsolf(i*stp)
  Tair_sun<- Tairf(i*stp)
  Trad_sun<- Tradf(i*stp)
  Qsol_shd<- Qsol_sun*(1-shade)
  Tair_shd<- Tashdf(i*stp)
  Trad_shd<- Tradfshd(i*stp)
  
  
  # get Te in sun normal for general reference
  Tc_init<-Tair_shd 
  Qsol<-Qsol_sun
  Tair<-Tair_sun
  Trad<-Trad_sun
  posture<-'n'
  if(colourchanger==1){
    abs<-abs_max
  }
  # run with time-invariant model using starting hour conditions + 1/2 hour
  times1<-step # sequence of seconds for a day
  input<-list(kflesh=kflesh,q=q,cp=cp,emis=emis,Fo_e=Fo_e,rho=rho,abs=abs,lometry=lometry,customallom=customallom,shape_a=shape_a,shape_b=shape_b,shape_c=shape_c,posture=posture,FATOSK=FATOSK,FATOSB=FATOSB,mass=mass,sub_reflect=sub_reflect,pctdif=pctdif,Qsol=Qsol,vel=vel,Tair=Tair,Trad=Trad)
  thresh<-max(vtmin,Tair_shd+0.5) # threshold body temperature at which time is required (deg C)
  Tbs<-onelump(times1,Tc_init,thresh,input)
  dayresults[i+1,8]<-Tbs$Tcf
  
  if(i==0){ # starting at night in the shade, normal to the sun ready for it to come up
    Tc_init<-Tair_shd 
    Qsol<-Qsol_shd
    Tair<-Tair_shd
    Trad<-Trad_shd
    posture<-'b'
    if(colourchanger==1){
      abs<-abs_max
    }
    # run with time-invariant model using starting hour conditions + 1/2 hour
    times1<-step # sequence of seconds for a day
    input<-list(kflesh=kflesh,q=q,cp=cp,emis=emis,Fo_e=Fo_e,rho=rho,abs=abs,lometry=lometry,customallom=customallom,shape_a=shape_a,shape_b=shape_b,shape_c=shape_c,posture=posture,FATOSK=FATOSK,FATOSB=FATOSB,mass=mass,sub_reflect=sub_reflect,pctdif=pctdif,Qsol=Qsol,vel=vel,Tair=Tair,Trad=Trad)
    thresh<-max(vtmin,Tair_shd+0.5) # threshold body temperature at which time is required (deg C)
    Tbs<-onelump(times1,Tc_init,thresh,input)
    dayresults[i+1,1]<-i*step/60
    dayresults[i+1,2]<-Tbs$Tc
    dayresults[i+1,3]<-0
    dayresults[i+1,4]<-0
    dayresults[i+1,5]<-0
    dayresults[i+1,6]<-Tbs$tau/60
    dayresults[i+1,7]<-Tbs$dTc/60
    dayresults[i+1,9]<-abs

  }else{
    if(Zen==90){ # sun is down, cooling in shade
      Tc_init<-dayresults[i,2] 
      Qsol<-Qsol_shd
      Tair<-Tair_shd
      Trad<-Trad_shd
      posture<-'b'
      if(colourchanger==1){
        abs<-abs_max
      }
      # run with time-invariant model using starting hour conditions + 1/2 hour
      times1<-step # sequence of seconds for a day
      input<-list(kflesh=kflesh,q=q,cp=cp,emis=emis,Fo_e=Fo_e,rho=rho,abs=abs,lometry=lometry,customallom=customallom,shape_a=shape_a,shape_b=shape_b,shape_c=shape_c,posture=posture,FATOSK=FATOSK,FATOSB=FATOSB,mass=mass,sub_reflect=sub_reflect,pctdif=pctdif,Qsol=Qsol,vel=vel,Tair=Tair,Trad=Trad)
      thresh<-max(vtmin,Tair_shd+0.5) # threshold body temperature at which time is required (deg C)
      Tbs<-onelump(times1,Tc_init,thresh,input)
      dayresults[i+1,1]<-i*step/60
      dayresults[i+1,2]<-Tbs$Tc
      dayresults[i+1,3]<-0
      dayresults[i+1,4]<-0
      dayresults[i+1,5]<-0
      dayresults[i+1,6]<-Tbs$tau/60
      dayresults[i+1,7]<-Tbs$dTc/60
      dayresults[i+1,9]<-abs
    }else{
      if(dayresults[i,2] < baskthresh){ # too cold to emerge and bask
        Tc_init<-dayresults[i,2] 
        Qsol<-Qsol_shd
        Tair<-Tair_shd
        Trad<-Trad_shd
        posture<-'b'
        if(colourchanger==1){
          abs<-abs_max
        }
        # run with time-invariant model using starting hour conditions + 1/2 hour
        times1<-step # sequence of seconds for a day
        input<-list(kflesh=kflesh,q=q,cp=cp,emis=emis,Fo_e=Fo_e,rho=rho,abs=abs,lometry=lometry,customallom=customallom,shape_a=shape_a,shape_b=shape_b,shape_c=shape_c,posture=posture,FATOSK=FATOSK,FATOSB=FATOSB,mass=mass,sub_reflect=sub_reflect,pctdif=pctdif,Qsol=Qsol,vel=vel,Tair=Tair,Trad=Trad)
        thresh<-max(vtmin,Tair_shd+0.5) # threshold body temperature at which time is required (deg C)
        Tbs<-onelump(times1,Tc_init,thresh,input)
        dayresults[i+1,1]<-i*step/60
        dayresults[i+1,2]<-Tbs$Tc
        dayresults[i+1,3]<-0
        dayresults[i+1,4]<-0
        dayresults[i+1,5]<-0
        dayresults[i+1,6]<-Tbs$tau/60
        dayresults[i+1,7]<-Tbs$dTc/60
        dayresults[i+1,9]<-abs
      }else{
        
        if(dayresults[i,2]<max(vtmin,Tair_shd+0.5) & out==0 & dayresults[i,2]<vtmax ){ # hasn't yet emerged, but sun is up so start to bask (really, there would be a threshold before they come out to bask)
          Tc_init<-dayresults[i,2]
          Qsol<-Qsol_sun
          Tair<-Tair_sun
          Trad<-Trad_sun
          posture<-'n'
          if(colourchanger==1){
            abs<-abs_max
          }
          # run with time-invariant model using starting hour conditions + 1/2 hour
          times1<-step # sequence of seconds for a day
          input<-list(kflesh=kflesh,q=q,cp=cp,emis=emis,Fo_e=Fo_e,rho=rho,abs=abs,lometry=lometry,customallom=customallom,shape_a=shape_a,shape_b=shape_b,shape_c=shape_c,posture=posture,FATOSK=FATOSK,FATOSB=FATOSB,mass=mass,sub_reflect=sub_reflect,pctdif=pctdif,Qsol=Qsol,vel=vel,Tair=Tair,Trad=Trad)
          thresh<-max(vtmin,Tair_shd+0.5) # threshold body temperature at which time is required (deg C)
          Tbs<-onelump(times1,Tc_init,thresh,input)
          dayresults[i+1,1]<-i*step/60
          dayresults[i+1,2]<-Tbs$Tc
          dayresults[i+1,3]<-1
          dayresults[i+1,4]<-0
          dayresults[i+1,5]<-1
          dayresults[i+1,6]<-Tbs$tau/60
          dayresults[i+1,7]<-Tbs$dTc/60
          dayresults[i+1,9]<-abs
          bask<-1
        }else{ # has emerged, check to see if needs to cool or bask again
          if(dayresults[i,2]<vtmax & dayresults[i,2] >= max(vtmin,Tair_shd+0.5)){
            if(dayresults[i,5]==3 & dayresults[i,2]>max(vtmin,Tair_shd+0.5) & dayresults[i,8]>max(vtmin,Tair_shd+0.5)){ # stay cooling
              Tc_init<-dayresults[i,2]
              Qsol<-Qsol_shd
              Tair<-Tair_shd
              Trad<-Trad_shd
              posture<-'b'
              if(colourchanger==1){
                abs<-abs_min #done
              }
              # run with time-invariant model using starting hour conditions + 1/2 hour
              times1<-step # sequence of seconds for a day
              input<-list(kflesh=kflesh,q=q,cp=cp,emis=emis,Fo_e=Fo_e,rho=rho,abs=abs,lometry=lometry,customallom=customallom,shape_a=shape_a,shape_b=shape_b,shape_c=shape_c,posture=posture,FATOSK=FATOSK,FATOSB=FATOSB,mass=mass,sub_reflect=sub_reflect,pctdif=pctdif,Qsol=Qsol,vel=vel,Tair=Tair,Trad=Trad)
              thresh<-max(vtmin,Tair_shd+0.5) # threshold body temperature at which time is required (deg C)
              Tbs<-onelump(times1,Tc_init,thresh,input)
              dayresults[i+1,1]<-i*step/60
              dayresults[i+1,2]<-Tbs$Tc
              dayresults[i+1,3]<-0
              dayresults[i+1,4]<-0
              dayresults[i+1,5]<-3
              dayresults[i+1,6]<-Tbs$tau/60
              dayresults[i+1,7]<-Tbs$dTc/60
              dayresults[i+1,9]<-abs
            }else{ # go foraging
              Tc_init<-dayresults[i,2]
              Qsol<-Qsol_sun
              Tair<-Tair_sun
              Trad<-Trad_sun
              posture<-'b'
              if(colourchanger==1){
                abs<-abs_min
              }
              # run with time-invariant model using starting hour conditions + 1/2 hour
              times1<-step # sequence of seconds for a day
              input<-list(kflesh=kflesh,q=q,cp=cp,emis=emis,Fo_e=Fo_e,rho=rho,abs=abs,lometry=lometry,customallom=customallom,shape_a=shape_a,shape_b=shape_b,shape_c=shape_c,posture=posture,FATOSK=FATOSK,FATOSB=FATOSB,mass=mass,sub_reflect=sub_reflect,pctdif=pctdif,Qsol=Qsol,vel=vel,Tair=Tair,Trad=Trad)
              thresh<-max(vtmin,Tair_shd+0.5) # threshold body temperature at which time is required (deg C)
              Tbs<-onelump(times1,Tc_init,thresh,input)
              if(colourchanger==1){
                while(Tbs$Tc<vtmin+0.25){
                 abs<-abs+0.01
                 if(abs>abs_max) break
                 input<-list(kflesh=kflesh,q=q,cp=cp,emis=emis,Fo_e=Fo_e,rho=rho,abs=abs,lometry=lometry,customallom=customallom,shape_a=shape_a,shape_b=shape_b,shape_c=shape_c,posture=posture,FATOSK=FATOSK,FATOSB=FATOSB,mass=mass,sub_reflect=sub_reflect,pctdif=pctdif,Qsol=Qsol,vel=vel,Tair=Tair,Trad=Trad)
                 thresh<-max(vtmin,Tair_shd+0.5) # threshold body temperature at which time is required (deg C)
                 Tbs<-onelump(times1,Tc_init,thresh,input)
                }
              }
              dayresults[i+1,1]<-i*step/60
              dayresults[i+1,2]<-Tbs$Tc
              dayresults[i+1,3]<-0
              if(Tbs$Tc>=vtmin+0.2){
              dayresults[i+1,4]<-1
              dayresults[i+1,5]<-2
              }else{
              dayresults[i+1,4]<-0
              dayresults[i+1,5]<-1
              }
              dayresults[i+1,6]<-Tbs$tau/60
              dayresults[i+1,7]<-Tbs$dTc/60
              dayresults[i+1,9]<-abs
              out<-1
            }
          }else{ # cooling in shade
            if(dayresults[i,2]<max(vtmin,Tair_shd+0.5) & dayresults[i,2]<vtmax & dayresults[i,8]>=max(vtmin,Tair_shd+0.5)){ # animal was forced inactive, perhaps in the shift form normal to average silhouette area, so check that the hour before Te value was >= vtmin and if so try basking/foraging again
              Tc_init<-dayresults[i,2]
              Qsol<-Qsol_sun
              Tair<-Tair_sun
              Trad<-Trad_sun
              abs1<-abs
              abs2<-abs
              if(dayresults[i,8]<vtmax){
              posture<-'n'
              if(colourchanger==1){
                abs<-abs_max
                abs1<-abs
              }
              }else{
              posture<-'b'
              if(colourchanger==1){
                abs<-abs_min#done
                abs1<-abs
              }
              }
              # run with time-invariant model using starting hour conditions + 1/2 hour
              times1<-step # sequence of seconds for a day
              input<-list(kflesh=kflesh,q=q,cp=cp,emis=emis,Fo_e=Fo_e,rho=rho,abs=abs,lometry=lometry,customallom=customallom,shape_a=shape_a,shape_b=shape_b,shape_c=shape_c,posture=posture,FATOSK=FATOSK,FATOSB=FATOSB,mass=mass,sub_reflect=sub_reflect,pctdif=pctdif,Qsol=Qsol,vel=vel,Tair=Tair,Trad=Trad)
              thresh<-max(vtmin,Tair_shd+0.5) # threshold body temperature at which time is required (deg C)
              Tbs<-onelump(times1,Tc_init,thresh,input)
              
              Qsol<-Qsol_sun
              Tair<-Tair_sun
              Trad<-Trad_sun
              posture<-'b'
              if(colourchanger==1){
                abs<-abs_max
                abs2<-abs
              }
              # run with time-invariant model using starting hour conditions + 1/2 hour
              times1<-step # sequence of seconds for a day
              input<-list(kflesh=kflesh,q=q,cp=cp,emis=emis,Fo_e=Fo_e,rho=rho,abs=abs,lometry=lometry,customallom=customallom,shape_a=shape_a,shape_b=shape_b,shape_c=shape_c,posture=posture,FATOSK=FATOSK,FATOSB=FATOSB,mass=mass,sub_reflect=sub_reflect,pctdif=pctdif,Qsol=Qsol,vel=vel,Tair=Tair,Trad=Trad)
              thresh<-max(vtmin,Tair_shd+0.5) # threshold body temperature at which time is required (deg C)
              Tbs2<-onelump(times1,Tc_init,thresh,input)
              
              if(Tbs$Tc<vtmax){ # this works so bask
                dayresults[i+1,1]<-i*step/60
                dayresults[i+1,2]<-Tbs$Tc
                dayresults[i+1,3]<-1 #normal
                dayresults[i+1,4]<-0 #inactive
                dayresults[i+1,5]<-1 #basking
                dayresults[i+1,6]<-Tbs$tau/60
                dayresults[i+1,7]<-Tbs$dTc/60
                dayresults[i+1,9]<-abs1
                # insert option to try forage if not above vtmax and greater than vtmin
                if(Tbs2$Tc<vtmax & Tbs2$Tc>max(vtmin,Tair_shd+0.5) & dayresults[i,8]<vtmax){ # this works so forage
                  dayresults[i+1,1]<-i*step/60
                  dayresults[i+1,2]<-Tbs2$Tc
                  dayresults[i+1,3]<-0 #average
                  dayresults[i+1,4]<-1 #active
                  dayresults[i+1,5]<-2 #forageing
                  dayresults[i+1,6]<-Tbs2$tau/60
                  dayresults[i+1,7]<-Tbs2$dTc/60
                  dayresults[i+1,9]<-abs2
                }
                # end insert option to try forage if not above vtmax
              }else{ # too hot, go to shade
                Tc_init<-dayresults[i,2]
                Qsol<-Qsol_shd
                Tair<-Tair_shd
                Trad<-Trad_shd
                posture<-'b'
                if(colourchanger==1){
                  abs<-abs_min #done
                }
                # run with time-invariant model using starting hour conditions + 1/2 hour
                times1<-step # sequence of seconds for a day
                input<-list(kflesh=kflesh,q=q,cp=cp,emis=emis,Fo_e=Fo_e,rho=rho,abs=abs,lometry=lometry,customallom=customallom,shape_a=shape_a,shape_b=shape_b,shape_c=shape_c,posture=posture,FATOSK=FATOSK,FATOSB=FATOSB,mass=mass,sub_reflect=sub_reflect,pctdif=pctdif,Qsol=Qsol,vel=vel,Tair=Tair,Trad=Trad)
                thresh<-max(vtmin,Tair_shd+0.5) # threshold body temperature at which time is required (deg C)
                Tbs<-onelump(times1,Tc_init,thresh,input)
                dayresults[i+1,1]<-i*step/60
                dayresults[i+1,2]<-Tbs$Tc
                dayresults[i+1,3]<-0
                dayresults[i+1,4]<-0
                dayresults[i+1,5]<-3
                dayresults[i+1,6]<-Tbs$tau/60
                dayresults[i+1,7]<-Tbs$dTc/60
                dayresults[i+1,9]<-abs
              } 
            }else{ #cooling in the shade
              Tc_init<-dayresults[i,2]
              Qsol<-Qsol_shd
              Tair<-Tair_shd
              Trad<-Trad_shd
              posture<-'b'
              if(colourchanger==1){
                if(dayresults[i,8]>vtmax){
                abs<-abs_min#done
                }else{
                abs<-abs_max
                }
              }
              # run with time-invariant model using starting hour conditions + 1/2 hour
              times1<-step # sequence of seconds for a day
              input<-list(kflesh=kflesh,q=q,cp=cp,emis=emis,Fo_e=Fo_e,rho=rho,abs=abs,lometry=lometry,customallom=customallom,shape_a=shape_a,shape_b=shape_b,shape_c=shape_c,posture=posture,FATOSK=FATOSK,FATOSB=FATOSB,mass=mass,sub_reflect=sub_reflect,pctdif=pctdif,Qsol=Qsol,vel=vel,Tair=Tair,Trad=Trad)
              thresh<-max(vtmin,Tair_shd+0.5) # threshold body temperature at which time is required (deg C)
              Tbs<-onelump(times1,Tc_init,thresh,input)
              dayresults[i+1,1]<-i*step/60
              dayresults[i+1,2]<-Tbs$Tc
              dayresults[i+1,3]<-0
              dayresults[i+1,4]<-0
              dayresults[i+1,5]<-3
              dayresults[i+1,6]<-Tbs$tau/60
              dayresults[i+1,7]<-Tbs$dTc/60
              dayresults[i+1,9]<-abs
            }
          }
        }
      }
    }
  }
  #cat(paste('hour ',hour,' done','\n'))
}


hrs<-dayresults[,1]/60
dates4<-seq(ISOdate(paste(substr(ystart,1,2),substr(daystart,1,2),sep=''),substr(daystart,4,5),substr(daystart,7,8),tz=tzone)-3600*12, ISOdate(paste(substr(ystart,1,2),substr(dayfin,1,2),sep=''),substr(dayfin,4,5),substr(dayfin,7,8),tz=tzone)-3600*12+3600*24, stp)
dates4<-seq(as.POSIXct(micro_sun[1,1]),as.POSIXct(micro_sun[1,1]+3600*24), stp)
dates4<-dates4[1:length(dates4)-1]
  dayresults<-cbind(dayresults,dates4)


        # now get metabolic rates
        #MRT (ml O2 per h) = 0.110 M 0.768 x 10(T – 20) x log10(Q10)/10, from Craig White emial 11/8/2014
        Q10<-2.44
        mrate.reptile<-(0.110*mass^0.768 * 10^((dayresults[,2]-20) * log10(Q10)/10))*0.0056*(24/interval)*3600/1000 # 0.0056 converts to Watts, then convert to kJ
        dayresults<-cbind(dayresults,mrate.reptile)
        mrate.sum<-sum(dayresults[,11])
        inactive<-subset(dayresults,dayresults[,5]==0)
        active<-subset(dayresults,dayresults[,5]==1)
        mrate.sum.inactive<-sum(inactive[,11])
        mrate.sum.active<-sum(active[,11])
       
        
        # now summarize to hourly activity times and max foraging bouts
        Hour<-trunc(dayresults[,1]/60)
        dayresults<-cbind(Hour,dayresults)
        active<-aggregate(dayresults[,5], by=list(dayresults[,1]),sum)
        active<-active$x/(interval/24)*60
        y <- rle(dayresults[,5])
        maxrun<-max((y$lengths[y$values==1]))/(interval/24)*60

        active.bouts<-y$lengths[y$values==1]
        total.bouts<-length(active.bouts)
        if(total.bouts>1){
          morning.bout<-y$lengths[2]/(interval/24)*60
          arvo.bout<-y$lengths[length(y$lengths)-1]/(interval/24)*60
        }else{
          morning.bout<-maxrun
          arvo.bout<-maxrun
        }
        if(total.bouts>2){
          midday.bouts<-active.bouts[2:(length(active.bouts)-1)]/(interval/24)*60
          midday.bout1<-midday.bouts[1]
          mean.midday.bout<-mean(midday.bouts)
        }else{
          midday.bout1<-maxrun
          mean.midday.bout<-maxrun
        }
        
        if(maxrun=='-Inf'){
          maxrun<-0
          morning.bout<-0
          midday.bout1<-0
          mean.midday.bout<-0
          arvo.bout<-0
        }
        sumact<-sum(active)
        sumstat<-t(c(micro_sun[1,2],maxrun,sumact,total.bouts,morning.bout,midday.bout1,mean.midday.bout,arvo.bout))



    sumstats[simday,]<-sumstat
  
        for(i in 0:23){
          run<-subset(dayresults,Hour==i)
          y <- rle(run[,5])
          run<-max((y$lengths[y$values==1]))/(interval/24)*60
          if(run=='-Inf'){run<-0}
          if(i==0){
            runs<-run
          }else{
            runs<-c(runs,run)
          }
        }
        contour<-cbind(micro_sun[1,2],seq(0,23,1),active,runs,micro_sun$ZEN)

          contourplot[(24*(simday-1)+1):(24*(simday-1)+24),]<-contour
  
if(plotxy==1){
plotdayresults<-as.data.frame(dayresults)
colnames(plotdayresults)<-c("Hour","Time", "Tb", "posture","active","state","tau","dTc","Te","Abs","datetime","mrate.reptile")
#plotdayresults$datetime<-as.POSIXct(plotdayresults$datetime,format=c("%Y-%m-%d %H:%M:%S"),origin="1970-10-01")
plotdates<-as.POSIXct(micro_shd$dates+3600,format=c("%Y-%m-%d %H:%M:%S"))
plot(plotdayresults$Tb~dates4,ylim=c(-5,70),type='l',col="dark green",main=as.Date(micro_shd[24,1],format=c("%Y-%m-%d")))
#points(plotdayresults$Te~plotdayresults$datetime,type='l',col="black")
points(micro_shd$TALOC~plotdates,type='l',col='blue')
points(plotdayresults$Abs*10~plotdayresults$datetime,type='l',col="orange")
#points(plotdayresults$active*5~plotdayresults$datetime,type='l',col="red")
abline(vtmax,0,col='red',lty=2)
abline(vtmin,0,col='light blue',lty=2)
#points(plotdayresults$state*5~plotdayresults$datetime,type='l',col="brown")
text(micro_shd[3,1],70,paste("bouts ",round(sumstat[,4],0),sep=""))
text(micro_shd[3,1],65,paste("maxrun ",round(sumstat[,2],0)," mins",sep=""))
text(micro_shd[3,1],60,paste("sumact ",round(sumstat[,3],0)," mins",sep=""))
text(micro_shd[3,1],55,paste("morn ",round(sumstat[,5],0)," mins",sep=""))
text(micro_shd[3,1],50,paste("mid1 ",round(sumstat[,6],0)," mins",sep=""))  
text(micro_shd[3,1],45,paste("meanmid ",round(sumstat[,7],0)," mins",sep="")) 
text(micro_shd[3,1],40,paste("arvo ",round(sumstat[,8],0)," mins",sep="")) 
}
  cat(paste('day ',(ystart-1990)*365+simday,' done \n'),sep="")
}
        
contourplot<-as.data.frame(contourplot)
sumstats<-as.data.frame(sumstats)
dates2<-seq(ISOdate(ystart,1,1,tz=tzone)-3600*12, ISOdate((ystart+nyears),1,1,tz=tzone)-3600*13, by="days")
sumstats<-cbind(dates2,sumstats)
contourplot<-cbind(dates,contourplot)

colnames(contourplot)<-c("dates","DOY","hour","forage.time.minute","forage.bout.minute","zen")
sumstats<-cbind(longlat[1],longlat[2],sumstats)
colnames(sumstats)<-c("long","lat","date","doy","maxrun","sumact","bouts","morn","mid1","meanmid","arvo")

foraging<-subset(contourplot,forage.time.minute>0)

night<-subset(contourplot,zen==90)
#with(night,plot(hour~DOY,pch=15,cex=2,col='dark blue'))
#with(foraging,points(hour~DOY,pch=15,cex=forage.time.minute/50,col='orange'))
#with(foraging,points(hour~DOY,pch=15,cex=forage.bout.minute/50,col='red'))
#sumstats

write.table(sumstats, file = paste(ystart,"_sumstats.csv",sep=""), sep = ",", col.names = F, qmethod = "double", append = T)
#write.csv(sumstats,'sumstats.csv')
#write.csv(contourplot,'MitchellPlot.csv')
}
}}
