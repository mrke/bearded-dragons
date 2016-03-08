############# ectotherm model parameters ################################
library(deSolve)
library(NicheMapR)
source('/git/NicheMapR/R/onelump_varenv.R') # load the analytical one lump model
source('/git/NicheMapR/R/onelump_varenv_ode.R') # load source for ode solver version without evaporation and Tskin

##################### microclimate simulation #####################################################

sites<-read.csv('beardie_sites.csv',stringsAsFactors=FALSE)
#            site       lat     long
# 1 Alice Springs -23.75141 133.9174
# 2      Menindee -32.36215 142.5093
# 3         Burke -29.90803 145.7415
# 4    Long Reach -23.76256 144.0962
# 5 Murray Bridge -35.59282 139.4773
# 6     Lake Eyre -29.65207 137.7043
# 7   Coober Pedy -29.00012 134.7284
# 8       Whyalla -33.10814 137.2712
# 9       Walpeup -35.13652 142.0233
mm=9
longlat<-c(sites[mm,3],sites[mm,2])
source("../micro_australia/get.soil.R")

loc<-longlat
ystart <- 2013# start year
yfinish <- 2013# end year
nyears<-yfinish-ystart+1# integer, number of years for which to run the microclimate model

DEP <- c(0., 1.,  3, 5.,  10,  15,  30.,  60.,  100.,  200.) # Soil nodes (cm) - keep spacing close near the surface, last value is where it is assumed that the soil temperature is at the annual mean air temperature
soil.hydro<-get.soil(SLGA = 1, soilpro = 0) # extract soil parameters
PE<-soil.hydro$PE
BB<-soil.hydro$BB
BD<-soil.hydro$BD
KS<-soil.hydro$KS
PE[1:9]<-CampNormTbl9_1[3,4] #air entry potential J/kg 
KS[1:9]<-CampNormTbl9_1[3,6] #saturated conductivity, kg s/m3
BB[1:9]<-CampNormTbl9_1[3,5] #soil 'b' parameter
PE[10:13]<-CampNormTbl9_1[4,4] #air entry potential J/kg 
KS[10:13]<-CampNormTbl9_1[4,6] #saturated conductivity, kg s/m3
BB[10:13]<-CampNormTbl9_1[4,5] #soil 'b' parameter
BulkDensity <- BD[seq(1,19,2)]*1000 #soil bulk density, kg/m3

# run microclimate model to get microclimate for ectotherm model and soil temps for predicting egg development and food availability
micro<-micro_aust(loc = longlat, ystart = ystart, yfinish = yfinish, PE = PE, BB = BB, BD = 
    BD, KS = KS, BulkDensity = BulkDensity, maxshade = 90, Usrhyt = 0.03, DEP = DEP, REFL = 0.2)
#save(micro,file = 'micro.Rda')

#load('micro.Rda')
metout<-as.data.frame(micro$metout)
soil<-as.data.frame(micro$soil)
shadmet<-as.data.frame(micro$shadmet)
shadsoil<-as.data.frame(micro$shadsoil)

##################### lizard heat budget simulation #####################################################

# lizard parameters 
mass<-319 # lizard mass, grams
vtmin<-31.8 # voluntary minimum Tb for foraging
baskthresh<-15.1 # min temp before animal will move from deep shade to a basking spot
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
sub_reflect<-0.2 # solar reflectance of substrate
pctdif<-0.1 # proportion of solar energy that is diffuse (rather than direct beam)
q<-0 # metabolic rate (W/m3)
shade<-1 # fractional shade cover (to correct solar radiation by)

# subsetting appropriate microclimate conditions 
simstart<-1 # day of year to start simulation
simfinish<-365 # day of year to finish simulation
daystart<-paste(substr(ystart,3,4),'/01/01',sep="") # y/m/d
dayfin<-paste(substr(ystart,3,4),'/12/31',sep="") # y/m/d # y/m/d
tzone<-paste("Etc/GMT-",10,sep="") # doing it this way ignores daylight savings!
dates<-seq(ISOdate(ystart,1,1,tz=tzone)-3600*12, ISOdate((ystart+nyears),1,1,tz=tzone)-3600*13, by="hours")
dates<-subset(dates, format(dates, "%m/%d")!= "02/29") # remove leap years
metout<-cbind(dates,metout)
shadmet<-cbind(dates,shadmet)
shadsoil<-cbind(dates,shadsoil)
soil<-cbind(dates,soil)
days<-simfinish-simstart
metout<-subset(metout, format(metout$dates, "%Y")== ystart)
soil<-subset(soil, format(soil$dates, "%Y")== ystart)
shadmet<-subset(shadmet, format(shadmet$dates, "%Y")== ystart)
shadsoil<-subset(shadsoil, format(shadsoil$dates, "%Y")== ystart)
# combine relevant input fields from sun and shade outputs
micro_sun_all<-cbind(metout[,1:5],metout[,8],soil[,4],metout[,13:15],metout[,6])
colnames(micro_sun_all)<-c('dates','JULDAY','TIME','TALOC','TA1.2m','VLOC','TS','ZEN','SOLR','TSKYC','RHLOC')
micro_shd_all<-cbind(shadmet[,1:5],shadmet[,8],shadsoil[,4],shadmet[,13:15],shadmet[,6])
colnames(micro_shd_all)<-c('dates','JULDAY','TIME','TALOC','TA1.2m','VLOC','TS','ZEN','SOLR','TSKYC','RHLOC')
time<-seq(0,60*24,60) #60 minute intervals from microclimate output
time3<-seq(0,60*24,20)
times2<-seq(0,60*24,2) #two minute intervals for prediction
time<-time*60 # minutes to seconds
times2<-times2*60 # minutes to seconds
time3<-time3*60
times_sec<-seq(0,3600*24*1,3600) # hours of day in seconds
lastt<-0 # last time value

# empty matrix for results
sumstats<-matrix(data = NA, nrow = simfinish-simstart+1, ncol = 6, byrow = FALSE, dimnames = NULL)

# Functions and events for desolve
emerge <- function (t, y, pars) { # if sun is up and body temperature greater than threshold for basking, then trigger emerge event
  if(Zenf(t)!=90 & y>baskthresh){y<-0}
  return(y)
}
toocold <- function (t, y, pars) { # if Tb exceeds voluntary min foraging temp
  return(y - vtmin)
}
morning<-function(){
  Tbs_ode<-as.data.frame(ode(y=Tc_init,times=subtime,func=onelump_varenv_ode,parms=indata,events = list(func = eventfun, root = TRUE, terminalroot = 1),
    rootfun = emerge,method='lsoda'))
  colnames(Tbs_ode)<-c('time','Tb','Tcfinal','tau','dTc','ABS')  
  return(Tbs_ode)
}
cooling<-function(){  
  Tbs_ode<-as.data.frame(ode(y=Tc_init,times=subtime,func=onelump_varenv_ode,parms=indata,events = list(func = eventfun, root = TRUE, terminalroot = 1),
    rootfun = toocold,method='lsoda'))
  colnames(Tbs_ode)<-c('time','Tb','Tcfinal','tau','dTc','ABS')  
  return(Tbs_ode)  
}
eventfun <- function(t, y, pars) {
  return(y = 1)
}

pdf("basking_plots.pdf",paper="A4r",width=15,height=11) # doing this means you're going to make a pdf - comment this line out if you want to see them in R
par(mfrow = c(2,2)) # set up for 5 plots in 2 columns
par(oma = c(2,2,2,2) + 0.1) # margin spacing stuff
par(mar = c(3,3,2,2) + 0.1) # margin spacing stuff 
par(mgp = c(3,1,0) ) # margin spacing stuff 

# begin loop through days
for(simday in simstart:simfinish){
  
  # subset the day's microclimate conditions
  micro_sun<-subset(micro_sun_all, micro_sun_all$JULDAY==simday)
  micro_shd<-subset(micro_shd_all,micro_shd_all$JULDAY==simday)
  
  # use approxfun to create interpolations for the required environmental variables
  Qsolf_sun<- approxfun(time, c(micro_sun[,9],(micro_sun[1,9]+micro_sun[24,9])/2), rule = 2)
  Tradf_sun<- approxfun(time, rowMeans(cbind(c(micro_sun[,7],(micro_sun[1,7]+micro_sun[24,7])/24),c(micro_sun[,10],(micro_sun[1,10]+micro_sun[24,10])/24)),na.rm=TRUE), rule = 2) 
  Qsolf_shd<- approxfun(time, c(micro_shd[,9],(micro_shd[1,9]+micro_shd[24,9])/2)*(1-shade), rule = 2)
  Tradf_shd<- approxfun(time, rowMeans(cbind(c(micro_shd[,7],(micro_shd[1,7]+micro_shd[24,7])/24),c(micro_shd[,10],(micro_shd[1,10]+micro_shd[24,10])/24)),na.rm=TRUE), rule = 2) 
  velf<- approxfun(time, c(micro_sun[,6],(micro_sun[1,6]+micro_sun[24,6])/2), rule = 2)
  Tairf_sun<- approxfun(time, c(micro_sun[,4],(micro_sun[1,4]+micro_sun[24,4])/2), rule = 2)
  Tairf_shd<- approxfun(time, c(micro_shd[,4],(micro_shd[1,4]+micro_shd[24,4])/2), rule = 2)  
  Zenf<- approxfun(time, c(micro_sun[,8],90), rule = 2)
  
  Tc_init<-Tairf_shd(0) # start with Tb at shaded air temp
  
  absorbs<-c(0.77,0.92)
  for(j in 1:length(absorbs)){ 
    
    ABS<-absorbs[j] # set thsi simulation's solar absorptivity
    
    indata<-list(Tc_init=Tc_init,thresh=vtmin,q=q,Spheat=cp,EMISAN=emis,rho=rho,ABS=ABS,lometry=lometry,customallom=customallom,shape_a=shape_a,shape_b=shape_b,shape_c=shape_c,posture=posture,FATOSK=FATOSK,FATOSB=FATOSB,AMASS=mass,sub_reflect=sub_reflect,PCTDIF=pctdif,colchange=0,lastt=lastt,ABSMAX=ABS,ABSMIN=ABS)
    
    times<-seq(0,3600*17,10) # sequence of seconds for a day - just to 5pm
    hours<-times/3600
    times_orig<-times
    out<-0 # initial foraging state
    bask<-1 # initial basking state
    daybreak<-0 # initialise daybreak even counter
    if(exists('dayresults')){rm(dayresults)} # clear the results, if any already in the memory
    subtime<-times # starting times to work with
    
    # start simulation, in the shade waiting for daybreak and basking threshold
    indata$posture<-'n'
    indata$colchange<-0
    indata$lastt<-subtime[1]
    Tairf<-Tairf_shd # choose shaded environment
    Tradf<-Tradf_shd
    Qsolf<-Qsolf_shd
    Tbs<-morning() # get Tbs until sun rises and basking threshold is reached
    ABS<-Tbs$ABS
    indata$ABS<-ABS[length(ABS)]
    Tbs<-Tbs[,1:5]
    Tbs$posture<-0
    Tbs$active<-0 
    Tbs$state<-0  
    Tbs$ABS<-ABS  
    if(exists('dayresults')){dayresults<-rbind(dayresults,Tbs)}else{dayresults<-Tbs}  
    Tc_init<-Tbs[nrow(Tbs),2] # get initial temp for next behavioural phase
    subtime<-subset(times,times>Tbs[nrow(Tbs),1]) # get times post basking event, for the next behavioural phase
    
    
    if(length(subtime)>0){
      daybreak<-1 # sun has now risen
      # now basking, waiting until hits voluntary minimum foraging temperature
      indata$posture<-'n' # change posture to be normal to the sun - basking
      indata$colchange<-0
      indata$lastt<-subtime[1]
      Tairf<-Tairf_sun # choose full sun environment
      Tradf<-Tradf_sun
      Qsolf<-Qsolf_sun
      Tbs<-cooling() # simulate Tb until it reaches VTmin - i.e. until it can forage
      ABS<-Tbs$ABS
      indata$ABS<-ABS[length(ABS)]
      Tbs<-Tbs[,1:5]
      Tbs$posture<-1
      Tbs$active<-0 
      Tbs$state<-1 
      Tbs$ABS<-ABS
      if(exists('dayresults')){dayresults<-rbind(dayresults,Tbs)}else{dayresults<-Tbs}  
    }
    # simulation finished, do some data processing
    dayresults$state[dayresults$Tb<vtmin-0.1 & dayresults$state!=1] <- 0
    dayresults$active[dayresults$Tb<vtmin-0.1] <- 0
    dayresults$state[dayresults$Tb<vtmin+0.15 & dayresults$Tb>vtmin-0.15] <- 1
    dayresults$active[dayresults$Tb<vtmin+0.15 & dayresults$Tb>vtmin-0.15] <- 0  
    dayresults<-subset(dayresults,dayresults$time %in% times_orig)
    interval<-length(times_orig)
    hrs<-dayresults[,1]/3600
    dates4<-seq(ISOdate(paste(substr(ystart,1,2),substr(daystart,1,2),sep=''),substr(daystart,4,5),substr(daystart,7,8),tz=tzone)-3600*12, ISOdate(paste(substr(ystart,1,2),substr(dayfin,1,2),sep=''),substr(dayfin,4,5),substr(dayfin,7,8),tz=tzone)-3600*12+3600*24, 10)
    dates4<-seq(as.POSIXct(micro_sun[1,1]),as.POSIXct(micro_sun[1,1]+3600*17), 10)
    dayresults<-cbind(dayresults,dates4[1:nrow(dayresults)])
    interval<-length(times_orig)
    subtime<-subset(times,times>Tbs[nrow(Tbs),1]) # get times post basking event, for the next behavioural phase
    
    # now summarize lenght of morning basking bout
    Hour<-trunc(dayresults[,1]/3600)
    dayresults<-cbind(Hour,dayresults)
    z <- rle(dayresults[,9])
    if(length(subtime)>0){
      morning.bask<-z$lengths[z$values==1][1]/(interval/24)*60
    }else{
      morning.bask=NA
    }
    
    # save time to bask, max Tb, and basking time saved
    if(ABS[1]==0.77){
      sumstats[simday-simstart+1,1]=simday
      sumstats[simday-simstart+1,2]=morning.bask
      sumstats[simday-simstart+1,4]=max(dayresults$Tb)
    }else{
      sumstats[simday-simstart+1,3]=morning.bask
      sumstats[simday-simstart+1,5]=max(dayresults$Tb)
      # only save basking time if both got to vtmin
      if(sumstats[simday-simstart+1,4]>(vtmin-0.5) & sumstats[simday-simstart+1,5]>(vtmin-0.5)){ 
        sumstats[simday-simstart+1,6]=sumstats[simday-simstart+1,2]-sumstats[simday-simstart+1,3]
      }else{
        sumstats[simday-simstart+1,6]=NA
      }
    }
    
    # plot results
    plotdayresults<-as.data.frame(dayresults)
    if(ABS[1]==0.77){
      plot(micro_shd$TALOC[4:17]~micro_shd$dates[4:17],type='l',ylab='',xlab='time of day',ylim=c(5,32),col='white',xaxt = "n",main=as.Date(micro_shd[24,1],format=c("%Y-%m-%d")))
      axis.POSIXct(side = 1, x = micro_shd$dates,
        at = seq(micro_shd$dates[4], micro_shd$dates[17], "hours"), format = "%H:%M",
        las = 2)
      points(plotdayresults$Tb~plotdayresults$dates4,type='l',col="orange")
      sunrise=subset(micro_shd,TIME<60*12)
      sunrise=as.data.frame(subset(sunrise,ZEN==90))
      sunrise=sunrise[nrow(sunrise),]
      abline(v=sunrise$dates[1],col='grey',lty=2) # add line to show sunrise
      abline(vtmin,0,col='light blue',lty=2)
      text(sunrise$dates[1],vtmin-1.5,"VTmin",col="light blue",pos=4,cex=1.5)
      text(sunrise$dates[1],30,"sunrise",col="grey",srt=90,pos=2,cex=1.5)
      text(micro_shd[14,1],8.2,paste("mornbask 77% abs ",round(morning.bask,0)," mins",sep=""),col='orange',cex=1)
    }else{
      points(plotdayresults$Tb~plotdayresults$dates4,type='l',col="dark grey",lty=2)
      text(micro_shd[14,1],6.7,paste("mornbask 92% abs ",round(morning.bask,0)," mins",sep=""),col='dark grey',cex=1)
      text(micro_shd[14,1],4.9,paste("saving of ",round(sumstats[simday-simstart+1,6],0)," mins",sep=""),col='black',cex=1)
    }
    cat(paste('day ',simday,' done \n'),sep="")
  } # end loop through absorptivities
} # end loop through days
dev.off()

sumstats<-as.data.frame(sumstats)
colnames(sumstats)<-c('doy','lightbask','darkbask','lightTb','darkTb','diff')
plot(sumstats$doy,sumstats$diff,type='h',lwd=2, ylab='time saved, min', xlab='day of year',col='grey',cex.lab=1.5,cex.axis=1.3)
plot(sumstats$doy,sumstats$lightbask,type='h',lwd=2, ylab='time saved, min', xlab='day of year',col='orange')
plot(sumstats$doy,sumstats$darkbask,type='h',lwd=2, ylab='time saved, min', xlab='day of year',col='brown')
nrow(subset(sumstats,darkbask>0))
nrow(subset(sumstats,lightbask>0))
mean(sumstats$diff,na.rm=TRUE)

write.csv(sumstats,'sumstats.csv')


seasons<-c(rep('summer',59),rep('autumn',92),rep('winter',92),rep('spring',91),rep('summer',31))

summer=subset(sumstats, sumstats[,1]<60 | sumstats[,1]>334)
autumn=subset(sumstats, sumstats[,1]>60 & sumstats[,1]<=151)
winter=subset(sumstats, sumstats[,1]>151 & sumstats[,1]<=243)
spring=subset(sumstats, sumstats[,1]>243 & sumstats[,1]<=334)
boxplot(summer[,6])
