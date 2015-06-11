############# ectotherm model parameters ################################
library(deSolve)
microin<-"microclimate" # subfolder containing the microclimate input data

# colour change options
# constant 1 0.85
# constant 2 0.72
# constant 3 0.785
# variable 1 20% change declining from a starting abs of 0.9 down to 0.82 in 1% intervals
# variable 2 increase change by 1% jumps from 1% (0.5% either side) to 10%, from base of 0.785
# variable 3 increase change by 1% jumps from 1% (0.5% either side) to 18%, from base of 0.76


# variable 1: 0.9

# get input microclimate files and read them in
 file.copy('/git/micro_australia/metout.csv',paste(microin,'/metout.csv',sep=""),overwrite=TRUE)
 file.copy('/git/micro_australia/shadmet.csv',paste(microin,'/shadmet.csv',sep=""),overwrite=TRUE)
 file.copy('/git/micro_australia/soil.csv',paste(microin,'/soil.csv',sep=""),overwrite=TRUE)
 file.copy('/git/micro_australia/shadsoil.csv',paste(microin,'/shadsoil.csv',sep=""),overwrite=TRUE)
 file.copy('/git/micro_australia/rainfall.csv',paste(microin,'/rainfall.csv',sep=""),overwrite=TRUE)
 file.copy('/git/micro_australia/ectoin.csv',paste(microin,'/ectoin.csv',sep=""),overwrite=TRUE)
 file.copy('/git/micro_australia/DEP.csv',paste(microin,'/DEP.csv',sep=""),overwrite=TRUE)
 file.copy('/git/micro_australia/MAXSHADES.csv',paste(microin,'/MAXSHADES.csv',sep=""),overwrite=TRUE)

# subset microclimate output files for relevant dates
ystart<-1990
yfinish<-ystart
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

metout<-read.csv(paste(microin,'/metout.csv',sep=""))[,-1]
shadmet<-read.csv(paste(microin,'/shadmet.csv',sep=""))[,-1]
soil<-read.csv(paste(microin,'/soil.csv',sep=""))[,-1]
shadsoil<-read.csv(paste(microin,'/shadsoil.csv',sep=""))[,-1]
rainfall<-read.csv(paste(microin,'/rainfall.csv',sep=""))[,-1]
tzone<-paste("Etc/GMT-",10,sep="") # doing it this way ignores daylight savings!
dates<-seq(ISOdate(ystart,1,1,tz=tzone)-3600*12, ISOdate((ystart+nyears),1,1,tz=tzone)-3600*13, by="hours")
dates<-subset(dates, format(dates, "%m/%d")!= "02/29") # remove leap years
#dates<-dates+3600*1.5
metout<-cbind(dates,metout)
shadmet<-cbind(dates,shadmet)
shadsoil<-cbind(dates,shadsoil)
soil<-cbind(dates,soil)

days<-as.numeric(as.POSIXlt(dayfin)-as.POSIXlt(daystart))
#metout<-subset(metout, format(metout$dates, "%Y")== ystart & as.numeric(format(metout$dates, "%m"))==month)
#soil<-subset(soil, format(soil$dates, "%Y")== ystart & as.numeric(format(soil$dates, "%m"))<=month)
#shadmet<-subset(shadmet, format(shadmet$dates, "%Y")== ystart & as.numeric(format(shadmet$dates, "%m"))==month)
#shadsoil<-subset(shadsoil, format(shadsoil$dates, "%Y")== ystart & as.numeric(format(shadsoil$dates, "%m"))==month)
metout<-subset(metout, format(metout$dates, "%Y")== ystart)
soil<-subset(soil, format(soil$dates, "%Y")== ystart)
shadmet<-subset(shadmet, format(shadmet$dates, "%Y")== ystart)
shadsoil<-subset(shadsoil, format(shadsoil$dates, "%Y")== ystart)

# combine relevant input fields
micro_sun_all<-cbind(metout[,1:5],metout[,8],soil[,4],metout[,13:15],metout[,6])
colnames(micro_sun_all)<-c('dates','JULDAY','TIME','TALOC','TA1.2m','VLOC','TS','ZEN','SOLR','TSKYC','RHLOC')
micro_shd_all<-cbind(shadmet[,1:5],shadmet[,8],shadsoil[,4],shadmet[,13:15],shadmet[,6])
colnames(micro_shd_all)<-c('dates','JULDAY','TIME','TALOC','TA1.2m','VLOC','TS','ZEN','SOLR','TSKYC','RHLOC')




time<-seq(0,(days+1)*60*24,60) #60 minute intervals from microclimate output
time3<-seq(0,(days+1)*60*24,20)
times2<-seq(0,(days+1)*60*24,2) #two minute intervals for prediction
time<-time*60 # minutes to seconds
times2<-times2*60 # minutes to seconds
time3<-time3*60


source('/git/OneLumpTrans/OneLumpAnalytical.R') # load the analytical one lump model
source('/git/OneLumpTrans/OneLump_varenv_noskin.R') # load source for ode solver version without evaporation and Tskin

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
sub_reflect<-0.2 # solar reflectance of substrate
pctdif<-0.1 # proportion of solar energy that is diffuse (rather than direct beam)
q<-0 # metabolic rate (W/m3)

elevation<-read.csv(paste(microin,'/ectoin.csv',sep=""))[1,2] # elevation
pressure<-101325 # air pressure

plotxy<-1
airoff<-0.5

times_sec<-seq(0,3600*24*1,3600) # hours of day in seconds

shade<-0.9

sumstats<-matrix(data = NA, nrow = nrow(metout)/24, ncol = 9, byrow = FALSE, dimnames = NULL)
contourplot<-matrix(data = NA, nrow = nrow(metout), ncol = 5, byrow = FALSE, dimnames = NULL)

for(simday in 1:7){#(nrow(metout)/24)){
micro_sun<-subset(micro_sun_all, micro_sun_all$JULDAY==simday)
micro_shd<-subset(micro_shd_all,micro_shd_all$JULDAY==simday)
#micro_shd<-subset(micro_shd_all, as.numeric(format(as.POSIXlt(micro_shd_all$dates), "%d"))==simday)

# use approxfun to create interpolations for the required environmental variables
Qsolf_sun<- approxfun(time, c(micro_sun[,9],(micro_sun[1,9]+micro_sun[24,9])/2), rule = 2)
Tradf_sun<- approxfun(time, rowMeans(cbind(c(micro_sun[,7],(micro_sun[1,7]+micro_sun[24,7])/24),c(micro_sun[,10],(micro_sun[1,10]+micro_sun[24,10])/24)),na.rm=TRUE), rule = 2) 
Qsolf_shd<- approxfun(time, c(micro_shd[,9],(micro_shd[1,9]+micro_shd[24,9])/2)*(1-shade), rule = 2)
Tradf_shd<- approxfun(time, rowMeans(cbind(c(micro_shd[,7],(micro_shd[1,7]+micro_shd[24,7])/24),c(micro_shd[,10],(micro_shd[1,10]+micro_shd[24,10])/24)),na.rm=TRUE), rule = 2) 
velf<- approxfun(time, c(micro_sun[,6],(micro_sun[1,6]+micro_sun[24,6])/2)*windfact, rule = 2)
Tairf_sun<- approxfun(time, c(micro_sun[,4],(micro_sun[1,4]+micro_sun[24,4])/2), rule = 2)
Tairf_shd<- approxfun(time, c(micro_shd[,4],(micro_shd[1,4]+micro_shd[24,4])/2), rule = 2)  
Zenf<- approxfun(time, c(micro_sun[,8],90), rule = 2)

#times<-seq(0,3600*24*(days+1),10) # sequence of seconds for a day
#hours<-times/3600

colourchanger<-0
if(colourchanger==1){
  abs<-abs_max # hottest possible
}else{
  abs<-0.85
}

indata<-list(thresh=vtmax,q=q,cp=cp,emis=emis,Fo_e=Fo_e,rho=rho,abs=abs,lometry=lometry,customallom=customallom,shape_a=shape_a,shape_b=shape_b,shape_c=shape_c,posture=posture,FATOSK=FATOSK,FATOSB=FATOSB,mass=mass,sub_reflect=sub_reflect,pctdif=pctdif)

emerge <- function (t, y, pars) {
  if(Zenf(t)!=90 & y>baskthresh){y<-0}
  return(y)
}
retreat <- function (t, y, pars) {
  if(Zenf(t)==90 | y<baskthresh){y<-0}
  return(y)
}
toohot <- function (t, y, pars) {
  if(y>=vtmax | y<vtmin){y<-0}
  return(y)
}
toocold <- function (t, y, pars) {
  return(y - max(vtmin,if(Tairf_shd(t)>vtmax-2){0}else{Tairf_shd(t)+0.5}))
}
eventfun <- function(t, y, pars) {
  return(y = 1)
}

morning<-function(){
Tbs_ode<-as.data.frame(ode(y=Tc_init,times=subtime,func=onelump_varenv,parms=indata,events = list(func = eventfun, root = TRUE, terminalroot = 1),
           rootfun = emerge,method='lsoda'))
colnames(Tbs_ode)<-c('time','Tb','Tcfinal','tau','dTc')  
return(Tbs_ode)
}

afternoon<-function(){
Tbs_ode<-as.data.frame(ode(y=Tc_init,times=subtime,func=onelump_varenv,parms=indata,events = list(func = eventfun, root = TRUE, terminalroot = 1),
           rootfun = retreat,method='lsoda'))
colnames(Tbs_ode)<-c('time','Tb','Tcfinal','tau','dTc')  
return(Tbs_ode)
}

warming<-function(){
Tbs_ode<-as.data.frame(ode(y=Tc_init,times=subtime,func=onelump_varenv,parms=indata,events = list(func = eventfun, root = TRUE, terminalroot = 1),
           rootfun = toohot,method='lsoda'))
colnames(Tbs_ode)<-c('time','Tb','Tcfinal','tau','dTc')  
return(Tbs_ode)
}

cooling<-function(){  
Tbs_ode<-as.data.frame(ode(y=Tc_init,times=subtime,func=onelump_varenv,parms=indata,events = list(func = eventfun, root = TRUE, terminalroot = 1),
           rootfun = toocold,method='lsoda'))
colnames(Tbs_ode)<-c('time','Tb','Tcfinal','tau','dTc')  
return(Tbs_ode)  
}

Tc_init<-Tairf_shd(0)

times<-seq(0,3600*24,10) # sequence of seconds for a day
times<-times[1:(length(times)-1)]
hours<-times/3600
times_orig<-times
out<-0
bask<-1  
daybreak<-0
posture<-'n'
rm(dayresults)
arvo<-times[(length(times)/2):length(times)]
zeniths<-as.data.frame(cbind(arvo,Zenf(arvo)))
colnames(zeniths)<-c('time','zen')
evening<-subset(zeniths,zen==90)
sunset<-evening[1,1]
times<-times[times<sunset]
subtime<-times
while(length(subtime)>0){
  if(daybreak==0){
   indata$posture<-'b'
   Tairf<-Tairf_shd
   Tradf<-Tradf_shd
   Qsolf<-Qsolf_shd
   Tbs<-morning()
   Tbs$posture<-0
   Tbs$active<-0 
   Tbs$state<-0  
   Tbs$abs<-abs  
   if(exists('dayresults')){dayresults<-rbind(dayresults,Tbs)}else{dayresults<-Tbs}  
   Tc_init<-Tbs[nrow(Tbs),2]
   subtime<-subset(times,times>Tbs[nrow(Tbs),1])
   daybreak<-1
  }
  while(bask==1 & length(subtime)>0){
   indata$posture<-'n' 
   Tairf<-Tairf_sun
   Tradf<-Tradf_sun
   Qsolf<-Qsolf_sun
   Tbs<-cooling()
   Tbs$posture<-1
   Tbs$active<-0 
   Tbs$state<-1 
   Tbs$abs<-abs
   if(exists('dayresults')){dayresults<-rbind(dayresults,Tbs)}else{dayresults<-Tbs}  
   Tc_init<-Tbs[nrow(Tbs),2]
   subtime<-subset(times,times>Tbs[nrow(Tbs),1])
   if(length(subtime)==0){break}  
   indata$posture<-'b' 
   Tbs<-warming()
   Tbs$posture<-0
   Tbs$active<-1
   Tbs$state<-2
   Tbs$abs<-abs  
   if(exists('dayresults')){dayresults<-rbind(dayresults,Tbs)}else{dayresults<-Tbs}  
   Tc_init<-Tbs[nrow(Tbs),2]
   subtime<-subset(times,times>Tbs[nrow(Tbs),1]) 
   if(length(subtime)==0){break}  
   if(Tc_init>vtmin){
     bask<-0
     Tbs$posture<-0
     Tbs$active<-1
     Tbs$state<-2
     Tbs$abs<-abs  
   }
  }
  if(length(subtime)==0){break}
  Tairf<-Tairf_shd
  Tradf<-Tradf_shd
  Qsolf<-Qsolf_shd
  Tbs<-cooling()
  Tbs$posture<-0
  Tbs$active<-0
  Tbs$state<-3 
  Tbs$abs<-abs  
  if(exists('dayresults')){dayresults<-rbind(dayresults,Tbs)}else{dayresults<-Tbs}
  Tc_init<-Tbs[nrow(Tbs),2]
  subtime<-subset(times,times>Tbs[nrow(Tbs),1])
  Tairf<-Tairf_sun
  Tradf<-Tradf_sun
  Qsolf<-Qsolf_sun
  Tbs<-warming()
  Tbs$posture<-0
  Tbs$active<-1
  Tbs$state<-2 
  Tbs$abs<-abs  
  if(exists('dayresults')){dayresults<-rbind(dayresults,Tbs)}else{dayresults<-Tbs}
  Tc_init<-Tbs[nrow(Tbs),2]
  subtime<-subset(times,times>Tbs[nrow(Tbs),1])
  Tc_init
  if(Tc_init<vtmin){
   indata$posture<-'n' 
   Tairf<-Tairf_sun
   Tradf<-Tradf_sun
   Qsolf<-Qsolf_sun
   Tbs<-afternoon()
   Tbs$posture<-1
   Tbs$active<-0 
   Tbs$state<-1   
   Tbs$abs<-abs  
   if(exists('dayresults')){dayresults<-rbind(dayresults,Tbs)}else{dayresults<-Tbs} 
   Tc_init<-Tbs[nrow(Tbs),2]
   subtime<-subset(times,times>Tbs[nrow(Tbs),1])  
   if(length(subtime)==0){break}  
   indata$posture<-'b'
   Tairf<-Tairf_shd
   Tradf<-Tradf_shd
   Qsolf<-Qsolf_shd
   Tbs<-morning()
   Tbs$posture<-0
   Tbs$active<-0 
   Tbs$state<-0  
   Tbs$abs<-abs   
   if(exists('dayresults')){dayresults<-rbind(dayresults,Tbs)}else{dayresults<-Tbs}  
   Tc_init<-Tbs[nrow(Tbs),2]
   subtime<-subset(times,times>Tbs[nrow(Tbs),1])
   if(length(subtime)==0){break}
  }
}
if(length(subtime)==0){
subtime<-evening[,1]
}
indata$posture<-'b'
Tairf<-Tairf_shd
Tradf<-Tradf_shd
Qsolf<-Qsolf_shd
Tbs<-morning()
Tbs$posture<-0
Tbs$active<-0 
Tbs$state<-0  
Tbs$abs<-abs  
if(exists('dayresults')){dayresults<-rbind(dayresults,Tbs)}else{dayresults<-Tbs}  

dayresults$state[dayresults$Tb<vtmin-0.1 & dayresults$state!=1] <- 0
dayresults$active[dayresults$Tb<vtmin-0.1] <- 0
dayresults$state[dayresults$Tb<vtmin+0.15 & dayresults$Tb>vtmin-0.15] <- 1
dayresults$active[dayresults$Tb<vtmin+0.15 & dayresults$Tb>vtmin-0.15] <- 0  
dayresults<-subset(dayresults,dayresults$time %in% times_orig)
plottime<-dayresults$time/3600
# with(dayresults,plot(Tb~plottime,type='l',col='dark green',ylim=c(-0,70),xlim=c(0,24)))
# abline(vtmax,0,col='red',lty=2)
# abline(vtmin,0,col='light blue',lty=2)
# points(Tairf_shd(times_orig)~hours,type='l',col='blue')
# with(dayresults,points(state*2~plottime,type='l',col="brown"))
# with(dayresults,points(Tcfinal~time,type='l',lty=2,col="light grey"))  
  
hrs<-dayresults[,1]/3600
dates4<-seq(ISOdate(paste(substr(ystart,1,2),substr(daystart,1,2),sep=''),substr(daystart,4,5),substr(daystart,7,8),tz=tzone)-3600*12, ISOdate(paste(substr(ystart,1,2),substr(dayfin,1,2),sep=''),substr(dayfin,4,5),substr(dayfin,7,8),tz=tzone)-3600*12+3600*24, 10)
dates4<-seq(as.POSIXct(micro_sun[1,1]),as.POSIXct(micro_sun[1,1]+3600*24), 10)
dates4<-dates4[1:length(dates4)-1]
  dayresults<-cbind(dayresults,dates4)
interval<-length(times_orig)

        # now get metabolic rates
        #MRT (ml O2 per h) = 0.110 M 0.768 x 10(T – 20) x log10(Q10)/10, from Craig White emial 11/8/2014
        Q10<-2.44
        mrate.reptile<-(0.110*mass^0.768 * 10^((dayresults[,2]-20) * log10(Q10)/10))*0.0056*(24/interval)*3600/1000 # 0.0056 converts to Watts, then convert to kJ
        dayresults<-cbind(dayresults,mrate.reptile)
        mrate.sum<-sum(dayresults[,11])
        inactive<-subset(dayresults,dayresults[,7]==0)
        active<-subset(dayresults,dayresults[,7]==1)
        mrate.sum.inactive<-sum(inactive[,11])
        mrate.sum.active<-sum(active[,11])
       
        
        # now summarize to hourly activity times and max foraging bouts
        Hour<-trunc(dayresults[,1]/3600)
        dayresults<-cbind(Hour,dayresults)
        active<-aggregate(dayresults[,7], by=list(dayresults[,1]),sum)
        active<-active$x/(interval/24)*60
        y <- rle(dayresults[,8])
        maxrun<-max((y$lengths[y$values==1]))/(interval/24)*60
        z <- rle(dayresults[,9])
        morning.bask<-z$lengths[z$values==1][1]/(interval/24)*60
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
          morning.bask<-0
        }
        sumact<-sum(active)
        sumstat<-t(c(micro_sun[1,2],maxrun,sumact,total.bouts,morning.bask,morning.bout,midday.bout1,mean.midday.bout,arvo.bout))



    sumstats[simday,]<-sumstat
  
        for(i in 0:23){
          run<-subset(dayresults,Hour==i)
          y <- rle(run[,8])
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
#colnames(plotdayresults)<-c("Hour","Time", "Tb", "posture","active","state","tau","dTc","Te","Abs","datetime","mrate.reptile")
#plotdayresults$datetime<-as.POSIXct(plotdayresults$datetime,format=c("%Y-%m-%d %H:%M:%S"),origin="1970-10-01")
plot(plotdayresults$Tb~dates4,ylim=c(-5,70),type='l',col="dark green",main=as.Date(micro_shd[24,1],format=c("%Y-%m-%d")))
#points(plotdayresults$Te~plotdayresults$datetime,type='l',col="black")
points(micro_shd$TALOC~micro_shd$dates,type='l',col='blue')
points(plotdayresults$abs*10~plotdayresults$dates4,type='l',col="orange")
#points(plotdayresults$active*5~plotdayresults$datetime,type='l',col="red")
abline(vtmax,0,col='red',lty=2)
abline(vtmin,0,col='light blue',lty=2)
points(plotdayresults$state*5~plotdayresults$dates4,type='l',col="brown")
text(micro_shd[3,1],70,paste("bouts ",round(sumstat[,4],0),sep=""))
text(micro_shd[3,1],65,paste("maxrun ",round(sumstat[,2],0)," mins",sep=""))
text(micro_shd[3,1],60,paste("sumact ",round(sumstat[,3],0)," mins",sep=""))
text(micro_shd[3,1],55,paste("mornbask ",round(sumstat[,5],0)," mins",sep=""))
text(micro_shd[3,1],50,paste("mornfor ",round(sumstat[,6],0)," mins",sep=""))
text(micro_shd[3,1],45,paste("mid1 ",round(sumstat[,7],0)," mins",sep=""))  
text(micro_shd[3,1],40,paste("meanmid ",round(sumstat[,8],0)," mins",sep="")) 
text(micro_shd[3,1],35,paste("arvo ",round(sumstat[,9],0)," mins",sep="")) 
}
  cat(paste('day ',simday,' done \n'),sep="")
}
        
contourplot<-as.data.frame(contourplot)
sumstats<-as.data.frame(sumstats)
dates2<-seq(ISOdate(ystart,1,1,tz=tzone)-3600*12, ISOdate((ystart+nyears),1,1,tz=tzone)-3600*13, by="days")
sumstats<-cbind(dates2,sumstats)
contourplot<-cbind(dates,contourplot)

colnames(contourplot)<-c("dates","DOY","hour","forage.time.minute","forage.bout.minute","zen")
colnames(sumstats)<-c("date","doy","maxrun","sumact","bouts","mornbask","mornfor","mid1","meanmid","arvo")

foraging<-subset(contourplot,forage.time.minute>0)

night<-subset(contourplot,zen==90)
with(night,plot(hour~DOY,pch=15,cex=2,col='dark blue'))
with(foraging,points(hour~DOY,pch=15,cex=forage.time.minute/50,col='orange'))
with(foraging,points(hour~DOY,pch=15,cex=forage.bout.minute/50,col='red'))
#sumstats



write.csv(sumstats,'sumstats.csv')
write.csv(contourplot,'MitchellPlot.csv')  
  
  
  
  
  
  