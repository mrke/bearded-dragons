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
yfinish<-1990
nyears<-yfinish-ystart+1

# chose period to simulate
daystart<-'90/01/01' # y/m/d
dayfin<-'90/01/31' # y/m/d

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

# combine relevant input fields
micro_sun_all<-cbind(metout[,1:5],metout[,8],soil[,4],metout[,13:15],metout[,6])
colnames(micro_sun_all)<-c('dates','JULDAY','TIME','TALOC','TA1.2m','VLOC','TS','ZEN','SOLR','TSKYC','RHLOC')
micro_shd_all<-cbind(shadmet[,1:5],shadmet[,8],shadsoil[,4],shadmet[,13:15],shadmet[,6])
colnames(micro_shd_all)<-c('dates','JULDAY','TIME','TALOC','TA1.2m','VLOC','TS','ZEN','SOLR','TSKYC','RHLOC')


days<-as.numeric(as.POSIXlt(dayfin)-as.POSIXlt(daystart))
metout_sub<-subset(metout, format(metout$dates, "%y/%m/%d")>= daystart & format(metout$dates, "%y/%m/%d")<=dayfin)
soil_sub<-subset(soil, format(soil$dates, "%y/%m/%d")>= daystart & format(soil$dates, "%y/%m/%d")<=dayfin)
shadmet_sub<-subset(shadmet, format(shadmet$dates, "%y/%m/%d")>= daystart & format(shadmet$dates, "%y/%m/%d")<=dayfin)
shadsoil_sub<-subset(shadsoil, format(shadsoil$dates, "%y/%m/%d")>= daystart & format(shadsoil$dates, "%y/%m/%d")<=dayfin)

micro_sun<-subset(micro_sun_all, format(as.POSIXlt(micro_sun_all$dates), "%y/%m/%d")>=daystart & format(as.POSIXlt(micro_sun_all$dates), "%y/%m/%d")<=dayfin)
micro_shd<-subset(micro_shd_all, format(as.POSIXlt(micro_shd_all$dates), "%y/%m/%d")>=daystart & format(as.POSIXlt(micro_shd_all$dates), "%y/%m/%d")<=dayfin)


time<-seq(0,(days+1)*60*24,60) #60 minute intervals from microclimate output
time<-time[-1]
time3<-seq(0,(days+1)*60*24,20)
time3<-time3[-1]
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



times_sec<-seq(0,3600*23*days,3600) # hours of day in seconds

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

summary stats:
  - 

dayresults<-as.data.frame(dayresults)
colnames(dayresults)<-c("Time", "Tb", "posture","active","state","tau","dTc","Te","Abs")
hrs<-dayresults$Time/60
dates4<-seq(ISOdate(paste(substr(ystart,1,2),substr(daystart,1,2),sep=''),substr(daystart,4,5),substr(daystart,7,8),tz=tzone)-3600*12, ISOdate(paste(substr(ystart,1,2),substr(dayfin,1,2),sep=''),substr(dayfin,4,5),substr(dayfin,7,8),tz=tzone)-3600*12+3600*24, stp)
dayresults$datetime<-dates4[1:length(dates4)-1]

plotdates<-as.POSIXct(shadmet_sub$dates+3600,format=c("%Y-%m-%d %H:%M:%S"))
plot(dayresults$Tb~dayresults$datetime,ylim=c(-5,70),type='l',col="dark green")
#points(dayresults$Te~dayresults$datetime,type='l',col="black")
points(shadmet_sub$TALOC~plotdates,type='l',col='blue')
points(dayresults$Abs*10~dayresults$datetime,type='l',col="orange")
#points(dayresults$active*5~dayresults$datetime,type='l',col="red")
abline(vtmax,0,col='red',lty=2)
abline(vtmin,0,col='light blue',lty=2)
points(dayresults$state*5~dayresults$datetime,type='l',col="brown")
sum(dayresults$active/interval*24)

        # now get metabolic rates
        #MRT (ml O2 per h) = 0.110 M 0.768 x 10(T – 20) x log10(Q10)/10, from Craig White emial 11/8/2014
        Q10<-2.44
        mrate.reptile<-(0.110*mass^0.768 * 10^((dayresults$Tb-20) * log10(Q10)/10))*0.0056*(24/interval)*3600/1000 # 0.0056 converts to Watts, then convert to kJ
        dayresults<-cbind(dayresults,mrate.reptile)
        mrate.sum<-sum(dayresults$mrate.reptile)
        inactive<-subset(dayresults,active==0)
        active<-subset(dayresults,active==1)
        mrate.sum.inactive<-sum(inactive$mrate.reptile)
        mrate.sum.active<-sum(active$mrate.reptile)
        
        
        # now summarize to hourly activity times and max foraging bouts
        Hour<-trunc(dayresults$Time/60)
        dayresults<-cbind(Hour,dayresults)
        active<-aggregate(dayresults$active, by=list(dayresults$Hour),sum)
        active<-active$x/(interval/24)*60
        y <- rle(dayresults$active)
        maxrun<-max((y$lengths[y$values==1]))/(interval/24)*60
        if(maxrun=='-Inf'){maxrun<-0}
        sumact<-sum(active)
        sumstat<-t(c(month,maxrun,sumact,mrate.sum,mrate.sum.inactive,mrate.sum.active))

          sumstats<-sumstat
          colnames(sumstats)<-c("month","maxrun","sumact","mrate.sum","mrate.sum.inactive","mrate.sum.active")
        for(i in 0:23){
          run<-subset(dayresults,Hour==i)
          y <- rle(run$active)
          run<-max((y$lengths[y$values==1]))/(interval/24)*60
          if(run=='-Inf'){run<-0}
          if(i==0){
            runs<-run
          }else{
            runs<-c(runs,run)
          }
        }
        contour<-cbind(month,seq(0,23,1),active,runs,ZEN)

          contourplot<-contour

        
      colnames(contourplot)<-c("month","hour","forage.time.minute","forage.bout.minute","zen")
      contourplot<-as.data.frame(contourplot)
      
      foraging<-subset(contourplot,forage.time.minute>0)
      night<-subset(contourplot,zen==90)
#sumstats
text(2,60,paste("maxrun ",round(sumstats[,2],0)," mins",sep=""))
text(2,55,paste("sumact ",round(sumstats[,3],0)," mins",sep=""))
text(2,50,paste("mrate ",round(sumstats[,3],0)," kJ",sep=""))
text(2,45,paste("mrate inactive ",round(sumstats[,4]/sumstats[,3]*100,0),"%",sep=""))








