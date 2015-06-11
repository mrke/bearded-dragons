############# ectotherm model parameters ################################
library(deSolve)
microin<-"microclimate" # subfolder containing the microclimate input data

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
daystart<-'90/01/10' # y/m/d
dayfin<-'90/01/14' # y/m/d

# key parameters to play with
mass<-319 # grams
vtmax<-40 # voluntary maximum Tb
vtmin<-32 # voluntary minimum Tb
baskthresh<-10 # min temp before animal will move to a basking spot
abs_min<-0.73 # minimum animal solar absorptivity
abs_max<-0.92 # maximum animal solar absorptivity
abs_ref<-0.85 # animal solar absorptivity, no colour change
colourchanger<-1 # if this is 1, then animal will choose between abs_min and abs_max otherwise just stays at value for abs


windfact<-1 # factor to multiply predicted wind by

metout<-read.csv(paste(microin,'/metout.csv',sep=""))[,-1]
shadmet<-read.csv(paste(microin,'/shadmet.csv',sep=""))[,-1]
soil<-read.csv(paste(microin,'/soil.csv',sep=""))[,-1]
shadsoil<-read.csv(paste(microin,'/shadsoil.csv',sep=""))[,-1]
rainfall<-read.csv(paste(microin,'/rainfall.csv',sep=""))[,-1]
tzone<-paste("Etc/GMT-",11,sep="") # doing it this way ignores daylight savings!
dates<-seq(ISOdate(ystart,1,1,tz=tzone)-3600*12, ISOdate((ystart+nyears),1,1,tz=tzone)-3600*13, by="hours")
dates<-subset(dates, format(dates, "%m/%d")!= "02/29") # remove leap years
dates<-dates+3600*1.5
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


colourchanger<-1
if(colourchanger==1){
  abs<-abs_max # hottest possible
}else{
  abs<-abs_ref
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
#   if(i==0){
#   debresults[i+1,]<-DEBmodel(dayresults[i+1,2],dayresults[i+1,5],E_pres,V_pres,E_H_pres,q_pres,hs_pres,surviv_pres,ms_pres,0,0)#,E_pres,V_pres,E_H_pres,q_pres,hs_pres,surviv_pres,ms_pres,0,0
#   }else{
#     E_init<-debresults[i,2]
#     v_init<-debresults[i,3]
#     E_H_init<-debresults[i,4]
#     q_init<-debresults[i,5]
#     hs_init<-debresults[i,6]
#     surviv_init<-debresults[i,7]
#     ms_init<-debresults[i,8]
#     cumrepro_prev<-debresults[i,9]
#     cumbatch_prev<-debresults[i,10]
#     debresults[i+1,]<-DEBmodel(dayresults[i+1,2],dayresults[i+1,4],E_init,v_init,E_H_init,q_init,hs_init,surviv_init,ms_init,cumrepro_prev,cumbatch_prev)
#   }
}



dayresults<-as.data.frame(dayresults)
colnames(dayresults)<-c("Time", "Tb", "posture","active","state","tau","dTc","Te","Abs")
hrs<-dayresults$Time/60
dates4<-seq(ISOdate(paste(substr(ystart,1,2),substr(daystart,1,2),sep=''),substr(daystart,4,5),substr(daystart,7,8),tz=tzone)-3600*12, ISOdate(paste(substr(ystart,1,2),substr(dayfin,1,2),sep=''),substr(dayfin,4,5),substr(dayfin,7,8),tz=tzone)-3600*12+3600*24, stp)
dayresults$datetime<-dates4[1:length(dates4)-1]


#with(Tbs_ode,plot(Tb~datetime,type='l',col='1',ylim=c(-5,70)))
#with(Tbs,points(Tcfinal~time,type='l',col='red',ylim=c(-10,70)))

plot(dayresults$Tb~dayresults$datetime,ylim=c(-5,70),type='l',col="dark green")
#points(dayresults$Te~dayresults$datetime,type='l',col="black")
points(shadmet$TALOC~metout$dates,type='l',col='blue')
points(dayresults$Abs*10~dayresults$datetime,type='l',col="orange")
points(dayresults$active*5~dayresults$datetime,type='l',col="red")
abline(vtmax,0,col='red',lty=2)
abline(vtmin,0,col='light blue',lty=2)
#points(dayresults$state*5~dayresults$datetime,type='l',col="brown")

#DEBresults<-as.data.frame(debresults)
#colnames(DEBresults)<-c("hour","E_pres","V_pres","E_H_pres","q_pres","hs_pres","surviv_pres","ms_pres","cumrepro","cumbatch","O2FLUX","CO2FLUX","MLO2","GH2OMET","DEBQMET","DRYFOOD","FAECES","NWASTE","wetgonad","wetstorage","wetfood","wetmass","gutfreemass")
#DEBresults$datetime<-dayresults$datetime
#plot(DEBresults$ms_pres,type='l')




















skinwet<-0.5 # estimated from data in Bently 1959 at 23 degrees and 34.5 degrees #0.2#0.35 # %, of surface area acting like a free water surface (e.g. most frogs are 100% wet, many lizards less than 5% wet)
extref<-20. # %, oxygen extraction efficiency (need to check, but based on 35 deg C for a number of reptiles, from Perry, S.F., 1992. Gas exchange strategies in reptiles and the origin of the avian lung. In: Wood, S.C., Weber, R.E., Hargens, A.R., Millard, R.W. (Eds.), Physiological Adaptations in Vertebrates: Respiration, Circulation, andMetabo -  lism. Marcel Dekker, Inc., New York, pp. 149-167.)
PFEWAT<-73. # %, fecal water (from Shine's thesis, mixed diet 75% clover, 25% mealworms)
PTUREA<-0. # %, water in excreted nitrogenous waste
FoodWater<-82#82 # 82%, water content of food (from Shine's thesis, clover)
minwater<-9.5 # %, minimum tolerated dehydration (% of wet mass) - prohibits foraging if greater than this
raindrink<-5. # daily rainfall (mm) required for animal to rehydrate from drinking (zero means standing water always available)
gutfill<-75. # % gut fill at which satiation occurs - if greater than 100%, animal always tries to forage




fract<-1
f<-1.
MsM<-186.03*6. # produces a stomach volume of 5.3 cm3/100 g, as measured for Disosaurus dorsalis
z<-7.174*fract
delta<- 0.217
kappa_X<-0.85#0.85
v_dotref<-0.05591/interval
kappa<-0.8501 
p_Mref<-45.14/interval
E_G<-7189
k_R<-0.95
k_J<-0.00628/interval
E_Hb<-6.533e+04*fract^3
E_Hj<-E_Hb*fract^3
E_Hp<-1.375e+05*fract^3
h_aref<-3.61e-11/(interval^2) 
s_G<-0.01

E_Egg<-1.04e+06*fract^3# J, initial energy of one egg # this includes the residual yolk, which is eaten upon hatching
E_m<-(p_Mref*z/kappa)/v_dotref
p_Xm<-12420/interval # J/h.cm2, maximum intake rate when feeding
K<-10 # half-saturation constant
X<-3265 # food density J/cm2

# these next five parameters control the thermal response, effectively generating a thermal response curve
T_REF<-20 # degrees C, reference temperature - correction factor is 1 for this temperature
TA<-7130
TAL<-5.305e+04
TAH<-9.076e+04
TL<-288.
TH<-315.

arrhenius<-matrix(data = 0, nrow = 8, ncol = 5)
arrhenius[,1]<-TA # critical thermal minimum
arrhenius[,2]<-TAL # critical thermal maximum
arrhenius[,3]<-TAH # voluntary thermal minimum
arrhenius[,4]<-TL # voluntary thermal maximum
arrhenius[,5]<-TH # basking threshold 


water_stages<-matrix(data = 0, nrow = 8, ncol = 8)

water_stages[,1]<-skinwet
water_stages[,2]<-extref
water_stages[,3]<-PFEWAT
water_stages[,4]<-PTUREA
water_stages[,5]<-FoodWater
water_stages[,6]<-minwater
water_stages[,7]<-raindrink
water_stages[,8]<-gutfill

# composition related parameters
andens_deb<-1. # g/cm3, density of structure 
d_V<-0.3 # density of structure (reflects fraction of mass that is dry)
d_E<-0.3 # density of reserve (reflects fraction of mass that is dry)
eggdryfrac<-0.3 # decimal percent, dry mass of eggs
mu_X<-525000 # J/cmol, chemical potential of food
mu_E<-585000 # J/cmol, chemical potential of reserve
mu_V<-500000 # J/cmol, chemical potential of structure 
mu_P<-480000 # J/cmol, chemical potential of product (faeces)
kappa_X_P<-0.1 # fraction of food energy into faeces

          # elemental maxtrix of organics  
nX<-c(1,1.8,0.5,.15) # composition of food (atoms per carbon atoms for CHON)
nE<-c(1,1.8,0.5,.15) # composition of reserve (atoms per carbon atoms for CHON)
nV<-c(1,1.8,0.5,.15) # composition of structure (atoms per carbon atoms for CHON)
nP<-c(1,1.8,0.5,.15) # composition of product/faeces (atoms per carbon atoms for CHON)
N_waste<-c(5,4,3,4) # chemical formula for nitrogenous waste product, CHON, e.g. Urea c(0,3,0,1), Uric acid c(5,4,3,4)


# breeding life history
clutchsize<-30. # clutch size
eggmass<-3.787 # initial dry mass of an egg (g)
viviparous<-0 # 1=yes, 0=no
batch<-1 # invoke Pequerie et al.'s batch laying model?

# the following four parameters apply if batch = 1, i.e. animal mobilizes
breedrainthresh<-0 # rain dependent breeder? 0 means no, otherwise enter rainfall threshold in mm
# photoperiod response triggering ovulation, none (0), summer solstice (1), autumnal equinox (2),  
# winter solstice (3), vernal equinox (4), specified daylength thresholds (5)
photostart<- 3 # photoperiod initiating breeding
photofinish<- 1 # photoperiod terminating breeding
daylengthstart<- 12.5 # threshold daylength for initiating breeding
daylengthfinish<- 13.8 # threshold daylength for terminating breeding
photodirs <- 1 # is the start daylength trigger during a decrease (0) or increase (1) in day length?
photodirf <- 1 # is the finish daylength trigger during a decrease (0) or increase (1) in day length?
startday<-1 # make it 90 for T. rugosa loop day of year at which DEB model starts
breedtempthresh<-200 # body temperature threshold below which breeding will occur
breedtempcum<-24*7 # cumulative time below temperature threshold for breeding that will trigger breeding

reset<-0 # reset options, 0=quit simulation upon death, 1=restart at emergence, 2=restart at first egg laid, 3=restart at end of breeding season, 4=reset at death

# frog breeding mode 0 is off, 
# 1 is exotrophic aquatic (eggs start when water present in container and within breeding season)
# 2 is exotrophic terrestrial/aquatic (eggs start at specified soil node within breeding season, 
# diapause at birth threshold, start larval phase if water present in container)
# 3 endotrophic terrestrial (eggs start at specified soil node within breeding season and continue
# to metamorphosis on land)
# 4 turtle mode (eggs start at specified soil node within breeding season, hatch and animals enter
# water and stay there for the rest of their life, but leave the water if no water is present)
frogbreed<-0 # frog breeding mode
frogstage<-0 # 0 is whole life cycle, 1 is just to metamorphosis ({ reset and start again)

# metabolic depression
aestivate<-0
depress<-0.2



#*********************************** DEB model initial conditions **************************************


v_init<-3e-9
E_init<-E_Egg/v_init
E_H_init<-0
stage<-0
v_init<-(3.82^3)*fract^3 #hatchling
E_init<-E_m
E_H_init<-E_Hb+5
stage<-1
v_init<-(7.063^3)*fract^3*0.85
E_init<-E_m
E_H_init<-E_Hp+1
stage<-3
ma<-1e-4  # hourly active mortality rate (probability of mortality per hour)
mi<-0  # hourly inactive mortality rate (probability of mortality per hour)
mh<-0.5   # survivorship of hatchling in first year
# DEB model initial conditions
V_init_baby<-3e-9
E_init_baby<-E_Egg/V_init_baby
E_baby_init<-E_init_baby
V_baby_init<-V_init_baby
ms_init<-0.
cumrepro_init<-0.
q_init<-0.
hs_init<-0.
cumbatch_init<-0.
pregnant<-0
E_m<-(p_Mref*z/kappa)/v_dotref

# conversions from percent to proportion
PTUREA1<-PTUREA/100
PFEWAT1<-PFEWAT/100
FoodWater1<-FoodWater/100
water_stages[,3]<-water_stages[,3]/100
water_stages[,4]<-water_stages[,4]/100
water_stages[,5]<-water_stages[,5]/100
eggmass<-0 # initial dry mass of an egg (g) - no longer used so delete



#******************** DEB mass balance calculations ************************

nO<-cbind(nX,nV,nE,nP) # matrix of composition of organics, i.e. food, structure, reserve and faeces
CHON<-c(12,1,16,14)
wO<-CHON%*%nO
w_V=wO[3]
M_V<-d_V/w_V
yEX<-kappa_X*mu_X/mu_E # yield of reserve on food
yXE<-1/yEX # yield of food on reserve
yVE<-mu_E*M_V/E_G  # yield of structure on reserve
yPX<-kappa_X_P*mu_X/mu_P # yield of faeces on food
yXP<-1/yPX # yield of food on faeces
yPE<-yPX/yEX # yield of faeces on reserve  0.143382353
nM<-matrix(c(1,0,2,0,0,2,1,0,0,0,2,0,N_waste),nrow=4)
N_waste_inv<-c(-1*N_waste[1]/N_waste[4],(-1*N_waste[2])/(2*N_waste[4]),(4*N_waste[1]+N_waste[2]-2*N_waste[3])/(4*N_waste[4]),1/N_waste[4])
nM_inv<-matrix(c(1,0,-1,0,0,1/2,-1/4,0,0,0,1/2,0,N_waste_inv),nrow=4)
JM_JO<--1*nM_inv%*%nO
etaO<-matrix(c(yXE/mu_E*-1,0,1/mu_E,yPE/mu_E,0,0,-1/mu_E,0,0,yVE/mu_E,-1/mu_E,0),nrow=4)
w_N<-CHON%*%N_waste


############################# end input #########################################

# fixes Mike made to get it to work in R, including case-sensitive issues and other symbol changes

w_X=wO[1]
w_E=wO[3]
w_V=wO[2]
w_P=wO[4]

T_A<-TA
T_ref<-T_REF
E_egg<-E_Egg
ANDENS_deb<-andens_deb
k_Jref<-k_J
zfact<-z
vdotref<-v_dotref
p_Xmref<-p_Xm
waiting<-0
hour<-1
daycount<-1
lambda=6./12.
breeding<-1
surviv_init<-1
delta_deb<-delta
halfsat<-K


funct<-f
fecundity<-0
clutches<-0

cumrepro_prev<-0
cumbatch_prev<-0
cumbatch_init<-0
cumrepro_init<-0


stage<-3
dead<-0

E_pres<-E_init
V_pres<-v_init
E_H_pres<-E_H_init
q_pres<-q_init
hs_pres<-hs_init
surviv_pres<-surviv_init
ms_pres<-ms_init
ms<-ms_pres
cumrepro<-0
E_H<-E_H_pres
hs<-hs_pres
q<-q_pres
p_B_past<-0
# end fixes Mike made to get it to work in R

# *********************************************** end DEB MODEL  ***********************************************
# **************************************************************************************************************




DEBmodel<-function(Tbody,Active,E_init,v_init,E_H_init,q_init,hs_init,surviv_init,ms_init,cumrepro_prev,cumbatch_prev){

  acthr<-Active# Reports true if turtle is in food. For cum. feeding bouts: acthr<- NLReport("feedcount") (# hours of activity. This would change with input from IBM) 
  X_food<-X

  
  E_pres<-E_init
  V_pres<-v_init
  E_H_pres<-E_H_init
  q_pres<-q_init
  hs_pres<-hs_init
  surviv_pres<-surviv_init
  ms_pres<-ms_init
  
#    Arrhenius temperature correction factor
Tcorr = exp(T_A*(1/(273+T_ref)-1/(273+Tbody)))/(1+exp(TAL*(1/(273+Tbody)-1/TL))+exp(TAH*(1/TH-1/(273+Tbody))))

clutchenergy = E_egg*clutchsize

#f=funct

M_V = d_V/w_V
p_Mv = p_Mref*Tcorr
k_Mdot = p_Mv/E_G
k_J = k_Jref*Tcorr

p_Am = p_Mv*zfact/kappa
vdot = vdotref*Tcorr
E_M = p_Am/vdot
p_Xm = p_Xmref*Tcorr
g = E_G/(kappa*E_M)
E_scaled=E_pres/E_m
V_max=(kappa*p_Am/p_Mv)^(3.)
h_a = h_aref*Tcorr
L_T = 0.
L_pres = V_pres^(1./3.)
L_max = V_max^(1./3.)
scaled_l = L_pres/L_max
kappa_G = (d_V*mu_V)/(w_V*E_G)
yEX=kappa_X*mu_X/mu_E
yXE=1/yEX
yPX=kappa_X_P*mu_X/mu_P
mu_AX=mu_E/yXE
eta_PA=yPX/mu_AX

#    now checking to see if starting with embryo, and if so setting the appropriate reserve density
if(hour==1){
if(daycount==1){
if(E_H_pres<=E_Hb){
#       E_pres=E_egg/V_pres
E_pres=E_init
}
}
#    checking to see if animal died recently and needs to start again as an embryo
if((daycount>1) & (dead==1)){
if(E_H_pres<=E_Hb){
E_pres=E_egg/debfirst(3)
}
}
}

if(E_H_pres<=E_Hb){
#     use embryo equation for length, from Kooijman 2009 eq. 2
if(waiting==1){
dLdt = 0
V_temp=(V_pres^(1./3.)+dLdt)^3
dVdt = 0
rdot=0
}else{
  dLdt=(vdot*E_scaled-k_Mdot*g*V_pres^(1./3.))/(3*(E_scaled+g))
V_temp=(V_pres^(1./3.)+dLdt)^3
dVdt = V_temp-V_pres
rdot=vdot*(e_scaled/L_pres-(1+L_T/L_pres)/L_max)/(E_scaled+g)
}
}else{
#    equation 2.21 from DEB3
rdot=vdot*(E_scaled/L_pres-(1+L_T/L_pres)/L_max)/(E_scaled+g)
dVdt = V_pres*rdot
if(dVdt<0){
dVdt=0
}
}


# Calculating body length at first hour of simulation 
if(hour==1){
if(E_H_pres<=E_Hb){
#      use embryo equation for scaled reserve, U_E, from Kooijman 2009 eq. 1
Sc = L_pres^2*(g*e_scaled)/(g+E_scaled)*(1+((k_Mdot*L_pres)/vdot))
dUEdt = -1*Sc 
E_temp=((E_pres*V_pres/p_Am)+dUEdt)*p_Am/(v_pres+dvdt)
dEdt=E_temp-E_pres
}else{
  if(ms_init>0.0000001*MsM*V_pres){
#        Equation 2.10 DEB3
dEdt = (p_Am*f-E_pres*vdot)/L_pres
}else{
  dEdt = (p_Am*0-E_pres*vdot)/L_pres
}
}
}else{
  if(E_H_pres<=E_Hb){
#      use embryo equation for scaled reserve, U_E, from Kooijman 2009 eq. 1
Sc = L_pres^2*(g*e_scaled)/(g+E_scaled)*(1+((k_Mdot*L_pres)/vdot))
dUEdt = -1*Sc 
E_temp=((E_pres*V_pres/p_Am)+dUEdt)*p_Am/(v_pres+dvdt)
dEdt=E_temp-E_pres
}else{
  if(ms>0.0000001*MsM*V_pres){
dEdt = (p_Am*f-E_pres*vdot)/L_pres
}else{
  dEdt = (p_Am*0-E_pres*vdot)/L_pres
}
}
}

p_M = p_Mv*V_pres
p_J = k_J*E_H_pres

#    powers
if(hour==1){
if(ms_init>0.0000001*MsM*V_pres){
p_A = V_pres^(2./3.)*p_Am*f
}else{
  p_A = 0
}
}else{
  if(ms>0.0000001*MsM*V_pres){
p_A = V_pres^(2./3.)*p_Am*f
}else{
  p_A = 0
}
}
#    J food eaten per hour
p_X = p_A/kappa_X




#    equation 2.20 DEB3
p_C = (E_m*(vdot/L_pres+k_Mdot*(1+L_T/L_pres))*(E_scaled*g)/(E_scaled+g))*V_pres


p_R = (1.-kappa)*p_C-p_J


if((E_H_pres<=E_Hp) | (pregnant==1)){
p_B = 0.
}else{
  if(batch==1){
batchprep=(k_R/lambda)*((1-kappa)*(E_m*(vdot*V_pres^(2./3.)+k_Mdot*V_pres)/(1+(1/g)))-p_J)
if(breeding==0){
p_B =0.
}else{ 
  if(hour==1){
#        if the repro buffer is lower than what p_B would be(see below), p_B is p_R
if(cumrepro_init<batchprep){        
p_B = p_R
}else{
#          otherwise it is a faster rate, as specified in Pecquerie et. al JSR 2009 Anchovy paper, 
#         with lambda (the fraction of the year the animals breed if food/temperature not limiting) = 0.583 or 7 months of the year
p_B = batchprep
}
}else{
#       if the repro bufffer is lower than what p_B would be(see below), p_B is p_R
if(cumrepro<batchprep){        
p_B = p_R
}else{
#         otherwise it is a faster rate, as specified in Pecquerie et. al JSR 2009 Anchovy paper, 
#        with lambda (the fraction of the year the animals breed if food/temperature not limiting) = 0.583 or 7 months of the year
p_B = batchprep
}
}
}
}else{
  p_B=p_R
#     end check for whether batch mode is operating
}
#    end check for immature or mature
}

#maturity
if(E_H_pres<E_Hp){
if(E_H_pres<=E_Hb){
#       use embryo equation for scaled maturity, U_H, from Kooijman 2009 eq. 3
if(waiting==1){
U_H_pres=E_H_pres/p_Am    
dUHdt=0
dE_Hdt=dUHdt*p_Am
}else{
  U_H_pres=E_H_pres/p_Am    
dUHdt=(1-kappa)*Sc-k_J*U_H_pres
dE_Hdt=dUHdt*p_Am
}
}else{
  dE_Hdt = (1-kappa)*p_C-p_J
}
}else{
  dE_Hdt = 0
}


if(E_H_pres>=E_Hp){
p_D = p_M+p_J+(1-k_R)*p_R
}else{
  p_D = p_M+p_J+p_R
}

p_G = p_C-p_M-p_J-p_R


if(hour==1){
E_H = E_H_init + dE_Hdt
}else{
  E_H = E_H+dE_Hdt
}

#    aging
dqdt = (q_pres*(V_pres/V_max)*s_G+h_a)*(E_pres/E_m)*((vdot/L_pres)-rdot)-rdot*q_pres

if(E_H_pres>E_Hb){
if(hour==1){
q = q_init + dqdt
}else{
  q = q+dqdt
}
}else{
  q = 0
}

#    dhsds = h_a*q(hour)/V_pres
dhsds = q_pres-rdot*hs_pres

if(E_H_pres>E_Hb){
if(hour==1){
hs = hs_init + dhsds
}else{
  hs = hs+dhsds
}
}else{
  hs = 0
}

h_w = ((h_a*(E_pres/E_m)*vdot)/(6*V_pres^(1./3.)))^(1./3.)
dsurvdt = -1*surviv_pres*hs
surviv = surviv_pres+dsurvdt

#     accumulate energy/matter in reproduction buffer
#    if it is the beginning of the day
if(E_H_pres>E_Hp){
if(hour==1){
#      if the buffer ran out in the previous hour
if(cumrepro_init<0){
#       keep it empty
cumrepro=0
}else{
  cumrepro = cumrepro_init
}
}else{
  #     it is not the first first hour and it is not the first day
#      if the buffer ran out in the previous hour
if(cumrepro<0){
#       keep it empty
cumrepro=0
}else{
  #       otherwise start it filling up according to p_R but subtract anything that goes to the batch
cumrepro = cumrepro_prev+p_R*k_R-p_B_past
}
}
}

#     accumulate energy/matter in egg batch buffer
#    if it is the beginning of the day
if(hour==1){
#     { if it is the first day of the simulation
#     if(day+365*(iyear-1)==1){
#      nothing in the buffer yet      
cumbatch = cumbatch_init
}else{
  #     it is not the first first hour of the first day
#       otherwise start it filling up 
cumbatch = cumbatch_prev+p_B
}

if(stage==2){
if(cumbatch<0.1*clutchenergy){
stage=3
}
}

V=V_pres+dVdt

if(V<0){
V=0
}

ED = E_pres+dEdt
#    make sure ED doesn't go below zero
if(ED<0){
ED=0
}


#    svl in mm
svl = V^(0.3333333333333)/delta_deb*10

if(E_H<=E_Hb){
stage=0
}else{
if(E_H<E_Hj){
stage=1
}else{
if(E_H<E_Hp){
stage=2
}else{
stage=3
}
}
}

if(cumbatch>0){
if(E_H>E_Hp){
stage=4
}else{
stage=stage
}
}

if((cumbatch>clutchenergy) | (pregnant==1)){
#       for variable clutch size from repro and batch buffers
#    if((cumbatch(hour)>clutchenergy) | (pregnant==1).or
#     &.((viviparous==1) & (cumbatch(hour)+cumrepro(hour)>
#     &clutchenergy))){
#     batch is ready so if viviparous, start gestation, }else{ dump it
if(viviparous==1){
if((pregnant==0) & (breeding==1)){
v_baby=v_init_baby
e_baby=e_init_baby
EH_baby=0.
pregnant=1
testclutch=floor(cumbatch/E_egg)
#       for variable clutch size from repro and batch buffers
#       testclutch=floor((cumbatch(hour)+cumrepro(hour))/E_egg)
#       testclutch=real(testclutch)
if(testclutch>clutchsize){
clutchsize=testclutch
clutchenergy = E_egg*clutchsize
}
#       for variable clutch size from repro and batch buffers
if(cumbatch<clutchenergy){
#        needs to draw from repro buffer - temporarily store current repro as cumrepro_temp, 
#        { remove what is needed from the repro buffer and add it to the batch buffer
cumrepro_temp=cumrepro
cumrepro=cumrepro+cumbatch-clutchenergy
cumbatch=cumbatch+cumrepro_temp-cumrepro
}
}
if(hour==1){
v_baby=v_baby_init
e_baby=e_baby_init
EH_baby=EH_baby_init
}
#if(pregnant==1){
#call deb_baby
#}
if(EH_baby>E_Hb){

cumbatch(hour) = cumbatch(hour)-clutchenergy
repro(hour)=1
pregnant=0
v_baby=v_init_baby
e_baby=e_init_baby
EH_baby=0
newclutch=clutchsize
         fecundity=fecundity+clutchsize
         clutches=clutches+1
         pregnant=0
        }
       }else{
#      not viviparous, so lay the eggs at next period of activity
        if(breedrainthresh>0){
         if(rainfall<breedrainthresh){
          #goto 898
         }
        }

         testclutch=floor(cumbatch/E_egg)
        if(testclutch>clutchsize){
         clutchsize=testclutch
        }
        cumbatch = cumbatch-clutchenergy
        repro=1
        fecundity=fecundity+clutchsize
        clutches=clutches+1
       }
      }



gutfull=ms_pres/(MsM*V_pres)  # Gut model
if(gutfull>1){
gutfull=1
}

if(E_H_pres>E_Hb){
if(acthr > 0){
# Regulates X_food dynamics
dMsdt = p_Xm*(X_food/(halfsat+X_food))*V_pres^(2./3.)*funct-1.*(p_Am/kappa_X)*V_pres^(2./3.)
}else{
dMsdt = -1.*(p_Am/kappa_X)*V_pres^(2./3.)
}
}else{
dMsdt = -1.*(p_Am/kappa_X)*V_pres^(2./3.)
}

if(V_pres==0){
dMsdt=0
}


ms = ms_init+dMsdt

if(ms<0){
ms=0
}

if(ms>MsM*V_pres){
ms=MsM*V_pres
}

gutfull=ms/(MsM*V_pres)


ms_past=ms    
p_B_past=p_B

#************    mass balances     *******************

JOJx=p_A*etaO[1,1]+p_D*etaO[1,2]+p_G*etaO[1,3]
JOJv=p_A*etaO[2,1]+p_D*etaO[2,2]+p_G*etaO[2,3]
JOJe=p_A*etaO[3,1]+p_D*etaO[3,2]+p_G*etaO[3,3]
JOJp=p_A*etaO[4,1]+p_D*etaO[4,2]+p_G*etaO[4,3]

JOJx_GM=p_D*etaO[1,2]+p_G*etaO[1,3]
JOJv_GM=p_D*etaO[2,2]+p_G*etaO[2,3]
JOJe_GM=p_D*etaO[3,2]+p_G*etaO[3,3]
JOJp_GM=p_D*etaO[4,2]+p_G*etaO[4,3]

JMCO2=JOJx*JM_JO[1,1]+JOJv*JM_JO[1,2]+JOJe*JM_JO[1,3]+JOJp*JM_JO[1,4]
JMH2O=JOJx*JM_JO[2,1]+JOJv*JM_JO[2,2]+JOJe*JM_JO[2,3]+JOJp*JM_JO[2,4]
JMO2=JOJx*JM_JO[3,1]+JOJv*JM_JO[3,2]+JOJe*JM_JO[3,3]+JOJp*JM_JO[3,4]
JMNWASTE=JOJx*JM_JO[4,1]+JOJv*JM_JO[4,2]+JOJe*JM_JO[4,3]+JOJp*JM_JO[4,4]

JMCO2_GM=JOJx_GM*JM_JO[1,1]+JOJv_GM*JM_JO[1,2]+JOJe_GM*JM_JO[1,3]+JOJp_GM*JM_JO[1,4]
JMH2O_GM=JOJx_GM*JM_JO[2,1]+JOJv_GM*JM_JO[2,2]+JOJe_GM*JM_JO[2,3]+JOJp_GM*JM_JO[2,4]
JMO2_GM=JOJx_GM*JM_JO[3,1]+JOJv_GM*JM_JO[3,2]+JOJe_GM*JM_JO[3,3]+JOJp_GM*JM_JO[3,4]
JMNWASTE_GM=JOJx_GM*JM_JO[4,1]+JOJv_GM*JM_JO[4,2]+JOJe_GM*JM_JO[4,3]+JOJp_GM*JM_JO[4,4]


#    mlO2/h, temperature corrected (including SDA)


O2FLUX = -1*JMO2/(T_ref/Tbody/24.4)*1000



CO2FLUX = JMCO2/(T_ref/Tbody/24.4)*1000
#    mlO2/h, stp
#    MLO2(hour) = -1*JMO2/(T_ref/Tbody/24.4)*1000
MLO2 = (-1*JMO2*(0.082058*(Tbody+273.15))/(0.082058*293.15))*24.06*1000
#    g metabolic water/h
GH2OMET = JMH2O*18.01528
#    metabolic heat production (Watts) - growth overhead plus dissipation power (maintenance, maturity maintenance, 
#    maturation/repro overheads) plus assimilation overheads - correct to 20 degrees so it can be temperature corrected
#    in MET.f for the new guessed Tbody
DEBQMET = ((1-kappa_G)*p_G+p_D+(p_X-p_A-p_A*mu_P*eta_PA))/3600/Tcorr

DRYFOOD=-1*JOJx*w_X
FAECES=JOJp*w_P
NWASTE=JMNWASTE*w_N
if(pregnant==1){
wetgonad = ((cumrepro/mu_E)*w_E)/eggdryfrac+((((v_baby*e_baby)/mu_E)*w_E)/d_V + v_baby)*clutchsize
}else{
wetgonad = ((cumrepro/mu_E)*w_E)/eggdryfrac+((cumbatch/mu_E)*w_E)/eggdryfrac
}
wetstorage = ((V*ED/mu_E)*w_E)/d_V
#    wetfood(hour) = ((ms(hour)/mu_E)*w_E)/d_V
wetfood = ms/21525.37/(1.-0.18)
wetmass = V*ANDENS_deb+wetgonad+wetstorage+wetfood
gutfreemass=V*ANDENS_deb+wetgonad+wetstorage
potfreemass=V*ANDENS_deb+(((V*E_m)/mu_E)*w_E)/d_V # this is the max potential mass if reserve density is at max value

# final results to go to next iteration

E_init<-ED
v_init<-V
E_H_init<-E_H
q_init<-q
hs_init<-hs
surviv_init<-surviv
ms_init<-ms


cumrepro_prev<-cumrepro
cumbatch_prev<-cumbatch

results_deb<-c(hour,E_init,v_init,E_H_init,q_init,hs_init,surviv_init,ms_pres,cumrepro,cumbatch,O2FLUX,CO2FLUX,MLO2,GH2OMET,DEBQMET,DRYFOOD,FAECES,NWASTE,wetgonad,wetstorage,wetfood,wetmass,gutfreemass)
  all_results_deb<-results_deb

} # end DEB function




