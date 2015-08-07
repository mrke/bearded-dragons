# Transient model simple environmental transition
source('OneLumpAnalytical_bearded.R')

# constants
cp<-4185 #specific heat of flesh, J/kg-C
emis<-0.95 #emissivity of skin, -
sigma<-0.0000000567 #Stefan-Boltzman, W/mK
Fo_e<-0.8 #config factor, object to IR environment, -
rho<-1000 #animal density, kg/m3
abs<-0.92 #animal solar absorptivity
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
posture<-'n' # pointing normal 'n' or parallel 'p' to the sun's rays?
FATOSK<-0.4 # configuration factor to sky
FATOSB<-0.4 # configuration factor to substrate

press<-101325 #atmospheric pressure, pa
sub_reflect<-0.2 # solar reflectance of substrate
pctdif<-0.1 # proportion of solar energy that is diffuse (rather than direct beam)

absorbs<-c(0.77,0.92)

data<-read.csv('environmental data.csv',stringsAsFactors=FALSE) # read in environmental data

for(j in 1:2){ # loop through each posture
for(k in 1:1){ # loop through each initial temperature
for(m in 1:2){ # loop through each absorptivity
for(i in 1:nrow(data)){ # loop through each set of environmental conditions

  if(m==1){abs<-absorbs[1]}else{abs<-absorbs[2]} #choose absorptivity, depending on the loop through m
  if(j==1){posture<-'n'}else{posture<-'p'} #choose the posture, depending on the loop through j
  if(k==1){Tc_init<-10}else{Tc_init<-40} #choose the initial temperature, depending on the loop through k

  
  # environment
  Qsol<-data[i,1] #solar radiation, W/m2
  Zen<-data[i,2] #zenith angle of sun (90 is below horizon), degrees
  vel<-data[i,3] #wind speed, m/s
  Tair<-data[i,4] #air temperature, C
  Tsurf<-data[i,5] #substrate temperature, C
  
  #get sky temp
  cloud<-data[i,7]
  shade<-data[i,8]
  RH<-data[i,9]
  VIEWF<-1
  QRADHL<-0
  TMAXK<-Tair+273.15
  loge<-TMAXK
  loge[loge>273.16]<- -7.90298*(373.16/TMAXK-1.)+5.02808*log10(373.16/TMAXK)-1.3816E-07*(10.^(11.344*(1.-TMAXK/373.16))-1.)+8.1328E-03*(10.^(-3.49149*(373.16/TMAXK-1.))-1.)+log10(1013.246)
  loge[loge<=273.16]<- -9.09718*(273.16/TMAXK-1.)-3.56654*log10(273.16/TMAXK)+.876793*(1.-TMAXK/273.16)+log10(6.1071)
  estar<-(10.^loge)*100.
  VAPRES<-RH/100*estar
  
  ARAD<-1.72*((VAPRES/1000.)/(Tair+273.16))^(1./7.)*0.0000000567*(Tair+273.16)**4*60./(4.185*10000.)
  
  SIGP<-.8126E-10
  SLEP<-1.
  #    APPROXIMATING CLOUD RADIANT TEMPERATURE AS 2 M SHADE TEMPERATURE 
  CRAD<-SIGP*SLEP*(Tair+273.)^4
  
  #    GROUND SURFACE RADIATION TEMPERATURE      
  SRAD<-SIGP*SLEP*(Tsurf+273.)^4
  CLR<-1.- (cloud/100.)
  CLEAR<-ARAD*CLR
  cloud<-CRAD*(cloud/100.)
  
  QRADSK<-(CLEAR + cloud)*((100.- shade)/100.)
  
  QRADVG<-(shade/100.)*CRAD
  
  QRADGR<-((100.-shade)/100.)*SRAD+(shade/100.)*CRAD
  QRAD <- (QRADSK + QRADVG)*VIEWF + QRADHL*(1-VIEWF) - QRADGR
  TSKY<-((QRAD+QRADGR)/(SIGP))**(1./4.)-273
  
  Trad<-mean(c(Tsurf,TSKY))
  
  
  mass<-500
  q<-0
  kflesh<-0.5
  input<-c(kflesh,q,cp,emis,sigma,Fo_e,rho,abs,lometry,customallom,shape_a,shape_b,shape_c,posture,FATOSK,FATOSB,mass,sub_reflect,pctdif,Qsol,vel,Tair,Trad,Zen)
  
  #Tbs<-as.data.frame(ode(y=Tc_init,times=seq(0,3600*2,1),func=transient,parms=input))
  times=seq(0,3600,10)
  thresh<-35
  results<-onelump(times, Tc_init, thresh, input)
  Tbs<-as.data.frame(cbind(times,results$Tc))
  colnames(Tbs)<-c('time','Tb')
  dTbs<-as.data.frame(cbind(times,results$dTc*60)) # rates of change in deg C per min
  colnames(dTbs)<-c('time','dTb')
  final_temp<-results$Tcf
  time.to.35<-results$timethresh
  tau<-results$tau
  Tbs$time<-Tbs$time/60 #convert to minutes
  rate<-if(Tc_init==10){max(dTbs$dTb)}else{min(dTbs$dTb)}
  #rate<-mean(dTbs$dTb)
  if(i==1){
  with(Tbs,{plot(Tb~time,type='l',col=i,ylim=c(0,40),main=paste('abs =',abs,' pos = ',posture,sep=""))})
  }else{
  with(Tbs,{points(Tb~time,type='l',col=i,ylim=c(0,40),main=paste('abs =',abs,' pos = ',posture,sep=""))})
  }

  #with(dTbs,{plot(dTb~time,type='l',col='red')})
  if(i==1 & j==1 & k==1 & m==1){
    Tpred<-cbind.data.frame(abs,final_temp,posture,Tc_init,time.to.35,tau,rate) # need to bind the predicted temperature, time to thresh and tau with the selected posture and initial temp, turn it into a data frame and then transpose it so it is a row
  }else{
    Tpred<-rbind.data.frame(Tpred,cbind.data.frame(abs,final_temp,posture,Tc_init,time.to.35,tau,rate)) # row bind it to the previous result
  }
  cat(i,'\n')
} #end loop through observations in 'data'
} #end loop through absorptivities  
} #end loop through initial temps
} #end loop through postures

data_final<-cbind.data.frame(rep(seq(8,16),4),Tpred,rbind.data.frame(data,data,data,data)) # put the final set of results (Tb and posture) together with the environmental data (but need two sets of environmental data, one for each posture)
colnames(data_final)[1:8]<-c('hour','absorb','Tcfinal','posture','Tc_init','Time.to.35','tau','rate')

plot(rate~hour,data=subset(data_final,absorb==0.92 & Tc_init==10 & posture=='n'),type='b',col='black',main='warming',ylab='max rate, C/min')
points(rate~hour,data=subset(data_final,absorb==0.77 & Tc_init==10 & posture=='n'),type='b',col='grey')
#plot(rate~hour,data=subset(data_final,absorb==0.92 & Tc_init==40 & posture=='n'),type='b',col='black',main='cooling',ylab='max rate, C/min')
#points(rate~hour,data=subset(data_final,absorb==0.77 & Tc_init==40 & posture=='n'),type='b',col='grey')

plot(Tcfinal~hour,data=subset(data_final,absorb==0.92 & Tc_init==10 & posture=='n'),type='b',col='black',ylab='Steady State Tb, C')
points(Tcfinal~hour,data=subset(data_final,absorb==0.77 & Tc_init==10 & posture=='n'),type='b',col='grey')


write.csv(data_final,'results.csv')


