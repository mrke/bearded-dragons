# Transient model function set up for a constant environment and a starting condition
# Michael Kearney & Warren Porter developed this R function on 28 July 2014.
# Michael modified it on 1st Aug to work as an analytical model, without evaporation
onelump<-function(t,y,thresh,input){
  unlist(input)
  sigma<-0.0000000567 #Stefan-Boltzman, W/(m.K)
  Zenith<-Zen*pi/180 # zenith angle in radians
  Tc<-y # core temperature, deg C
  Tskin<-y+0.1 # make skin temperature very close to core temperature
  if(vel<0.01){
    vel<-0.01 # don't let wind speed go too low - always some free convection
  }
  S2<-0.0001 # initializing
  DENSTY<-press/(287.04*(Tair+273)) # air density, kg/m3
  THCOND<-0.02425+(7.038*10^-5*Tair) # air thermal conductivity, W/(m.K)
  VISDYN<-(1.8325*10^-5*((296.16+120)/((Tair+273)+120)))*(((Tair+273)/296.16)^1.5) # dynamic viscosity of air, kg/(m.s)

  m<-mass/1000 # convert mass to kg
  C<-m*cp # thermal capacitance, J/K
  V<-m/rho # volume, m3
  Qgen<-q*V # total metabolic heat, J
  L<-V^(1./3.) # characteristic dimension, m
  
  # geometry section ############################################################
  # FLAT PLATE geometry
  if(lometry==0){
    ALENTH<-(V/shape_b*shape_c)^(1./3.) # length, m
    AWIDTH<-ALENTH*shape_b # width, m
    AHEIT<-ALENTH*shape_c # height, m    
    ATOT<-ALENTH*AWIDTH*2.+ALENTH*AHEIT*2.+AWIDTH*AHEIT*2. # total area, m2
    ASILN<-ALENTH*AWIDTH # max silhouette area, m2
    ASILP<-AWIDTH*AHEIT # min silhouette area, m2
    L<-AHEIT # characteristic dimension, m
    if(AWIDTH<=ALENTH){
      L<-AWIDTH
    }else{
      L<-ALENTH
    }
    R<-ALENTH/2. # 'radius', m
  } 
  
  # CYLINDER geometry      
  if(lometry==1){
    R1<-(V/(pi*shape_b*2))^(1./3.) # radius, m
    ALENTH<-2*R1*shape_b # length, m
    ATOT<-2*pi*R1^2+2*pi*R1*ALENTH # total surface area, m2
    AWIDTH<-2.*R1 # width, m
    ASILN<-AWIDTH*ALENTH # max silhouette area, m2
    ASILP<-pi*R1^2 # min silhouette area, m2
    L<-ALENTH # characteristic dimension, m
    R2<-L/2
    if(R1>R2){ # choose shortest dimension as R
      R<-R2
    }else{
      R<-R1
    }
  }
  
  # Ellipsoid geometry
  if(lometry==2){
    A1<-((3./4.)*V/(pi*shape_b*shape_c))^0.333 # axis A, m  
    B1<-A1*shape_b # axis B, m
    C1<-A1*shape_c # axis C, m
    P1<-1.6075 # a constant
    ATOT<-(4*pi*(((A1^P1*B1^P1+A1^P1*C1^P1+B1^P1*C1^P1))/3)^(1/P1)) # total surface area, m2
    ASILN<-max(pi*A1*C1,pi*B1*C1) # max silhouette area, m2
    ASILP<-min(pi*A1*C1,pi*B1*C1) # min silhouette area, m2
    S2<-(A1^2*B1^2*C1^2)/(A1^2*B1^2+A1^2*C1^2+B1^2*C1^2) # fraction of semi-major and minor axes, see Porter and Kearney 2009 supp1
    kflesh<-0.5 + 6.14*B1 + 0.439 # thermal conductivity of flesh as a function of radius, see Porter and Kearney 2009
  }              
  
  # Lizard geometry - DESERT IGUANA (PORTER ET AL. 1973 OECOLOGIA)
  if(lometry==3){
    ATOT<-(10.4713*mass^.688)/10000. # total surface area, m2
    AV<-(0.425*mass^.85)/10000. # ventral surface area, m2   
    # NORMAL AND POINTING @ SUN SILHOUETTE AREA: PORTER & TRACY 1984   
    ASILN<-(3.798*mass^.683)/10000. # Max. silhouette area (normal to the sun), m2
    ASILP<-(0.694*mass^.743)/10000. # Min. silhouette area (pointing toward the sun), m2
    R<-L
    dors_vent_ratio<-AV/ATOT
    ATOT<-595.25/10000 #using ratio for Dipsosaurus to get ventral from total
    AV<-dors_vent_ratio*ATOT# 283.75/10000.   
    # NORMAL AND POINTING @ SUN SILHOUETTE AREA: PORTER & TRACY 1984   
    # Max. silhouette area (normal to the sun)
    ASILN<-142.25/10000.
    # Min. silhouette area (pointing toward the sun)         
    ASILP<-85.25/10000.
    R<-L
  }
  
  # Frog geometry - LEOPARD FROG (C.R. TRACY 1976 ECOL. MONOG.)
  if(lometry==4){
    ATOT = (12.79*mass^.606)/10000. # total surface area, m2
    AV = (0.425*mass^.85)/10000. # ventral surface area, m2  
    # NORMAL AND POINTING @ SUN SILHOUETTE AREA: EQ'N 11 TRACY 1976
    ZEN<-0.
    PCTN<-1.38171E-06*ZEN^4-1.93335E-04*ZEN^3+4.75761E-03*ZEN^2-0.167912*ZEN+45.8228  
    ASILN<-PCTN*ATOT/100. # Max. silhouette area (normal to the sun), m2
    ZEN<-90. 
    PCTP<-1.38171E-06*ZEN^4-1.93335E-04*ZEN^3+4.75761E-03*ZEN^2-0.167912*ZEN+45.8228  
    ASILP<-PCTP*ATOT/100. # Min. silhouette area (pointing toward the sun), m2
    R<-L
  }
  
  # user defined geometry
  if(lometry==5){
    ATOT = (customallom[1]*mass^customallom[2])/10000. # total surface area, m2 
    AV = (customallom[3]*mass^customallom[4])/10000. # ventral surface area, m2   
    # NORMAL AND POINTING @ SUN SILHOUETTE AREA: PORTER & TRACY 1984   
    # User must define Max. silhouette area (normal to the sun)
    ASILN = (customallom[5]*mass^customallom[6])/10000. # Max. silhouette area (normal to the sun), m2
    # User must define Min. silhouette area (pointing toward the sun)         
    ASILP = (customallom[7]*mass^customallom[8])/10000. # Min. silhouette area (pointing toward the sun), m2
    R<-L
  }
  # end geometry section ############################################################
  
  if(Zen>=90){
    Qnorm<-0
  }else{
    Qnorm <- (Qsol / cos(Zenith)) 
  }
  if(Qnorm>1367){
    Qnorm<-1367 #making sure that low sun angles don't lead to solar values greater than the solar constant
  }
  if(posture=='p'){
    Qabs<-(Qnorm*(1-pctdif)*ASILP+Qsol*pctdif*FATOSK*ATOT+Qsol*sub_reflect*FATOSB*ATOT)*abs
  }
  if(posture=='n'){
    Qabs<-(Qnorm*(1-pctdif)*ASILN+Qsol*pctdif*FATOSK*ATOT+Qsol*sub_reflect*FATOSB*ATOT)*abs
  }
  if(posture=='b'){
    Qabs<-(Qnorm*(1-pctdif)*(ASILN+ASILP)/2+Qsol*pctdif*FATOSK*ATOT+Qsol*sub_reflect*FATOSB*ATOT)*abs
  }
  
  Rrad<-((Tskin+273)-(Trad+273))/(emis*sigma*Fo_e*ATOT*((Tskin+273)^4-(Trad+273)^4)) # radiation resistance
  Re<-DENSTY*vel*L/VISDYN # Reynolds number
  PR<-1005.7*VISDYN/THCOND # Prandlt number
  
  if(lometry==0){
    NUfor<-0.102*Re^0.675*PR^(1./3.)
  }
  if(lometry==3|lometry ==5){
    NUfor<-0.35*Re^0.6
  }
  if(lometry==1){ 
    #       FORCED CONVECTION OF A CYLINDER    
    #       ADJUSTING NU - RE CORRELATION FOR RE NUMBER (P. 260 MCADAMS,1954) 
    if(Re<4.){ 
      NUfor=.891*Re**.33  
    }else{
      if(Re<40.){
        NUfor=.821*Re**.385 
      }else{
        if(Re<4000.){
          NUfor=.615*Re**.466 
        }else{ 
          if(Re<40000.){  
            NUfor=.174*Re**.618 
          }else{
            if(Re<400000.){ 
              NUfor=.0239*Re**.805
            }else{
              NUfor=.0239*Re**.805
            }}}}} 
  }
  if(lometry==2|lometry==4){
    NUfor<-0.35*Re^(0.6) # Nusselt number, forced convection
  }
  hc_forced<-NUfor*THCOND/L # convection coefficent, forced
  
  GR<-abs(DENSTY^2*(1/(Tair+273.15))*9.80665*L^3*(Tskin-Tair)/VISDYN^2) # Grashof number
  Raylei<-GR*PR # Rayleigh number
  
  # get Nusselt for Free Convect
  if(lometry==0){
    NUfre=0.55*Raylei^0.25
  }
  if(lometry==1|lometry==3|lometry==5){
    if(Raylei<1.0e-05){  
      NUfre=0.4
    }else{
      if(Raylei<0.1){   
        NUfre=0.976*Raylei^0.0784
      }else{
        if(Raylei<100){  
          NUfre=1.1173*Raylei^0.1344   
        }else{  
          if(Raylei<10000.){  
            NUfre=0.7455*Raylei^0.2167   
          }else{  
            if(Raylei<1.0e+09){  
              NUfre=0.5168*Raylei^0.2501
            }else{ 
              if(Raylei<1.0e+12){  
                NUfre=0.5168*Raylei^0.2501
              }}}}}}
  }
  
  if(lometry==2|lometry==4){
    Raylei=(GR^.25)*(PR^.333)
    NUfre=2.+0.60*Raylei
  }
  hc_free<-NUfre*THCOND/L # convection coefficent, forced
  hc_comb<-hc_free+hc_forced
  Rconv<-1/(hc_comb*ATOT)
  Nu<-hc_comb*L/THCOND # Nu combined
  hr<-4*emis*sigma*((Tc+Trad)/2+273)^3 # radiation resistance
  hc<-hc_comb
  
  if(lometry==2){
    j<-(Qabs+Qgen+hc*ATOT*((q*S2)/(2*kflesh)+Tair)+hr*ATOT*((q*S2)/(2*kflesh)+Trad))/C
  }else{
    j<-(Qabs+Qgen+hc*ATOT*((q*R^2)/(2*kflesh)+Tair)+hr*ATOT*((q*S2)/(2*kflesh)+Trad))/C
  }
  kTc<-ATOT*(Tc*hc+Tc*hr)/C
  k<-ATOT*(hc+hr)/C
  Tcf<-j/k # final Tc
  Tci<-Tc
  Tc<-(Tci-Tcf)*exp(-1*k*t)+Tcf # Tc at time t
  timethresh<-log((thresh-Tcf)/(Tci-Tcf))/(-1*k)
  tau<-(rho*V*cp)/(ATOT*(hc+hr)) # time constant
  dTc<-j-k*Tc
  return(list(Tc=Tc,Tcf=Tcf,tau=tau,dTc=dTc,timethresh=timethresh))
}
