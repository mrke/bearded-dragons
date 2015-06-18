library(raster)
library(ncdf)

setwd('c:/NicheMapR_Working/projects/bearded dragons/landscape sims/')
longlats<-read.csv('longlats.csv')[,2:3]
files<-list.files()
files1997<-files[grep(files,pattern = '1997')]
files1998<-files[grep(files,pattern = '1998')]
files1999<-files[grep(files,pattern = '1999')]
files2000<-files[grep(files,pattern = '2000')]
files2001<-files[grep(files,pattern = '2001')]

colnames<-c("n","long","lat","date","doy","maxrun","sumact","bouts","mornbask","mornfor","mid1","meanmid","arvo")

quadres<-0.6 # quadrangle resolution (degrees)
lat1<-min(longlats[,2])-quadres/2 # min latitude
lat2<-max(longlats[,2])+quadres/2 # max latitude
lon1<-min(longlats[,1])-quadres/2 # min longitude
lon2<-max(longlats[,1])+quadres/2 # max longitude
quadwid<-(lon2-lon1)/quadres
quadlen<-(lat2-lat1)/quadres
gridout <- raster(ncol=quadwid, nrow=quadlen, xmn=lon1, xmx=lon2, ymn=lat1, ymx=lat2)



sims<-c("var72_85","con72","con78","con85")
years<-seq(1997,2001,1)
tzone<-paste("Etc/GMT+",10,sep="")

# maxrun
for(m in 1:length(years)){
dates<-seq(ISOdate(years[m],1,1,tz=tzone)-3600*12, ISOdate(years[m]+1,1,1,tz=tzone)-3600*13, by="days")
dates<-subset(dates, format(dates, "%m/%d")!= "02/29") # remove leap years
for(i in 1:eval(parse(text=paste("length(files",years[m],")",sep='')))){
  data<-read.csv(eval(parse(text=paste("files",years[m],"[",i,"]",sep=""))),header=FALSE)
  colnames(data)<-colnames
  data[,5]<-data[,5]-data[1,5]+1 # get correct day of year
  for(j in 1:365){
  daily_sub<-subset(data,doy==j)
  x<-cbind(daily_sub[,2],daily_sub[,3])
  grid <- rasterize(x, gridout, daily_sub$maxrun)
  #plot(grid,main=j)
  if(j==1){s<-grid}else{s<-stack(s,grid)}
  cat(j,'of 365 ', years[m],'sim ',sims[i], '\n')  
  }
  
  filename<-paste("maxrun_",years[m],"_",sims[i],".nc",sep="")
  writeRaster(s, filename=filename, overwrite=TRUE)
  temp<-open.ncdf( filename, write=TRUE, readunlim=TRUE)
  put.var.ncdf( temp, varid='z', vals=dates)
  close.ncdf(temp)  
  # yearly data
  } #end loop through sims
} #end loop through years

#sumact
for(m in 1:length(years)){
dates<-seq(ISOdate(years[m],1,1,tz=tzone)-3600*12, ISOdate(years[m]+1,1,1,tz=tzone)-3600*13, by="days")
dates<-subset(dates, format(dates, "%m/%d")!= "02/29") # remove leap years
for(i in 1:eval(parse(text=paste("length(files",years[m],")",sep='')))){
  data<-read.csv(eval(parse(text=paste("files",years[m],"[",i,"]",sep=""))),header=FALSE)
  colnames(data)<-colnames
  data[,5]<-data[,5]-data[1,5]+1 # get correct day of year
  for(j in 1:365){
  daily_sub<-subset(data,doy==j)
  x<-cbind(daily_sub[,2],daily_sub[,3])
  grid <- rasterize(x, gridout, daily_sub$sumact)
  #plot(grid,main=j)
  if(j==1){s<-grid}else{s<-stack(s,grid)}
  cat(j,'of 365 ', years[m],'sim ',sims[i], '\n')  
  }
  
  filename<-paste("sumact_",years[m],"_",sims[i],".nc",sep="")
  writeRaster(s, filename=filename, overwrite=TRUE)
  temp<-open.ncdf( filename, write=TRUE, readunlim=TRUE)
  put.var.ncdf( temp, varid='z', vals=dates)
  close.ncdf(temp)  
  # yearly data
  } #end loop through sims
} #end loop through years

#bouts
for(m in 1:length(years)){
dates<-seq(ISOdate(years[m],1,1,tz=tzone)-3600*12, ISOdate(years[m]+1,1,1,tz=tzone)-3600*13, by="days")
dates<-subset(dates, format(dates, "%m/%d")!= "02/29") # remove leap years
for(i in 1:eval(parse(text=paste("length(files",years[m],")",sep='')))){
  data<-read.csv(eval(parse(text=paste("files",years[m],"[",i,"]",sep=""))),header=FALSE)
  colnames(data)<-colnames
  data[,5]<-data[,5]-data[1,5]+1 # get correct day of year
  for(j in 1:365){
  daily_sub<-subset(data,doy==j)
  x<-cbind(daily_sub[,2],daily_sub[,3])
  grid <- rasterize(x, gridout, daily_sub$bouts)
  #plot(grid,main=j)
  if(j==1){s<-grid}else{s<-stack(s,grid)}
  cat(j,'of 365 ', years[m],'sim ',sims[i], '\n')  
  }
  
  filename<-paste("bouts_",years[m],"_",sims[i],".nc",sep="")
  writeRaster(s, filename=filename, overwrite=TRUE)
  temp<-open.ncdf( filename, write=TRUE, readunlim=TRUE)
  put.var.ncdf( temp, varid='z', vals=dates)
  close.ncdf(temp)  
  # yearly data
  } #end loop through sims
} #end loop through years

#mornbask
for(m in 1:length(years)){
dates<-seq(ISOdate(years[m],1,1,tz=tzone)-3600*12, ISOdate(years[m]+1,1,1,tz=tzone)-3600*13, by="days")
dates<-subset(dates, format(dates, "%m/%d")!= "02/29") # remove leap years
for(i in 1:eval(parse(text=paste("length(files",years[m],")",sep='')))){
  data<-read.csv(eval(parse(text=paste("files",years[m],"[",i,"]",sep=""))),header=FALSE)
  colnames(data)<-colnames
  data[,5]<-data[,5]-data[1,5]+1 # get correct day of year
  for(j in 1:365){
  daily_sub<-subset(data,doy==j)
  x<-cbind(daily_sub[,2],daily_sub[,3])
  grid <- rasterize(x, gridout, daily_sub$mornbask)
  #plot(grid,main=j)
  if(j==1){s<-grid}else{s<-stack(s,grid)}
  cat(j,'of 365 ', years[m],'sim ',sims[i], '\n')  
  }
  
  filename<-paste("mornbask_",years[m],"_",sims[i],".nc",sep="")
  writeRaster(s, filename=filename, overwrite=TRUE)
  temp<-open.ncdf( filename, write=TRUE, readunlim=TRUE)
  put.var.ncdf( temp, varid='z', vals=dates)
  close.ncdf(temp)  
  # yearly data
  } #end loop through sims
} #end loop through years


#mornfor
for(m in 1:length(years)){
dates<-seq(ISOdate(years[m],1,1,tz=tzone)-3600*12, ISOdate(years[m]+1,1,1,tz=tzone)-3600*13, by="days")
dates<-subset(dates, format(dates, "%m/%d")!= "02/29") # remove leap years
for(i in 1:eval(parse(text=paste("length(files",years[m],")",sep='')))){
  data<-read.csv(eval(parse(text=paste("files",years[m],"[",i,"]",sep=""))),header=FALSE)
  colnames(data)<-colnames
  data[,5]<-data[,5]-data[1,5]+1 # get correct day of year
  for(j in 1:365){
  daily_sub<-subset(data,doy==j)
  x<-cbind(daily_sub[,2],daily_sub[,3])
  grid <- rasterize(x, gridout, daily_sub$mornfor)
  #plot(grid,main=j)
  if(j==1){s<-grid}else{s<-stack(s,grid)}
  cat(j,'of 365 ', years[m],'sim ',sims[i], '\n')  
  }
  
  filename<-paste("mornfor_",years[m],"_",sims[i],".nc",sep="")
  writeRaster(s, filename=filename, overwrite=TRUE)
  temp<-open.ncdf( filename, write=TRUE, readunlim=TRUE)
  put.var.ncdf( temp, varid='z', vals=dates)
  close.ncdf(temp)  
  # yearly data
  } #end loop through sims
} #end loop through years

#mid1
for(m in 1:length(years)){
dates<-seq(ISOdate(years[m],1,1,tz=tzone)-3600*12, ISOdate(years[m]+1,1,1,tz=tzone)-3600*13, by="days")
dates<-subset(dates, format(dates, "%m/%d")!= "02/29") # remove leap years
for(i in 1:eval(parse(text=paste("length(files",years[m],")",sep='')))){
  data<-read.csv(eval(parse(text=paste("files",years[m],"[",i,"]",sep=""))),header=FALSE)
  colnames(data)<-colnames
  data[,5]<-data[,5]-data[1,5]+1 # get correct day of year
  for(j in 1:365){
  daily_sub<-subset(data,doy==j)
  x<-cbind(daily_sub[,2],daily_sub[,3])
  grid <- rasterize(x, gridout, daily_sub$mid1)
  #plot(grid,main=j)
  if(j==1){s<-grid}else{s<-stack(s,grid)}
  cat(j,'of 365 ', years[m],'sim ',sims[i], '\n')  
  }
  
  filename<-paste("mid1_",years[m],"_",sims[i],".nc",sep="")
  writeRaster(s, filename=filename, overwrite=TRUE)
  temp<-open.ncdf( filename, write=TRUE, readunlim=TRUE)
  put.var.ncdf( temp, varid='z', vals=dates)
  close.ncdf(temp)  
  # yearly data
  } #end loop through sims
} #end loop through years

#meanmid
for(m in 1:length(years)){
dates<-seq(ISOdate(years[m],1,1,tz=tzone)-3600*12, ISOdate(years[m]+1,1,1,tz=tzone)-3600*13, by="days")
dates<-subset(dates, format(dates, "%m/%d")!= "02/29") # remove leap years
for(i in 1:eval(parse(text=paste("length(files",years[m],")",sep='')))){
  data<-read.csv(eval(parse(text=paste("files",years[m],"[",i,"]",sep=""))),header=FALSE)
  colnames(data)<-colnames
  data[,5]<-data[,5]-data[1,5]+1 # get correct day of year
  for(j in 1:365){
  daily_sub<-subset(data,doy==j)
  x<-cbind(daily_sub[,2],daily_sub[,3])
  grid <- rasterize(x, gridout, daily_sub$meanmid)
  #plot(grid,main=j)
  if(j==1){s<-grid}else{s<-stack(s,grid)}
  cat(j,'of 365 ', years[m],'sim ',sims[i], '\n')  
  }
  
  filename<-paste("meanmid_",years[m],"_",sims[i],".nc",sep="")
  writeRaster(s, filename=filename, overwrite=TRUE)
  temp<-open.ncdf( filename, write=TRUE, readunlim=TRUE)
  put.var.ncdf( temp, varid='z', vals=dates)
  close.ncdf(temp)  
  # yearly data
  } #end loop through sims
} #end loop through years

#arvo
for(m in 1:length(years)){
dates<-seq(ISOdate(years[m],1,1,tz=tzone)-3600*12, ISOdate(years[m]+1,1,1,tz=tzone)-3600*13, by="days")
dates<-subset(dates, format(dates, "%m/%d")!= "02/29") # remove leap years
for(i in 1:eval(parse(text=paste("length(files",years[m],")",sep='')))){
  data<-read.csv(eval(parse(text=paste("files",years[m],"[",i,"]",sep=""))),header=FALSE)
  colnames(data)<-colnames
  data[,5]<-data[,5]-data[1,5]+1 # get correct day of year
  for(j in 1:365){
  daily_sub<-subset(data,doy==j)
  x<-cbind(daily_sub[,2],daily_sub[,3])
  grid <- rasterize(x, gridout, daily_sub$arvo)
  #plot(grid,main=j)
  if(j==1){s<-grid}else{s<-stack(s,grid)}
  cat(j,'of 365 ', years[m],'sim ',sims[i], '\n')  
  }
  
  filename<-paste("arvo_",years[m],"_",sims[i],".nc",sep="")
  writeRaster(s, filename=filename, overwrite=TRUE)
  temp<-open.ncdf( filename, write=TRUE, readunlim=TRUE)
  put.var.ncdf( temp, varid='z', vals=dates)
  close.ncdf(temp)  
  # yearly data
  } #end loop through sims
} #end loop through years

var72_85_1997<-brick('maxrun_1997_var72_85.nc')
con72_1997<-brick('maxrun_1997_con72.nc')
con78_1997<-brick('maxrun_1997_con78.nc')
con85_1997<-brick('maxrun_1997_con85.nc')

diff_var_con72_1997<-(var72_85_1997-con72_1997)
plot(sum(diff_var_con72_1997)/(sum(con72_1997)+sum(var72_85_1997))*100)
diff_var_con78_1997<-var72_85_1997-con78_1997
plot(sum(diff_var_con78_1997)/(sum(con78_1997)+sum(var72_85_1997))*100)
diff_var_con85_1997<-var72_85_1997-con85_1997
plot(sum(diff_var_con85_1997)/(sum(con85_1997)+sum(var72_85_1997))*100)

diff_con72_con78_1997<-(con72_1997-con78_1997)
plot(sum(diff_con72_con78_1997)/sum(con72_1997)*100)
diff_con72_con85_1997<-(con72_1997-con85_1997)
plot(sum(diff_con72_con85_1997)/sum(con72_1997)*100)
diff_con78_con85_1997<-(con78_1997-con85_1997)
plot(sum(diff_con78_con85_1997)/sum(con78_1997)*100)



var72_85_1997<-brick('maxrun_1997_var72_85.nc')
var72_85_1997<-brick('sumact_1997_var72_85.nc')
var72_85_1997<-brick('mornbask_1997_var72_85.nc')
plot(var72_85_1997[[100:108]])

