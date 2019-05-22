
library(MASS)
#library(quantreg)
library(gdata)
library(lattice)

library(chron)

xlsfile <- file.path("Y:/SEQWater/wiv_QW.xls")
datadir <- file.path("Y:/SEQWater")
 pfile<- 'C:/perl/bin/perl.exe'
 allsite1<- read.xls(xlsfile,sheet=1,perl=pfile)

## allsite1<- read.table("Y:/SEQWater/wiv_QW.txt",h=T)


wih<-as.data.frame(allsite1)

names(wih) <- tolower(names(wih))
wih$site<-paste(wih$site)

wih$time=paste(wih$time)
wih$date=paste(wih$date)
d = strptime(paste(wih$date,wih$time),format="%Y-%m-%d %H:%M:%S")
nobs=dim(wih)[1]
 newdate<- as.data.frame(matrix(unlist(as.POSIXlt(d)),nr=nobs))
 names(newdate)<- names(unlist(as.POSIXlt(wih$date[1])))
 wih<-cbind(wih,newdate) 
T<- strptime(format(wih$date,format="%Y-%m-%d"),format="%Y-%m-%d")

# in seconds;
wih$days<- T - strptime("30/07/1997", format="%d/%m/%Y") 
wih$days<-as.numeric(wih$days)/(24*3600) +1
wih$times<- as.numeric(wih$hour + wih$min/60 + wih$sec/3600) 
wih$year<-as.numeric(format(T,format="%Y"))
wih$mon<-as.numeric(format(T,format="%m")) 
wih$day<-as.numeric(format(T,format="%d")) 
parname<-unique(wih$name)

wih$cens<- (wih$qualifier=="<")
 
 xx<- tapply(wih$y,wih$name,median)
 
summary(wih)

tapply(wih$cens,wih$name,mean)

dn=as.numeric(row.names(xx))
xx=as.data.frame(xx)
xx$depth=dn
xx= (wih$name %in% parname[1:6])
wih1=wih[xx,]
wih2=wih[!xx,]

wih=wih[order(wih$site),]


par(mfrow=c(3,2))
plot()



for (i in 1:2)
 {
#cat("Loop ", i, "\n")
 
i=8

xx=(wih$name==parname[i])

 sx= ( & (!(is.na(wih$y))))
wih0=wih[sx,]
tt <- parname[i]

xyplot(y ~ dateandtime | name, data=wih2,span=1,aspect=1, 
xlab=paste(tt),
scales = list(relation = "free"),panel=function(x,y,span){
panel.xyplot(x,y); panel.loess(x,y,span)} )

}


par(mfrow=c(3,4))
for (i in 1:12)
 {
 xx=(wih$name==parname[i])
 sx= (xx & (!(is.na(wih$y))))

 plot(wih$dateandtime[sx],wih$y[sx],ylab=paste(parname[i]), xlab="Date")
}


 x1=as.data.frame(t(matrix(unlist(xx),nc=12)))
 names(x1)=names(xx)
  x1$par<-names(xx)

 x1$pcens=tapply(wih$cens,wih$name,FUN=mean)
x1<-x1[,c(7,8,1:6)]



for (i in 1:length(parname))
 {
  data1<-wih[wih$name==parname[i],]
  print(parname[i])
}


plot(data1$days, log(data1$y))



wih$amm<-paste(wih$amm)
wih$ammlim<-(substr(wih$amm,1,2)== "<3")
wih$amm[wih$ammlim]<-c("3")
wih$amm<-as.numeric(paste(wih$amm))

wih$orth<-paste(wih$orth)
wih$orthlim<-(substr(wih$orth,1,2) =="<2")
wih$orth[wih$orthlim]<-c("2")
wih$orth<-as.numeric(paste(wih$orth))

wih$no23<-paste(wih$no23)
wih$no23lim<-(substr(wih$no23,1,2)=="<2")
wih$no23[wih$no23lim]<-c("2")
wih$no23<-as.numeric(paste(wih$no23))

wih$chl<-paste(wih$chl)
wih$chllim<-(substr(wih$chl,1,4)=="<0.1")
wih$chl[wih$chllim]<-c("0.1")
wih$chl<-as.numeric(paste(wih$chl))

wih$phae<-paste(wih$phae)
wih$phaelim<-(substr(wih$phae,1,4)=="<0.1")
wih$phae[wih$phaelim]<- "0.1"
wih$phae<- as.numeric(paste(wih$phae))
wih$loc <- substr(paste(wih$site),1,2)

wih<-as.data.frame(wih)

splom(~wih[,4:7] | wih$loc, pscales=0)

attach(wih)

plot(days,wih[,4],fin=c(2.3,2),type="n")


attach(wih)

par(mfrow=c(4,3),omi=c(0.3,0.3,0.3,0.3),mai=0.1*rep(0,4))
LOC<- unique(wih$loc)

yax<- c("s","n","n")
for (i in 1:4)
{
ymax<-max(wih[!is.na(wih[,i+2]),i+2])
  for (j in 1:3)
 
{ xx<- (loc==LOC[j])& (!is.na(wih[,i+2]))
  plot(T[xx],wih[xx,i+2],xlab="Days",ylab=names(wih)[i+2],col=2+(substr(wih$site[xx],3,3)=="R"),
  bty="c", ylim=c(0,ymax),yaxt=yax[j])
  text(T[100],0.8*ymax, LOC[j])
   text(T[100],0.65*ymax, names(wih)[i+2])
 }
}

detach(wih)

##title.string<-"Figure 1: Nutrient concentrations at three sites (corresponding to 3 columns)"
## mtext(title.string, side=3, outer=T) 
