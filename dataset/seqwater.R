
library(MASS)
#library(quantreg)
library(gdata)
library(lattice)

library(chron)


pfile<- 'C:/perl/bin/perl.exe'
allsite1<- read.delim("C:\\Users\\fuly737\\Dropbox\\GEEBook2016\\dataset\\wiv_QW.txt",header=TRUE,sep="\t")

## allsite1<- read.table("Y:/SEQWater/wiv_QW.txt",h=T)
dim(allsite1)
head(allsite1)

wih<-as.data.frame(allsite1)
unique(wih$Name)
length(unique(wih$Site))


datetxt <- as.Date(wih$Date,"%d/%m/%Y")
df <- data.frame(
  date = datetxt,
  year = as.numeric(format(datetxt, format = "%Y")),
  month = as.numeric(format(datetxt, format = "%m")),
  day = as.numeric(format(datetxt, format = "%d")))
wih$Year=df$year
wih$Month=df$month
wih$Date=as.Date(wih$Date,"%d/%m/%Y")
unique(wih$Year)
subset1=wih[wih$Name=="Total (CYANOPHYTES)",]
length(unique(subset1$Site))
table(subset1$Site)


xyplot(Val~Date,data=subset1,groups=Site)
subset1$Site=factor(subset1$Site)
xyplot(Val~Date|Site,data=subset1)
dat4=subset1[subset1$Site=="30004",]
dat16=subset1[subset1$Site=="30016",]
dat17=subset1[subset1$Site=="30017",]
dat18=subset1[subset1$Site=="30018",]
dat=rbind(dat4,dat16,dat17,dat18)
xyplot(Val~Date|Site, data=dat,type="b",xlab="Year",ylab="Values of Total Cyanophytes mg/L")
qqnorm(subset1$Val)
qqline(subset1$Val)

#########################################3333333




subset2=wih[wih$Name=="Chlorophyll a",]
length(unique(subset2$Site))
table(subset2$Site)


library(lattice)
subset2$Site=factor(subset2$Site)

xyplot(Val~Year|Site, data=subset2,type="b",xlab="Year",ylab="Values of Total Cyanophytes mg/L")

subset3=wih[wih$Name=="Total Nitrogen",]
length(unique(subset3$Site))
table(subset3$Site)
subset3$Site=factor(subset3$Site)

xyplot(Val~Year|Site, data=subset3,type="b",xlab="Year",ylab="Values of Total Cyanophytes mg/L")
