###################chapter 2

###################################
#Example 1 HIV 
##################################
rm(list=ls())
setwd("E:\\2016GEEBookEnglish\\Rcode")
library(lattice)
require(nlme)

######################### Beginning of Splus Code ########################
# Input Data: assume that the data at the end is stored in an ASCII file
# named "data".
workd<-read.table(file="ACTG315-Ding.dat",header=T,row.names=NULL)
# The measurements before treatment initiation are not used in the analysis.

#workd<-workd[workd$Day>0,]
# Impute below detectable viral measurement by 50.
workd$RNA[workd$RNA==100]<-50

workd<- groupedData(RNA~Day|ID, data=workd)

ID<-unique(workd$ID)

n<-length(ID)

with(workd,plot(Day,log10(RNA),type='n',xlab="Days"))
for (i in 1:n)
{xx<-(workd$ID==ID[i])
lines(workd$Day[xx],log10(workd$RNA[xx]),col=i)
}




################################################################
#Example 2 Progabide study 
################################################################

seiz<-as.data.frame(read.table(file="seizure.data.txt",header=FALSE, sep=""))

names(seiz)=c("id", "counts","visit","trt","age","weeks")
seiz$medicine=as.factor(rep(c('placebo','progabide'), as.numeric(table(seiz$trt))))
head(seiz)

#####################################fig boxplot
library(lattice)
seiz$visit=factor(seiz$visit)
bwplot(counts~visit|medicine,  layout=c(2,1), xlab="Visit",
ylab="The number of seizures",strip = strip.custom(strip.names = TRUE, strip.levels = TRUE),horizontal=FALSE,data=seiz, outline=TRUE, col=1)

savePlot(filename = "ex2-seizure-boxplot", type = "eps",device = dev.cur(),restoreConsole = TRUE)
savePlot(filename = "ex2-seizure-boxplot", type = "pdf",device = dev.cur(),restoreConsole = TRUE)

savePlot(filename = "ex2-seizure-boxplot", type = "jpg",device = dev.cur(),restoreConsole = TRUE)


############################blue for placebo; red for progabide
library(geepack)
data(seizure)
colnames(seizure)[c(6,1:4)]=c("baseline","visit 1","visit 2","visit 3","visit 4")
seizure$trt=as.factor(seizure$trt)
pairs(seizure[c(6,1:4)],pch=21,bg = c("blue","red")[unclass(seizure$trt)],main="Progabide study--two treatments")
savePlot(filename = "ex2-seizure-corplot", type = "eps",device = dev.cur(),restoreConsole = TRUE)
savePlot(filename = "ex2-seizure-corplot", type = "pdf",device = dev.cur(),restoreConsole = TRUE)


###############################################################
#Example 3 horme data
###############################################################
horme=as.data.frame(read.table(file="control.dat",header=FALSE))
colnames(horme)=c("id", "cycle","age","bmi","day", "logpdg" )

head(horme)
attach(horme)
fm=loess(logpdg~day,horme)

plot(logpdg~day,ylab="Log progesterone",xlab="Days in standardized menstrual cycle", xaxt = "n")
axis(side=1, at =seq(0,28,by=2), labels = TRUE, tick = TRUE,lty = "solid", col = NULL, col.ticks = NULL)
lines((spline(day,fitted(fm))),col="blue")
savePlot(filename = "ex3-hormeplot", type = "eps",device = dev.cur(),restoreConsole = TRUE)
savePlot(filename = "ex3-hormeplot", type = "pdf",device = dev.cur(),restoreConsole = TRUE)

detach(horme)


###############################################################
#Example 4 Teratology Study  data??
###############################################################
library(nlme)
workd<-read.table(file="ACTG315-Ding.dat",header=T,row.names=NULL)


#workd$RNA[workd$RNA==100]<-50

workd<- groupedData(RNA~Day|ID, data=workd)

ID<-unique(workd$ID)

n<-length(ID)

with(workd,plot(Day,log10(RNA),type='n',xlab="Days"))
for (i in 1:n)
{xx<-(workd$ID==ID[i])
lines(workd$Day[xx],log10(workd$RNA[xx]),col=i)
}

savePlot(filename = "ex4-teratologyplot", type = "eps")
savePlot(filename = "ex4-teratologyplot", type = "pdf")

###############################################################
#Example 5 Schizophrenia Study
###############################################################
mad<- read.table( file="madras.data.txt",sep="", header=TRUE)

# Name the variables according to madras.txt
colnames(mad) <- c(
  "ID",              # patient ID
  "Y",               # symptom indicator: binary
  "MONTH",           # month since hospitalization
  "AGE",             # age-at-onset (1= age<20 ; 0= age>=20)
  "GENDER",          # gender (1=female; 0=male)
  "MONTH.AGE",       # interaction term MONTH x AGE
  "MONTH.GENDER"     # interaction term MONTH x GENDER
)

head(mad,40)


with(mad, ftable(Y,AGE,GENDER,MONTH))


######################################################################
#Example 6 labor pain £¨Jung 1996£©
###################################################################

laborp=read.table(file="laborpain.txt",header=T,sep=" ")
head(laborp,10)
labor=reshape(laborp,varying=list(names(laborp)[-c(1,8)]), direction='long')

names(labor)=c("patient","trt","time","y","id")
labor$time=labor$time*30
labor=labor[!is.na(labor$y),]

labor=labor[order(labor$id),]
head(labor)
labor$group=rep(c("Active","Placebo"),times=c(189,169))
###############################boxplot
bwplot(y~time|group,  data=labor,pars = list(boxwex = 0.15, staplewex = 0.25, outwex = 0.5),
       horizontal =F, groups=as.factor(trt),xlab="Measurement times",ylab="Pain score",col=1,layout=c(1,2))

savePlot(filename = "ex6-laborpain-boxplot", type = "eps",device = dev.cur(),restoreConsole = TRUE)
savePlot(filename = "ex6-laborpain-boxplot", type = "pdf",device = dev.cur(),restoreConsole = TRUE)
###############################parallel line plot
xyplot(y~time|group, groups=id,data=labor,type="l",xlab="Measurement time (minutes)",ylab="Pain score")

savePlot(filename = "ex6-laborpain-parallelplot", type = "eps",device = dev.cur(),restoreConsole = TRUE)
savePlot(filename = "ex6-laborpain-Parallelplot", type = "pdf",device = dev.cur(),restoreConsole = TRUE)

#############################histogram plot
labora=labor[labor$trt==1,]
hist(labora$y,freq=F,main="Histogram for active group ",xlab="Pain score ")
savePlot(filename = "ex6-active-histplot", type = "eps",device = dev.cur(),restoreConsole = TRUE)
savePlot(filename = "ex6-active-histplot", type = "pdf",device = dev.cur(),restoreConsole = TRUE)


laborp=labor[labor$trt==0,]
hist(laborp$y,freq=F,main="Histogram for placebo group ",xlab="Pain score ")
savePlot(filename = "ex6-placebo-histplot", type = "eps",device = dev.cur(),restoreConsole = TRUE)
savePlot(filename = "ex6-placebo-histplot", type = "pdf",device = dev.cur(),restoreConsole = TRUE)


############################################
#Example 7 Labor Market Experience
###############################################

wage=as.data.frame(read.table(file="nlslme.csv",header = TRUE, sep=","))
wage=wage[!is.na(wage$age),]

length(unique(wage$idcode))

wage$south=as.factor(wage$south)
levels(wage$south)=c("Others","South")
#########################
plot(lnwage~grade,xlab="Education (year)",ylab="Log wage", data=wage, type="p", col=1 )


savePlot(filename = "ex7-wageEDplot", type = "eps",device = dev.cur(),restoreConsole = TRUE)
savePlot(filename = "ex7-wageEDplot", type = "pdf",device = dev.cur(),restoreConsole = TRUE)

bwplot(lnwage~grade, horizontal =F,xlab="Education (year)",ylab="Log wage", data=wage, type="p", col=1 )
savePlot(filename = "ex7-wageED-boxplot", type = "eps",device = dev.cur(),restoreConsole = TRUE)
savePlot(filename = "ex7-wageED-boxplot", type = "pdf",device = dev.cur(),restoreConsole = TRUE)



xyplot(lnwage~grade|south,data=wage, horizontal =F, xlab="Education (year)",ylab="Log wage",type="p",col=1)
     
savePlot(filename = "ex7-southplot", type = "eps",device = dev.cur(),restoreConsole = TRUE)
savePlot(filename = "ex7-southplot", type = "pdf",device = dev.cur(),restoreConsole = TRUE)


xyplot(lnwage~age|south,data=wage, horizontal =F, xlab="Age",ylab="Log wage",type="p",col=1)
  
savePlot(filename = "ex7-wageAgeplot", type = "eps",device = dev.cur(),restoreConsole = TRUE)
savePlot(filename = "ex7-wageAgeplot", type = "pdf",device = dev.cur(),restoreConsole = TRUE)


###############################################
#Example 8 water quality data
##################################################
###############################################################################
