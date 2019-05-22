###################chapter 5 model selection#######################
rm(list=ls())
##########################functions########################
setwd("E:\\2016GEEBookEnglish\\Rcode")

source("CriteriaFunction.R")

#################################main codes################################


library(geepack)
library(emplik)
library(MESS)
###################################
#Model selection Example 1: hormon data 
##################################


horm<-as.data.frame(read.table(file="control.dat",header=FALSE, sep=""))
colnames(horm)=c("id", "cycle","age","bmi","day", "logpdg" )
m=length(unique(horm$id))
ni=as.numeric(table(horm$id))
horm$id=rep(1:m,ni)
head(horm)
y=horm$logpdg


horm$TIME=(horm$day-mean(horm$day))/10
horm$AGE=(horm$age-mean(horm$age))/100
horm$BMI=(horm$bmi-mean(horm$bmi))/100

############################correlation structure is independence
model1.in=geeglm(logpdg~AGE+BMI+TIME,data=horm, id=id,family=gaussian, corstr="independence")
model2.in=update(model1.in,.~.+I(TIME^2))
model3.in=update(model2.in,.~.+I(TIME^3))
model4.in=update(model3.in,.~.+I(TIME^4))
model5.in=update(model4.in,.~.+I(TIME^5))
#############################################################calculate criteria   QIC and EQIC
model1=QIC.CIC.geeglm(datatype="continuous",model.geeglm=model1.in, model.independence=model1.in)
model2=QIC.CIC.geeglm(datatype="continuous",model.geeglm=model2.in, model.independence=model2.in)
model3=QIC.CIC.geeglm(datatype="continuous",model.geeglm=model3.in, model.independence=model3.in)
model4=QIC.CIC.geeglm(datatype="continuous",model.geeglm=model4.in, model.independence=model4.in)
model5=QIC.CIC.geeglm(datatype="continuous",model.geeglm=model5.in, model.independence=model5.in)
resin=rbind(model1,model2,model3,model4,model5)


######################################Gaussian criteria
p=3
pp=c(p,p+1,p+2,p+3,p+4)
gausin=sapply(list(model1.in,model2.in,model3.in,model4.in,model5.in),function(x) Gau.geeglm("continuous",x))
GAIC.in=gausin+2*pp
GBIC.in=gausin+log(m)*pp

hormon.gausind=cbind(GAIC.in,GBIC.in)
hormon.ind=cbind(resin,hormon.gausind)
print(hormon.ind, digits=6)
##################################################################################################
library(MESS)
rbind(QIC(model1.in),QIC(model2.in),QIC(model3.in),QIC(model4.in),QIC(model5.in))

############################correlation structure is exchangeable
model1.ex=update(model1.in, corstr="exchangeable")
model2.ex=update(model1.ex,.~.+I(TIME^2))
model3.ex=update(model2.ex,.~.+I(TIME^3))
model4.ex=update(model3.ex,.~.+I(TIME^4))
model5.ex=update(model4.ex,.~.+I(TIME^5))
#############################################################calculate criteria
model1ex=QIC.CIC.geeglm(datatype="continuous",model.geeglm=model1.ex, model.independence=model1.in)
model2ex=QIC.CIC.geeglm(datatype="continuous",model.geeglm=model2.ex, model.independence=model2.in)
model3ex=QIC.CIC.geeglm(datatype="continuous",model.geeglm=model3.ex, model.independence=model3.in)
model4ex=QIC.CIC.geeglm(datatype="continuous",model.geeglm=model4.ex, model.independence=model4.in)
model5ex=QIC.CIC.geeglm(datatype="continuous",model.geeglm=model5.ex, model.independence=model5.in)

resex=rbind(model1ex,model2ex,model3ex,model4ex,model5ex)
print(resex,digits=4)
rbind(QIC(model1.ex),QIC(model2.ex),QIC(model3.ex),QIC(model4.ex),QIC(model5.ex))

################################Gaussian criteria

gausex=sapply(list(model1.ex,model2.ex,model3.ex,model4.ex,model5.ex),function(x) Gau.geeglm("continuous",x))
GAIC.ex=gausex+2*pp
GBIC.ex=gausex+log(m)*pp

hormon.gauex=cbind(GAIC.ex,GBIC.ex)
hormon.exc=cbind(resex,hormon.gauex)

print(hormon.exc,digist=5)
rbind(QIC(model1.ex),QIC(model2.ex),QIC(model3.ex),QIC(model4.ex),QIC(model5.ex))


##########################################################################
# model selection: Madras Longitudinal Schizophrenia Study  correlation structure selection
############################################################################



madras<- read.table( file="madras.data.txt",sep="", header=TRUE)

# Name the variables according to madras.txt
colnames(madras) <- c(
  "ID",              # patient ID
  "Y",               # symptom indicator: binary
  "MONTH",           # month since hospitalization
  "AGE",             # age-at-onset (1= age<20 ; 0= age>=20)
  "GENDER",          # gender (1=female; 0=male)
  "MONTH.AGE",       # interaction term MONTH x AGE
  "MONTH.GENDER"     # interaction term MONTH x GENDER
)

head(madras,40)
ni=as.numeric(table(madras$ID))
m=length(unique(madras$ID)) #sample size
madras$id=rep(1:m,ni)
###############################geepack

madras.geeglmind=geeglm(Y~MONTH+AGE+GENDER+MONTH.AGE+MONTH.GENDER,data=madras,id=id,family=binomial,corstr="independence")
madras.geeglmexc=update(madras.geeglmind,corstr="exchangeable")
madras.geeglmar1=update(madras.geeglmind,corstr="ar1")

 summary(madras.geeglmind)
 summary(madras.geeglmexc)
 summary(madras.geeglmar1)
###################################QIC and CIC
madrasQICCIC=sapply(list(madras.geeglmind, madras.geeglmexc, madras.geeglmar1), function(x) QIC.CIC.geeglm(datatype="binary",x,madras.geeglmind))
colnames(madrasQICCIC)=c("IN","EX","AR")
print(madrasQICCIC,digits=4)

print(QIC(madras.geeglmind),digits=4)
print(QIC(madras.geeglmexc),digits=4)
print(QIC(madras.geeglmar1),digits=4)
##################################################GAIC and GBIC
madras.gauind=Gau.geeglm(datatype="binary",model.geeglm=madras.geeglmind)

madras.gauexc=Gau.geeglm(datatype="binary",model.geeglm=madras.geeglmexc)
madras.gauar1=Gau.geeglm(datatype="binary",model.geeglm=madras.geeglmar1)

madras.gau=c(madras.gauind,madras.gauexc,madras.gauar1)
p=6
np=c(p,p+1,p+1)
madrasGAIC=madras.gau+2*np
madrasGBIC=madras.gau+log(m)*np

madras.gau=rbind(madrasGAIC,madrasGBIC)
print(madras.gau,digits=5)
#####################################EAIC EBIC
ELind=elr.stall(datatype="binary",madras.geeglmind)
ELexc=elr.stall(datatype="binary",madras.geeglmexc)
ELar1=elr.stall(datatype="binary",madras.geeglmar1)

EL=c(ELind,ELexc,ELar1)
madras.EAIC=EL+2*np
madras.EBIC=EL+np*log(m)


madres=rbind(madrasQICCIC,madras.gau,madras.EAIC,madras.EBIC)
print(madres,digits=5)
#######################################################################################################
########################################example 3######################################################

data(seizure)
## Diggle, Liang, and Zeger (1994) pp166-168, compare Table 8.10
seizure<- reshape(seizure,
              varying=list(c("y1", "y2", "y3", "y4")),
              v.names="y", times=1:4, direction="long")


loc=with(seizure,id[which.max(y)])
#seizure=seizure[seizure$id!=loc,]  #remove the outliers
m=length(unique(seizure$id))

seizure<- seizure[order(seizure$id, seizure$time),]
seizure$visit=rep(c(0,0,0,1),m)
M=length(seizure[,1])

ni=as.numeric(table(seizure$id))

YY=with(seizure,matrix(y,nc=4,byrow=T))
cor(YY)

x=cbind(rep(1,M),log(seizure$base/4),seizure$trt,log(seizure$age),seizure$visit,log(seizure$base/4)*seizure$trt)
#############################

seiz.geeglmind=geeglm(y ~ log(base/4)+trt+log(age)+visit+log(base/4):trt, id = id,data=seizure, corstr="independence", family=poisson)
seiz.geeglmexc=update(seiz.geeglmind,corstr="exchangeable")
seiz.geeglmar1=update(seiz.geeglmind, corstr="ar1")
summary(seiz.geeglmind$geese)
summary(seiz.geeglmexc$geese)
summary(seiz.geeglmar1$geese)
#####################################################



######################################QIC and CIC
seizQIC=sapply(list(seiz.geeglmind, seiz.geeglmexc, seiz.geeglmar1), function(x) QIC.CIC.geeglm(datatype="count",x,seiz.geeglmind))
colnames(madrasQICCIC)=c("IN","EX","AR")

rbind(QIC(seiz.geeglmind),QIC(seiz.geeglmexc),QIC(seiz.geeglmar1))
##################################################GAIC and GBIC
seiz.gauind=Gau.geeglm(datatype="count",model.geeglm=seiz.geeglmind)

seiz.gauexc=Gau.geeglm(datatype="count",model.geeglm=seiz.geeglmexc)
seiz.gauar1=Gau.geeglm(datatype="count",model.geeglm=seiz.geeglmar1)

seiz.gau=c(seiz.gauind,seiz.gauexc,seiz.gauar1)
p=dim(x)[2]
np=c(p,p+1,p+1)
seizGAIC=seiz.gau+2*np
seizGBIC=seiz.gau+log(m)*np
#####################################EAIC EBIC

seiz.ELind=elr.stall("count",seiz.geeglmind)
seiz.ELexc=elr.stall("count",seiz.geeglmexc)
seiz.ELar1=elr.stall("count",seiz.geeglmar1)

seiz.ELC=c(seiz.ELind,seiz.ELexc,seiz.ELar1)


seiz.EAIC=seiz.ELC+2*np
seiz.EBIC=seiz.ELC+np*log(m)
##################################################CR
seiz.CRind=shof("count",seiz.geeglmind)
seiz.CRexc=shof("count",seiz.geeglmexc)
seiz.CRar1=shof("count",seiz.geeglmar)

seiz.res=rbind(seizQIC,seizGAIC,seizGBIC,seiz.EAIC,seiz.EBIC)
########################################################
print(seiz.res,digits=5)
