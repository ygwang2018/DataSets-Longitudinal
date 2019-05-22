dat=read.table("C:\\Users\\fuly737\\Dropbox\\GEEBook2016\\dataset\\hiv315.txt",header=TRUE)
head(dat)
dat=data.frame(dat)
n=length(unique(dat[,1]))
plot(dat$lgcopy~dat$day,type="n",xlab="Days",ylab=expression(log10("RNA")))
for(i in 1:n){
  lines(dat$day[dat$id==i],dat$lgcopy[dat$id==i],col=i)
}

