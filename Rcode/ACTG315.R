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



# Define functions representing Model (12) and Model (13)
# of Wu and Ding (1999).
exp.model12<-function(p0, p1, d1, X)
{
P0 <- exp(p0)
P1 <- exp(p1)
P0 + P1 * exp( - d1 * X)
}

exp.model13<-function(p1, d1, p2, d2, X)
{
P1 <- exp(p1)
P2 <- exp(p2)
P1 * exp( - d1 * X) + P2 * exp( - d2 * X)
}

# Fit Model (12) and Model (13) by NLME
actg315.model12<-nlme(Y ~ log10(exp.model12(p0, p1, d1, X)),
fixed = list(p0~., p1~., d1~.),
random = list(p0~., p1~., d1~.),
cluster = ~ Z,
data = data.frame(Y=log10(workd$RNA), X=workd$Day, Z=workd$ID),
start = list(fixed= c(5,11,0.5)),verbose=T)

actg315.model13<-nlme(Y ~ log10(exp.model13(p1, d1, p2, d2, X)),
fixed = list(p1~., d1~., p2~., d2~.),
random = list(p1~., d1~., p2~., d2~.),
cluster = ~ Z,
data = data.frame(Y=log10(workd$RNA), X=workd$Day, Z=workd$ID),
start = list(fixed= c(12,0.4,7,0.03)),verbose=T)

anova(actg315.model12,actg315.model13)
###################### End of Splus Code ######################################

RESULTS:
--------

Model Df AIC BIC Loglik Test Lik.Ratio P value
actg315.model12 1 10 433.03 470.69 -206.52
actg315.model13 2 15 257.98 314.46 -113.99 1 vs. 2 185.05 0

Fixed Effects Estimates from Model (12):
LogP0 LogP1 delta_p
5.666933 11.10344 0.2246409

Fixed Effects Estimates from Model (13):
LogP1 delta_p LogP2 lamba_l
12.32537 0.4759455 7.945189 0.04134005

**********************************************************************************
8. Please cite the following references if appropriate when you use the data in your paper:

REFERENCES:
----------
1) Connick E, Lederman MM, Kotzin BL, et al. (2000), "Immune reconstitution in the first year of potent antiretroviral therapy and its relationship to virologic response," Journal of Infectious Diseases, 181:358-63.

2) Ding, A.A. and Wu, H. (1999), "Relationships between Antiviral Treatment Effects and Biphasic Viral Decay Rates in Modeling HIV Dynamics," Mathematical Biosciences, 160. 63-82.

3) Ding, A.A. and Wu, H. (2000), "A Comparison Study of Models and Fitting Procedures for Biphasic Viral Dynamics in HIV-1 Infected Patients Treated with Antiviral Therapies," Biometrics, 56, 293-300.

4) Lederman MM, Connick E, Landay A, et al. (1998), "Immunologic responses associated with 12 weeks of combination antiretroviral therapy consisting of zidovudine, lamivudine and ritonavir: results of AIDS Clinical Trials Group Protocol 315," Journal of Infectious Diseases, 178: 70-79.

5) Wu, H. and Ding, A. (1999), "Population HIV-1 Dynamics in Vivo: Applicable Models and Inferential Tools for Virological Data from AIDS Clinical Trials," Biometrics, 55, 410-418.

6) Wu, H. and Ding, A. and DeGruttola. V. (1998), "Estimation of HIV Dynamic Parameters," Statistics in Medicine, 17, 2463-2485.

7) Wu, H., Kuritzkes, D.R., and McClernon, D.R. et al. (1999), "Characterization of Viral Dynamics in Human Immunodeficiency Virus Type 1-Infected Patients Treated with Combination Antiretroviral Therapy: Relationships to Host Factors, Cellular Restoration and Virological Endpoints," Journal of Infectious Diseases, 179(4):799-807.
---------------------------------------------------------------------------

 
