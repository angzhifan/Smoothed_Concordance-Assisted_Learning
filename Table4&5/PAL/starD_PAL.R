#STAR*D PAL 
#Date: 12/21/20
data2=read.csv(file = "stard.csv",header = T)


#PAL
treat = rep(0,length(data2$tx))
treat[which(data2$tx=="BUP")]=1
covariates = as.matrix(data2[,4:308])+matrix(data = 0.0, 319,305) 
library(ITRSelect)
fit <- PAL.fit(-data2$qctot, x1=covariates,a1=treat)
summary(fit)
fit$pi1.est
fit$beta1.est

#constants
n=dim(data2)[1]
d=dim(data2)[2]-3
A=treat
pi=sum(A)/n
Y=-data2[,3]


#estimate threshold
Xbeta=covariates%*%as.matrix(fit$beta1.est[2:306])+fit$beta1.est[1]
sample1=which(Xbeta>0)
sample0=which(Xbeta<=0)
length(intersect(sample1,which(A==1)))
length(intersect(sample1,which(A==0)))
length(intersect(sample0,which(A==1)))
length(intersect(sample0,which(A==0)))


#Estimated Value
bootmean=rep(0,1000)
bootn=1000
for (i in 1:1000) {
  boot=ceiling(runif(bootn,0,n))
  IPW=Y[boot]*is.element(boot,union(intersect(sample1,which(A==1)),intersect(sample0,which(A==0))))/(A[boot]*pi+(1-A[boot])*(1-pi))
  bootmean[i]=sum(IPW)/bootn
}
mean(bootmean)
sort(bootmean)[25]
sort(bootmean)[975]
