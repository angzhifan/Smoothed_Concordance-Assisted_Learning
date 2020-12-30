#Date: 12/25/20
library(ITRSelect)
library(MASS)

start=Sys.time()

#numsim is the times of simulations
numsim=100

#d is the dimension
d=50
#rho
rho=0.0
n=200


#sigma
if(rho==0){
  sig=diag(d)
}else{
  sig=matrix(data = rep(1:d,d),d,d)
  sig=abs(sig-t(sig))
  sig=rho^sig
}
testn=1000

#a matrix to store the simulation results
IPW=matrix(data = NA,100,1)
simresult=matrix(data=NA,nrow = numsim,ncol = 5) 
colnames(simresult)<-c("MSE","Incorr0(0)","Corr0(46)","Estimated Value","testPCD")

beta0=c(1,1,-1,1,0,0,rep(0,d-6))
beta0=beta0/sqrt((t(beta0)%*%beta0)[1,1])

#simulations
for (t in 1:numsim) {  
  
  #generate samples
  X=t(mvrnorm(n,rep(0,d),sig))
  A=rbinom(n,size = 1,prob = 0.5)
  Y=1+X[1,]-X[2,]+X[3,]+X[4,]+rnorm(n,0,0.5)+A*(exp(1+X[1,]+X[2,]-X[3,]+X[4,])-exp(1))
  v=mean(Y[which(A==0)])
  W=4*(Y-v)*(A-0.5)
  W1=matrix(data = rep(W,n),n,n)-t(matrix(data = rep(W,n),n,n))
  W1=W1-t(W1)
  
  
  #PAL
  fit <- PAL.fit(Y,x1=t(X),a1=A)
  
  #calculate MSE
  beta = fit$beta1.est[2:(d+1)]
  beta1=beta/sqrt((t(beta)%*%beta)[1,1])
  simresult[t,1]=t(beta1-beta0)%*%(beta1-beta0)
  

  #Incorr0(0)
  simresult[t,2]=length(intersect(which(beta==0),c(1,2,3,4)))
  
  #Corr0(48)
  simresult[t,3]=length(intersect(which(beta==0),c(5:d)))
  
  # threshold
  Xbeta=t(X)%*%beta
  b0=fit$beta1.est[1]
 
  #Estimated Value
  testX=t(mvrnorm(testn,rep(0,d),sig))
  testXbeta=t(testX)%*%beta
  testA=rep(0,testn) 
  testA[which(testXbeta>b0)]=1 
  testY=1+testX[1,]-testX[2,]+testX[3,]+testX[4,]+testA*(exp(1+testX[1,]+testX[2,]-testX[3,]+testX[4,])-exp(1))
  simresult[t,4]=mean(testY)
  
  #test PCD
  X2test=t(testX)%*%beta
  test1=which(X2test>b0)
  true1=which(testX[1,]+testX[2,]-testX[3,]+testX[4,]>0)
  simresult[t,5]=(length(intersect(test1,true1))+testn-length(union(test1,true1)))/testn
  
}

print(simresult)
print(colMeans(simresult))
print(sd(simresult[,4]))
print(sd(simresult[,5]))
Sys.time()-start

