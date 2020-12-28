#Low-dimensional Example
#Date: 12/22/20

library(ITRSelect)

start=Sys.time()

#numsim is the times of simulations
numsim=100

#d is the dimension
d=50

#the sample size is n, which takes the value of 30, 100 and 200.
n=200

#the simulation results
simresult=matrix(data=NA,nrow = numsim,ncol = 9)
colnames(simresult)<-c("MSE","Incorr0(0)","Corr0(48)","Within PCD","Estimated Value","lambda","alpha","beta[2]","testPCD")
testn=1000


#simulations
for (t in 1:numsim) {  
  
  #generate samples
  X=matrix(data=runif(d*n,-1,1),d,n)
  A<-rbinom(n,size = 1,prob = 0.5)
  Y=1+2*X[1,]+X[2,]+0.5*X[3,]+rnorm(n,0,1)+0.442*(rep(1,n)-X[1,]-X[2,])*(2*A-rep(1,n))
  v=mean(Y[which(A==0)])
  W=4*(Y-v)*(A-0.5)
  W1=matrix(data = rep(W,n),n,n)
  W1=W1-t(W1)
  
  
  #PAL
  fit <- PAL.fit(Y,x1=t(X),a1=A)
  
  #calculate MSE
  beta = fit$beta1.est[2:51]
  beta1=beta/sqrt((t(beta)%*%beta)[1,1])
  simresult[t,1]=t(beta1-rbind(-1/sqrt(2),-1/sqrt(2),matrix(rep(0,d-2))))%*%(beta1-rbind(-1/sqrt(2),-1/sqrt(2),matrix(rep(0,d-2))))
  
  #Incorr0(0)
  simresult[t,2]=length(intersect(which(beta==0),1:2))
  
  #Corr0(48)
  simresult[t,3]=length(intersect(which(beta==0),3:50))
  
  # threshold
  Xbeta=t(X)%*%beta
  b0=fit$beta1.est[1]
  
  #PCD within these n samples
  sample1=which((b0+Xbeta)>0)#
  truesample1=which(1-X[1,]-X[2,]>0)
  simresult[t,4]=(length(intersect(sample1,truesample1))+n-length(union(sample1,truesample1)))/n
  
  #Estimated value function
  testX=matrix(data=runif(d*testn,-1,1),d,testn)
  testXbeta=t(testX)%*%beta
  testA=rep(0,testn) 
  testA[which((b0+testXbeta)>0)]=1 
  testY=1+2*testX[1,]+testX[2,]+0.5*testX[3,]+rnorm(testn,0,1)+0.442*(1-testX[1,]-testX[2,])*(2*testA-1)
  simresult[t,5]=mean(testY)
  
  #test PCD
  X2test=t(testX)%*%beta
  test1=which((b0+X2test)>0)
  true1=which((rep(1,testn)-testX[1,]-testX[2,])>0)
  simresult[t,9]=(length(intersect(test1,true1))+testn-length(union(test1,true1)))/testn
  
  
}

print(simresult)
print(colMeans(simresult))
print(sd(simresult[,4]))
print(sd(simresult[,5]))
print(sd(simresult[,9]))
Sys.time()-start

