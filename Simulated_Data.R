######################################################################################
#############################POLI SPATIAL#############################################
######################################################################################


##This are the library used:
#library
library(invgamma)
library(matlib)
library(mvnfast)
library(mvtnorm)
library(spBayes)
library(classInt)
library(RColorBrewer)
library(spBayes)
library(MBA)
library(geoR)
library(fields)
library(sp)
library(maptools)
library(rgdal)
library(classInt)
library(lattice)
library("scatterplot3d")
library(ggplot2)
library(BAS)


#At the beginning, we just made the computation for simulated data. 

#This are the function that we implemented

# EUCLIDEAN DISTANCE
Eucl <- function(vect1, vect2) sqrt(sum((vect1 - vect2)^2))
#DECAY FUNCTION
kernel = function(d,phi) exp(-d^2/phi)


#We set a seed since we are taking random variable from different distribution
#and we want that our result can be repeated

set.seed(1)

###############################################################################
############################Environment########################################
###############################################################################


n0=200 #-> suppose to be the number of our monitors
s <- replicate(2,rnorm(n0,0,1.5))
#since we are working with simulated data, we decide to create patological case
s[1,] = c(4,4)
s[2,] = c(-3,-3)
s[3,] = c(2,2)
s[4,] = c(-5,-5)
s[5,] = c(-3,4)


plot(s,type='p',col=2,pch=19,xlab ='Latitude' ,ylab='Longitude', main='Monitor Location')

s1 <- s

#Then we created a grid in which our monitors lie

# coord for prediction
grid = seq(-5,5,by=0.2)
G = s2 = expand.grid(grid,grid) #51*51
s2=as.matrix(s2)

dim(s1)
dim(s2)
n1=2601 #51*51
n=n0

###############################################################################
#######################Construction of cell####################################
###############################################################################


#X=cbind(rep(1,n),rnorm(n,2.5,0.0001))##this is in the case if we are in just one cell

#we divided our grid in cell

x_B=rep(0,n)
for (i in 1:n){
  #primo quadrante
  if (s[i,1]>=2.5 && s[i,2]>=2.5) {x_B[i]=20}
  if (s[i,1]>=2.5 && 0<=s[i,2] && s[i,2]<=2.5) {x_B[i]=15}
  if (0<=s[i,1] && s[i,1]<=2.5 && s[i,2]>=2.5) {x_B[i]=12}
  if (0<=s[i,1]&& s[i,1]<=2.5 && 0<=s[i,2] &&s[i,2]<=2.5) {x_B[i]=11}
  
  #secondo quadrante
  if (s[i,1]<=(-2.5) && s[i,2]>=2.5) {x_B[i]=13}
  if (s[i,1]<=(-2.5) && 0<=s[i,2] && s[i,2]<=2.5) {x_B[i]=7}
  if ((-2.5)<=s[i,1] && s[i,1]<=0 && s[i,2]>=2.5) {x_B[i]=3}
  if ((-2.5)<=s[i,1] && s[i,1]<=0 && 0<=s[i,2] && s[i,2]<=2.5) {x_B[i]=10}
  
  #terzo quadrante
  if (s[i,1]<=(-2.5) && s[i,2]<=(-2.5)) {x_B[i]=30}
  if (s[i,1]<=(-2.5) && (-2.5) <=s[i,2] && s[i,2]<=0) {x_B[i]=25}
  if ((-2.5)<=s[i,1] && s[i,1]<=0 && s[i,2]<=(-2.5)) {x_B[i]=22}
  if ((-2.5)<=s[i,1] && s[i,1]<=0 && (-2.5)<=s[i,2] && s[i,2]<=0) {x_B[i]=11.5}
  
  #quarto quadrante
  if (s[i,1]>=2.5 && s[i,2]<=(-2.5)) {x_B[i]=54}
  if (s[i,1]>=2.5 && (-2.5)<=s[i,2] && s[i,2]<=0) {x_B[i]=37}
  if (0<=s[i,1] && s[i,1]<=2.5 && s[i,2]<=(-2.5)) {x_B[i]=32.5}
  if (0<=s[i,1] && s[i,1]<=2.5 && (-2.5)<=s[i,2] && s[i,2]<=0) {x_B[i]=10.5}}

v=cbind(s[,1],s[,2],x_B)
X=cbind(rep(1,n),x_B)





#X=cbind(rep(1,n),apply(s^2/10,1,sum))
#X=cbind(rep(1,n),2)# for now we keep the X fixed
X


###############################################################################
############################Prior Parameter####################################
###############################################################################



#PRIOR parameter for the coregionalization matrix


a00=exp(rnorm(1))
a10=rnorm(1)
a11=exp(rnorm(1))


#PRIOR for phi0-> discrete uniform


phi0_seq=as.vector(seq(0.0005,0.05,0.005))
phi0_seq
phi0_prior=rep(1/length(phi0_seq),length(phi0_seq))
phi0_prior
phi0=sample(phi0_seq,1,replace=FALSE,phi0_prior)


#PRIOR for phi1-> discrete uniform


phi1_seq=as.vector(seq(0.01,0.1,0.01))
phi1_seq
phi1_prior=rep(1/length(phi1_seq),length(phi1_seq))
phi1_prior
phi1=sample(phi1_seq,1,replace=FALSE,phi1_prior)

#Following the steps done in the Rasmussen book, we procede calculating the covariance
#matrix of the initial point
D0=matrix(NA,n,n)

for (i in 1:n) {
  for (j in 1:n) {
    D0[i,j] =Eucl(s1[i,], s1[j,])
  }
}
D0
COV00=kernel(D0,phi0)

#generate a list of COV00 for all the possible value of phi in the sequence
#this will be useful in the Gibbs sampler in order to take at every sample the 
#corrisponding phi

COV00_v=list()
for (i in 1:length(phi0_seq)){
COV00_v[[i]] = kernel(D0,phi0_seq[i])}

#COV00_2 = kernel(D0,phi0_seq[2])
#COV00_3 = kernel(D0,phi0_seq[3])
#COV00_4 = kernel(D0,phi0_seq[4])
#COV00_5 = kernel(D0,phi0_seq[5])
#COV00_6 = kernel(D0,phi0_seq[6])

#analogously, we did the same also for the new point
D1=matrix(NA,n,n)

for (i in 1:n) {
  for (j in 1:n) {
    D1[i,j] =Eucl(s1[i,], s1[j,])
  }
}
D1
COV11 = kernel(D1,phi1)

#generate a list of COV11 for all the possible value of phi in the sequence
COV11_v=list()
for (i in 1:length(phi1_seq)){
  COV11_v[[i]] = kernel(D1,phi1_seq[i])}

is_symmetric_matrix(COV11)


#PRIOR OF W0S


b0s=rep(0,n)
COV00
W0S= t(rmvnorm(1, rep(0,n), COV00+diag(n)*.00001))


#PRIOR OF W1S


b1s=rep(0,n)
COV11
W1S= t(rmvnorm(1, rep(0,n), COV11+diag(n)*.00001))


#PRIOR OF TAU2


a1=4
b1=0.2
tau2=rinvgamma(1,a1,b1)  


#PRIOR OF BETA => bivariate normal (mu, sigmaB)


b_mu_prior=c(0,0)
b_mu_prior
sigmaB= solve(0.5*diag(2))  # it is already the inverse
sigmaB

A=sigmaB*b_mu_prior

b=rmvnorm(1,b_mu_prior,0.5*diag(2)) ##they will be the first value of the Gibbs Sampler

b



###############################################################################
######################Value of our simulated y#################################
###############################################################################


y_mu=b[1]*rep(1,n)+b[2]*X[,2]+a00*W0S+a10*X[,2]*W0S+a11*X[,2]*W1S

y = rmvnorm(1, y_mu, tau2*diag(n))
y=t(y)


#z= b[1]*rep(1,n) + b[2]*X[,2]+ W0S+eps
#eps=rnorm(n,0,tau2)


#divido tra training and validation
nt=n-50
nv=50

sp=s1[(nt+1):n,]
yp=y[(nt+1):n,]

W0S=W0S[1:nt,]
W1S=W1S[1:nt,]
y=y[1:nt,]
X=X[1:nt,]
s1=s1[1:nt,]

n=nt

#Following the steps done in the Rasmussen book, we procede calculating the covariance
#matrix of the initial point
D0=matrix(NA,n,n)

for (i in 1:n) {
  for (j in 1:n) {
    D0[i,j] =Eucl(s1[i,], s1[j,])
  }
}
D0
COV00=kernel(D0,phi0)

#generate a list of COV00 for all the possible value of phi in the sequence
#this will be useful in the Gibbs sampler in order to take at every sample the 
#corrisponding phi

COV00_v=list()
for (i in 1:length(phi0_seq)){
  COV00_v[[i]] = kernel(D0,phi0_seq[i])}

#COV00_2 = kernel(D0,phi0_seq[2])
#COV00_3 = kernel(D0,phi0_seq[3])
#COV00_4 = kernel(D0,phi0_seq[4])
#COV00_5 = kernel(D0,phi0_seq[5])
#COV00_6 = kernel(D0,phi0_seq[6])

#analogously, we did the same also for the new point
D1=matrix(NA,n,n)

for (i in 1:n) {
  for (j in 1:n) {
    D1[i,j] =Eucl(s1[i,], s1[j,])
  }
}
D1
COV11 = kernel(D1,phi1)

#generate a list of COV11 for all the possible value of phi in the sequence
COV11_v=list()
for (i in 1:length(phi1_seq)){
  COV11_v[[i]] = kernel(D1,phi1_seq[i])}

is_symmetric_matrix(COV11)

b0s=rep(0,n)
b1s=rep(0,n)

###############################################################################
############################Gibbs Sampler######################################
###############################################################################


# Full conditional of beta parameter  =>  bivariate normal (muPost, sigmaBPost)
#muPost= solve(diag(1/tau2,dim(COV)[1])+solve(COV)) %*% ((( 1/tau2*(solve(t(X)%*%X)%*%t(X)%*%that)))
          #+ solve(COV)%*%b0s)
#muPost
#sigmaBpost=  solve(diag(1/tau2,length(COV))+COV)
##sigmaBpost

# Full conditional of tau^2 parameter   => inversegamma  (an, bn)


#an= a1+ n/2
#an
#mu=sum(y)/n
#bn= b1+ sum(y - mu)^2/2
#bn




S<-2000

#y=b0 +b1*x(b)+w(S)+ e(s)
#We initialize a matrix with 2*n+5 columns and S rows 

# it contains simulated values of the BIVARIATE Gibbs sampler MC 
THETA<-matrix(nrow=S,ncol=2*n+5)


# Initial point of the chain 
theta<-cbind(b[1] ,b[2],tau2,t(W0S),t(W1S),phi0,phi1) # Initial point of the chain 


THETA[1,]=theta



sost1=t(X)%*%X
C=diag(a00*rep(1,n)+a10*X[,2])
sost2=t(C)%*%C
C2=diag(a11*X[,2])
sost3=t(C2)%*%C2


################################
####### Gibbs sampling #########

##THETA (bo,b1,tau2,wos1,wos2,..wosN,w1s1,w1s2,...w1sN)

r=2
for(r in 2:2000) {
  
  # generate a new beta value from its full conditional
  if(r%%20==0){
    cat(paste(r,"\n"))
  }
  
  betahat= solve(sost1)%*%t(X)%*%(y-a00*THETA[r-1,4:(n+3)]-a10*X[,2]*THETA[r-1,4:(n+3)]-a11*X[,2]*THETA[r-1,(n+4):(2*n+3)]) #da mettere i valori dentro THETA
  muPost= solve((sost1)/THETA[r-1,3] +sigmaB) %*% ((( (sost1)/THETA[r-1,3])%*%betahat) + sigmaB%*%b_mu_prior)
  sigmaBpost=  solve((sost1)/THETA[r-1,3] +sigmaB)
  
  bpost= rmvnorm(1, muPost, sigmaBpost )
  theta[1]=bpost[1]
  theta[2]=bpost[2]
  
  # generate a new tau^2 from its full conditional
  mu=X%*%t(bpost)+a00*THETA[r-1,4:(n+3)]+a10*X[,2]*THETA[r-1,4:(n+3)]+a11*X[,2]*THETA[r-1,(n+4):(2*n+3)]
  an= a1+ n/2
  bn= b1+ sum((y - mu)^2)/2 
  
  theta[3]<- rinvgamma(1, an, bn)
  
  
  ## generate a new W0S from its full conditional
  
  
  that=solve(sost2)%*%t(C)%*%(y-X%*%t(bpost)-a11*X[,2]*THETA[r-1,(n+4):(2*n+3)]) 
  
  sigmaBpostW0=solve(((sost2)/ theta[3]+diag(n)*.00001)+solve(COV00+diag(n)*.00001))
  
  muPostW0= sigmaBpostW0%*% (((sost2)/ theta[3]+diag(n)*.00001)%*%that + solve(COV00+diag(n)*.00001)%*%b0s)
  
  
  bpostW0= rmvnorm(1, muPostW0, sigmaBpostW0+diag(n)*.00001)
  
  
  for (j in 1:n) {
    theta[j+3]=bpostW0[j]
  }
  
  ## generate a new W1S from its full conditional
  
  
  that2=solve(sost3)%*%t(C2)%*%(y-X%*%t(bpost)-C%*% theta[4:(n+3)])

  
  sigmaBpostW1=solve(((sost3)/ theta[3]+diag(n)*.00001)+solve(COV11+diag(n)*.00001))
  
  
  muPostW1= sigmaBpostW1 %*%  (((sost3)/ theta[3]+diag(n)*.00001)%*%that2 + solve(COV11+diag(n)*.00001)%*%b1s)
  
  
  bpostW1= rmvnorm(1, muPostW1, sigmaBpostW1+diag(n)*.00001)
  
  
  
  for (j in 1:n) {
    theta[j+n+3]=bpostW1[j]
  }
  
  ## generate a new phi0 from its full conditional
  
  phi0_post=rep(0,length(phi0_seq))
  for (i in 1:length(phi0_seq)) {
    Fernando=COV00_v[[i]]
    
    C=diag(a00*rep(1,n)+a10*X[,2])
    that=solve(sost2)%*%t(C)%*%(y-X%*%t(bpost)-a11*X[,2]*theta[(n+4):(2*n+3)]) 
    
    sigmaBpostW0=solve(((sost2)/ theta[3]+diag(n)*.00001)+solve(Fernando+diag(n)*.00001))
    
    
    muPostW0= sigmaBpostW0%*% (((sost2)/ theta[3]+diag(n)*.00001)%*%that + solve(Fernando+diag(n)*.00001)%*%b0s)
    
    bpostW0= rmvnorm(1, muPostW0, sigmaBpostW0+diag(n)*.00001)
    
    phi0_post[i]=log(1/length(phi0_seq))-log(sqrt(det(Fernando)))-((bpostW0)%*%solve(Fernando+diag(n)*.00001)%*%t(bpostW0))/(2)
    
  
  }
  
  phi0_post=exp(phi0_post)/sum(exp(phi0_post))
  theta[2*n+4]=sample(phi0_seq,1,replace=FALSE,phi0_post)
  
  
  COV00=COV00_v[[which(phi0_seq==theta[n*2+4])]]
  
  
    
   ## generate a new phi1 from its full conditional
  
  phi1_post=rep(0,length(phi1_seq))
  for (i in 1:length(phi1_seq)) {
    Fernando=COV11_v[[i]]
    
    that2=solve(sost3)%*%t(C2)%*%(y-X%*%t(bpost)-C%*% theta[4:(n+3)])
    
    
    sigmaBpostW1=solve(((sost3)/ theta[3]+diag(n)*.00001)+solve(Fernando+diag(n)*.00001))
    
    muPostW1= sigmaBpostW1 %*%  (((sost3)/ theta[3]+diag(n)*.00001)%*%that2 + solve(Fernando+diag(n)*.00001)%*%b1s)
    
    bpostW1= rmvnorm(1, muPostW1, sigmaBpostW1+diag(n)*.00001)
    
    phi1_post[i]=log(1/length(phi1_seq))-log(sqrt(det(Fernando+diag(n)*.00001)))-((bpostW1)%*%solve(Fernando+diag(n)*.00001)%*%t(bpostW1))/(2)
    
  }
  
  phi1_post=exp(phi1_post)/sum(exp(phi1_post))
  theta[2*n+5]=sample(phi1_seq,1,replace=FALSE,phi1_post)
  
  COV11=COV11_v[[which(phi1_seq==theta[n*2+5])]]
  
  
  theta
  THETA[r,]<-t(theta)
  THETA
}


plot(phi1_post)
plot(phi0_post)

THETA



##traceplots and histogram of the parameters of which we have assumed a prior

x11()
par(mfrow=c(2,3))
plot(ts(THETA[,1])) 
plot(ts(THETA[,2])) 
plot(ts(THETA[,3]))
hist(THETA[,2*n+4],breaks = 100)
hist(THETA[,2*n+5],breaks = 100)




###we decided to take the last 200 iterations and to take the average of them as final value

THETA_last=THETA[(S-1000):S,]

plot(ts(THETA_last[,1])) ###b0

plot(ts(THETA_last[,2])) ###b1
plot(ts(THETA_last[,3])) ###tau2

hist(THETA_last[,2*n+4],breaks = 100) ###phi0
hist(THETA_last[,2*n+5],breaks = 100) ###phi1


##this are our final result of the parameters
b0_full=mean(THETA_last[,1])
b0_full
b1_full=mean(THETA_last[,2])
b1_full
tau_full=mean(THETA_last[,3])
tau_full

w0s_monitor=colMeans(THETA_last[,4:(n+3)])

w1s_monitor=colMeans(THETA_last[,(n+4):(2*n+3)])



###############################################################################
############################Kriging############################################
###############################################################################

#compute the distances between points inside the grid 

##we compute the COV00 using the parameter with the highest frequency

D0=matrix(NA,n,n)

for (i in 1:n) {
  for (j in 1:n) {
    D0[i,j] =Eucl(s1[i,], s1[j,])
  }
}
D0
COV00=kernel(D0,mean(THETA_last[,2*n+4]))



#Then we created a grid in which our monitors lie

# coord for prediction
grid = seq(-5,5,by=0.2)
G = s2 = expand.grid(grid,grid) #51*51
s2=as.matrix(s2)
s2=rbind(s2,sp)


plot(s2)
points(sp,col=2)

dim(s1)
dim(s2)
n1=dim(s2)[1]

D01=matrix(NA,n1,n1)


for (i in 1:n1) {
  for (j in 1:n1) {
    D01[i,j] =Eucl(s2[i,], s2[j,])
  }
}
D01
COV01 = kernel(D01,mean(THETA_last[,2*n+4])) 
is_symmetric_matrix(COV01)

# compute the distances between new and old
D0mix=matrix(NA,n1,n)

for (i in 1:n1){
  for (j in 1:n) {
    D0mix[i,j]=Eucl(s2[i,], s1[j,])
    
  }
}
D0mix
COV0mix = kernel(D0mix,mean(THETA_last[,2*n+4]))


W0new=COV0mix%*%solve(COV00 + tau_full*diag(length(w0s_monitor)))%*%w0s_monitor


##I DO THE SAME WITH THE OTHER GP PROCESS


D1=matrix(NA,n,n)

for (i in 1:n) {
  for (j in 1:n) {
    D1[i,j] =Eucl(s1[i,], s1[j,])
  }
}
D1
COV11=kernel(D1,THETA[S,2*n+5])

D12=matrix(NA,n1,n1)


for (i in 1:n1) {
  for (j in 1:n1) {
    D12[i,j] =Eucl(s2[i,], s2[j,])
  }
}
D12
COV12 = kernel(D12,THETA[S,2*n+5]) 
is_symmetric_matrix(COV12)

# compute the distances between new and old
D1mix=matrix(NA,n1,n)

for (i in 1:n1){
  for (j in 1:n) {
    D1mix[i,j]=Eucl(s2[i,], s1[j,])
    
  }
}
D1mix
COV1mix = kernel(D1mix,THETA[S,2*n+5])



W1new=COV1mix%*%solve(COV11 + tau_full*diag(length(w1s_monitor)))%*%w1s_monitor


#I compute the new y with the parameters obtained through gibbs


x_B_new=rep(0,n1)
for (i in 1:n1){
  #primo quadrante
  if (s2[i,1]>=2.5 && s2[i,2]>=2.5) {x_B_new[i]=20}
  if (s2[i,1]>=2.5 && 0<=s2[i,2] && s2[i,2]<=2.5) {x_B_new[i]=15}
  if (0<=s2[i,1] && s2[i,1]<=2.5 && s2[i,2]>=2.5) {x_B_new[i]=12}
  if (0<=s2[i,1]&& s2[i,1]<=2.5 && 0<=s2[i,2] &&s2[i,2]<=2.5) {x_B_new[i]=11}
  
  #secondo quadrante
  if (s2[i,1]<=(-2.5) && s2[i,2]>=2.5) {x_B_new[i]=13}
  if (s2[i,1]<=(-2.5) && 0<=s2[i,2] && s2[i,2]<=2.5) {x_B_new[i]=7}
  if ((-2.5)<=s2[i,1] && s2[i,1]<=0 && s2[i,2]>=2.5) {x_B_new[i]=3}
  if ((-2.5)<=s2[i,1] && s2[i,1]<=0 && 0<=s2[i,2] && s2[i,2]<=2.5) {x_B_new[i]=10}
  
  #terzo quadrante
  if (s2[i,1]<=(-2.5) && s2[i,2]<=(-2.5)) {x_B_new[i]=30}
  if (s2[i,1]<=(-2.5) && (-2.5) <=s2[i,2] && s2[i,2]<=0) {x_B_new[i]=25}
  if ((-2.5)<=s2[i,1] && s2[i,1]<=0 && s2[i,2]<=(-2.5)) {x_B_new[i]=22}
  if ((-2.5)<=s2[i,1] && s2[i,1]<=0 && (-2.5)<=s2[i,2] && s2[i,2]<=0) {x_B_new[i]=11.5}
  
  #quarto quadrante
  if (s2[i,1]>=2.5 && s2[i,2]<=(-2.5)) {x_B_new[i]=54}
  if (s2[i,1]>=2.5 && (-2.5)<=s2[i,2] && s2[i,2]<=0) {x_B_new[i]=37}
  if (0<=s2[i,1] && s2[i,1]<=2.5 && s2[i,2]<=(-2.5)) {x_B_new[i]=32.5}
  if (0<=s2[i,1] && s2[i,1]<=2.5 && (-2.5)<=s2[i,2] && s2[i,2]<=0) {x_B_new[i]=10.5}}





x_B_new
Xnew=cbind(x_B_new[1:(n1)])

Xtot=cbind(rep(1,n1),Xnew)




ynew=b0_full*rep(1,n1)+b1_full*Xtot[,2]+a00*W0new+a10*Xtot[,2]*W0new+a11*Xtot[,2]*W1new


#Estimated Gaussian process with 101 monitors with the new GPs with the previous y
ggplot()+
  geom_tile(aes(x=s2[1:2601,1],y=s2[1:2601,2],fill=ynew[1:2601]))+
  scale_fill_gradient2(low = "green",mid="white",high = "red")+
  geom_point(aes(x=s1[,1],y=s1[,2],col=y,size=y))


#3d plot made with points , in red y

rgl::plot3d(x=s2[,1],y=s2[,2],z=ynew)
rgl::points3d(x=s1[,1],y=s1[,2],z=y,col=2)
rgl::points3d(x=sp[,1],y=sp[,2],z=yp,col=3)

#3d continuous plot, in red y
rgl::surface3d(x=grid,y=grid,z=matrix(ynew[1:2601],51,51),col="lightyellow",size=.1)
rgl::points3d(x=s1[,1],y=s1[,2],z=y,col="red",size=5)


ypval=ynew[2602:2651,]
MSE=(1/50)*sum((yp-ypval)^2)
MSE
MAE=(1/50)*sum(abs(yp-ypval))
MAE



######################things to be removed ... don t look 

















###I create new spatial data





###We add diag(length(y)) assuming that sigma2 is equal to 1

#This is the initial GPs, co in this case we should use COV00,COV01,COVmix with the prior phi
A0=COV0mix%*%solve(COV00 + 0.001*diag(length(y)))%*%y
B0=COV01-COV0mix%*%solve(COV00+ 0.001*diag(length(y)))%*%t(COV0mix)

mean_w0s=as.matrix(A0,51,51)

#Estimated Gaussian process with 101 monitors with the initial GPs

ggplot()+
  geom_tile(aes(x=G[,1],y=G[,2],fill=A0))+
  scale_fill_gradient2(low = "green",mid="white",high = "red")+
  geom_point(aes(x=s[,1],y=s[,2],col=y,size=y))

#3d plot made with points
rgl::plot3d(x=G[,1],y=G[,2],z=A0)
rgl::points3d(x=s[,1],y=s[,2],z=y,col=2)

#3d continuous plot
rgl::surface3d(x=grid,y=grid,z=matrix(A0,51,51),col="lightyellow",size=0.1)
rgl::points3d(x=s[,1],y=s[,2],z=y,col="red",size=5)


#Estimated Gaussian process with 101 monitors with the new GPs with the new y
ggplot()+
  geom_tile(aes(x=G[,1],y=G[,2],fill=A0new))+
  scale_fill_gradient2(low = "green",mid="white",high = "red")+
  geom_point(aes(x=s[,1],y=s[,2],col=ynew_finale,size=ynew_finale))


#3d plot made with points
rgl::plot3d(x=G[,1],y=G[,2],z=A0new)
rgl::points3d(x=s[,1],y=s[,2],z=y,col=2)

#3d continuous plot
rgl::surface3d(x=grid,y=grid,z=matrix(A0new,51,51),col="lightyellow",size=0.1)
rgl::points3d(x=s[,1],y=s[,2],z=y,col="red",size=5)
rgl::points3d(x=s[,1],y=s[,2],z=ynew_finale,col="blue",size=5)

#Estimated Gaussian process with 101 monitors with the new GPs with the previous y
ggplot()+
  geom_tile(aes(x=G[,1],y=G[,2],fill=A1new))+
  scale_fill_gradient2(low = "green",mid="white",high = "red")+
  geom_point(aes(x=s[,1],y=s[,2],col=y,size=y))


#3d plot made with points , in red y

rgl::plot3d(x=G[,1],y=G[,2],z=A1new)
rgl::points3d(x=s[,1],y=s[,2],z=y,col=2)

#3d continuous plot, in red y
rgl::surface3d(x=grid,y=grid,z=matrix(A1new,51,51),col="lightyellow",size=.1)
rgl::points3d(x=s[,1],y=s[,2],z=y,col="red",size=5)

