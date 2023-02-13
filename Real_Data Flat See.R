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

data=read.csv("Obs_pm25conc_region_NE_2018_training-1.csv")
populations=read.csv("us2021census.csv") #popolazione con coordinate centri città
density=read.csv("uscitypopdensityGIUSTO.csv",sep=';') #densità città and state, ma state abbreviati
#density=read.csv("uscitypopdensity.csv")

abb=read.csv("abbreviation.csv")
colnames(abb)=c("Esteso","Abb","State")
ok=merge(abb,populations,by="State") #ora state con acronico
ok=ok[,c(1:2,4:9)]
colnames(ok)=c("Code","State","City","Type","Counties","Population","Latitude","Longitude")
finale=merge(ok,density,by=c("State","City"))
finale=finale[,c(1:2,5,6:8,10)]
colnames(finale)=c("State","City","Counties","Population","Latitude","Longitude","Pop_density")
populations=finale

riga2=which(populations$Latitude>=min(data$Lat) & populations$Latitude<=max(data$Lat) & populations$Longitude>=min(data$Lon) & populations$Longitude<=max(data$Lon))
populations_area<-populations[riga2,]

riga=which(data$Date=='2018-01-01')
data_train<-data[riga,]
y=data_train$Conc


s1=as.matrix(data_train[,4:5])
plot(s1,col="red",pch=19,main="Monitor Location",xlab="Latitude",ylab="Longitude")coord_pop=as.matrix(populations_area[,5:6])
n_pop=as.numeric(dim(coord_pop)[1])
grid=seq(min(s1[,1]),max(s1[,1]),by=0.2)
grid2=seq(min(s1[,2]),max(s1[,2]),by=0.2)
G = s2 = expand.grid(grid,grid2)
s2=as.matrix(s2)
n=as.numeric(dim(s1)[1])
n1=as.numeric(dim(s2)[1])

distanza=matrix(NA,n,n_pop)
k=5
neighbor_ind=rep(0,k)
weights=rep(0,k)
density_x_B=rep(0,n)
for (i in 1:n){
  for (j in 1:n_pop){
    distanza[i,j]=Eucl(s1[i,],coord_pop[j,])}
  neighbor_ind = which(distanza[i,] %in% sort(distanza[i,])[1:k])
  weights=(1/(distanza[i,neighbor_ind]))/sum(1/distanza[i,neighbor_ind])
  density_x_B[i]=log(round(as.numeric(weights%*%populations_area$Pop_density[neighbor_ind])))
}

tot=as.data.frame(cbind(s1[,1],s1[,2],density_x_B))
x_B=as.numeric(tot$density_x_B)
X=cbind(rep(1,n),x_B)
DataconConc=cbind(tot,data_train$Conc)

distanza_new=matrix(NA,n1,n_pop)
k=5
neighbor_ind_new=rep(0,k)
weights_new=rep(0,k)
density_x_B_new=rep(0,n1)
for (i in 1:n1){
  for (j in 1:n_pop){
    distanza_new[i,j]=Eucl(s2[i,],coord_pop[j,])}
  neighbor_ind_new = which(distanza_new[i,] %in% sort(distanza_new[i,])[1:k])
  weights_new=(1/(distanza_new[i,neighbor_ind_new]))/sum(1/distanza_new[i,neighbor_ind_new])
  density_x_B_new[i]=log(round(as.numeric(weights_new%*%populations_area$Pop_density[neighbor_ind_new])))
}

tot_new=as.data.frame(cbind(s2[,1],s2[,2],density_x_B_new))
x_B_new=as.numeric(tot_new$density_x_B_new)
X_new=cbind(rep(1,n1),x_B_new)
Xtot=X_new

RR=0
mean_mare_new=mean(as.numeric(tot_new$density_x_B_new))
for (i in 1:(n1)){
  if (as.numeric((tot_new$V1)[i])<=(38)){
    if(as.numeric((tot_new$V2)[i])>=(-75)){
      tot_new$density_x_B_new[i]=mean_mare_new
      RR=rbind(RR,i)
    } 
  }
}
RR=as.vector(RR[2:length(RR)])
RR2=0
#as.numeric((tot_new$Var2)[i])<=(((max(grid2)-(-75))/(40.5-38))*(as.numeric(tot_new$Var1)[i]-38)+75)
for (i in 1:n1){
  if (as.numeric((tot_new$V1)[i])>=(38)){
    if (as.numeric((tot_new$V2)[i])>=((4/3)*(as.numeric((tot_new$V1)[i]))-(377/3))){
      tot_new$density_x_B_new[i]=mean_mare_new
      RR2=rbind(RR2,i)      
    }
  }
}
RR2=as.vector(RR2[2:length(RR2)])

x_B_new=as.numeric(tot_new$density_x_B_new)
X_new=cbind(rep(1,n1),x_B_new)
Xtot=X_new

###############################################################################
############################Prior Parameter####################################
###############################################################################



#PRIOR parameter for the coregionalization matrix


a00=exp(rnorm(1))
a10=rnorm(1)
a11=exp(rnorm(1))


#PRIOR for phi0-> discrete uniform


phi0_seq=as.vector(seq(0.1,0.9,0.1))
phi0_seq
phi0_prior=rep(1/length(phi0_seq),length(phi0_seq))
phi0_prior
phi0=sample(phi0_seq,1,replace=FALSE,phi0_prior)


#PRIOR for phi1-> discrete uniform

#phi1_seq=as.vector(seq(0.05,0.1,0.05))
phi1_seq=as.vector(seq(0.1,0.9,0.1))
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


#analogously, we did the same also for the second GP
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

nt=n-16
nv=16

sp=s1[(nt+1):n,]
yp=y[(nt+1):n]

W0S=W0S[1:nt]
W1S=W1S[1:nt]
y=y[1:nt]
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
    
    phi0_post[i]=log(1/length(phi0_seq))-log(sqrt(det(Fernando+diag(n)*.00001)))-((bpostW0)%*%solve(Fernando+diag(n)*.00001)%*%t(bpostW0))/(2)
    
    
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

#THETA_last=THETA[250:500,]
x11()
par(mfrow=c(2,3))
plot(ts(THETA_last[,1]),xlab='iterations',ylab='beta0',main='Traceplot of beta0') ###b0

plot(ts(THETA_last[,2]),xlab='iterations',ylab='beta1',main='Traceplot of beta1') ###b1
plot(ts(THETA_last[,3]),xlab='iterations',ylab='tau^2',main='Traceplot of tau^2') ###tau2

hist(THETA_last[,2*n+4],breaks=100,xlab='phi0',ylab='frequencies',main='Frequencies of phi0') ###phi0
hist(THETA_last[,2*n+5],breaks=100,xlab='phi1',ylab='frequencies',main='Frequencies of phi1') ###phi1


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

# coord for prediction
s2=rbind(s2,sp)
x11()
plot(s2)
points(sp,col=2)

dim(s1)
dim(s2)
n1=dim(s2)[1]

Xnew=c(x_B_new,x_B[(nt+1):136])
Xtot=cbind(rep(1,n1),Xnew)


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


## DO THE SAME WITH THE OTHER GP PROCESS


D1=matrix(NA,n,n)

for (i in 1:n) {
  for (j in 1:n) {
    D1[i,j] =Eucl(s1[i,], s1[j,])
  }
}
D1
COV11=kernel(D1,mean(THETA_last[,2*n+5]))

D12=matrix(NA,n1,n1)


for (i in 1:n1) {
  for (j in 1:n1) {
    D12[i,j] =Eucl(s2[i,], s2[j,])
  }
}
D12
COV12 = kernel(D12,mean(THETA_last[,2*n+5]))
is_symmetric_matrix(COV12)

# compute the distances between new and old
D1mix=matrix(NA,n1,n)

for (i in 1:n1){
  for (j in 1:n) {
    D1mix[i,j]=Eucl(s2[i,], s1[j,])
    
  }
}
D1mix
COV1mix = kernel(D1mix,mean(THETA_last[,2*n+5]))



W1new=COV1mix%*%solve(COV11 + tau_full*diag(length(w1s_monitor)))%*%w1s_monitor


#I compute the new y with the parameters obtained through gibbs


ynew=b0_full*rep(1,n1)+b1_full*Xtot[,2]+a00*W0new+a10*Xtot[,2]*W0new+a11*Xtot[,2]*W1new

#Estimated Gaussian process with 101 monitors with the new GPs with the previous y
y_new_centered =ynew[1:2760]-mean(ynew)
y_centered=y-mean(y)
y_centered_size=y-mean(y)

ggplot()+
  geom_tile(aes(x=G[,1],y=G[,2],fill=y_new_centered))+
  scale_fill_gradient2(low = "green",mid="white",high = "red")+
  geom_point(aes(x=s1[,1],y=s1[,2],col=y_centered,size=y_centered_size))


#3d continuous plot, in red y
z=as.numeric(ynew[1:2760])
ypval=ynew[2761:2776]
zz= (z-min(z))/diff(range(z))
x= colorRamp(c("lightgreen","white","red"))(zz)
#x= colorRamp(c("red","lightblue","purple"))(zz)
colori=rgb(x[,1],x[,2],x[,3],maxColorValue =255)

#rgl::surface3d(x=grid,y=grid2,z=matrix(ynew[1:2760],40,69),col="lightyellow",size=.1)
rgl::plot3d(x=s1[,1],y=s1[,2],z=y,col="red",size=5)
rgl::surface3d(x=grid,y=grid2,z=matrix(ynew[1:2760],40,69),col=colori,size=.1)
for(k in 1:length(yp)){
  rgl::lines3d(x=sp[k,1],y=sp[k,2],z=c(yp[k],ypval[k]),col="darkgreen",size=30)
}

rgl::points3d(x=sp[,1],y=sp[,2],z=yp,col=3,size=5)
rgl::points3d(x=sp[,1],y=sp[,2],z=ypval,size=5)


MSE=(1/16)*sum((yp-ypval)^2)
MSE
MAE=(1/16)*sum(abs(yp-ypval))
MAE
