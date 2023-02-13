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

#populations=read.csv("us2021census.csv")
data=read.csv("Obs_pm25conc_region_NE_2018_training-1.csv")
#riga2=which(populations$Latitude>=min(data$Lat) & populations$Latitude<=max(data$Lat) & populations$Longitude>=min(data$Lon) & populations$Longitude<=max(data$Lon))
#populations_area<-populations[riga2,]

riga=which(data$Date=='2018-01-01')
data_train<-data[riga,]
y=data_train$Conc


s1=as.matrix(data_train[,4:5])
x11()
plot(s1)
grid=seq(min(s1[,1]),max(s1[,1]),by=0.2)
grid2=seq(min(s1[,2]),max(s1[,2]),by=0.2)
G = s2 = expand.grid(grid,grid2)
s2=as.matrix(s2)
plot(s2)
n=as.numeric(dim(s1)[1])
n1=as.numeric(dim(s2)[1])

distanza=matrix(NA,n,n)

for (i in 1:n){
  for (j in 1:n){
    distanza[i,j]=Eucl(s1[i,],s1[j,])
    }
  }
distanza
D0=distanza

distanza_new=matrix(NA,n1,n1)
for (i in 1:n1){
  for (j in 1:n1){
    distanza_new[i,j]=Eucl(s2[i,],s2[j,])
  }
}


###############################################################################
############################Prior Parameter####################################
###############################################################################

#Y= b0 +b0(S) +er

#PRIOR parameter for the coregionalization matrix


a00=exp(rnorm(1))
a10=rnorm(1)
a11=exp(rnorm(1))


#PRIOR for phi0-> discrete uniform


phi0_seq=as.vector(seq(0.01,0.1,0.01))
phi0_seq
phi0_prior=rep(1/length(phi0_seq),length(phi0_seq))
phi0_prior
phi0=sample(phi0_seq,1,replace=FALSE,phi0_prior)
COV00=kernel(D0,phi0)

COV00_v=list()
for (i in 1:length(phi0_seq)){
  COV00_v[[i]] = kernel(D0,phi0_seq[i])}

#COV00_2 = kernel(D0,phi0_seq[2])
#COV00_3 = kernel(D0,phi0_seq[3])
#COV00_4 = kernel(D0,phi0_seq[4])
#COV00_5 = kernel(D0,phi0_seq[5])
#COV00_6 = kernel(D0,phi0_seq[6])

#analogously, we did the same also for the new point1)


#PRIOR OF W0S
b0s=rep(0,n)
COV00
W0S= t(rmvnorm(1, rep(0,n), COV00+diag(n)*.00001))


#PRIOR OF TAU2
a1=4
b1=0.2
tau2=rinvgamma(1,a1,b1)  


#PRIOR OF BETA => bivariate normal (mu, sigmaB)

b_mu_prior=0
sigma=2
sigmaB=sigma^-1



b=rnorm(1,b_mu_prior,sigma) ##they will be the first value of the Gibbs Sampler

b



###############################################################################
######################Value of our simulated y#################################
###############################################################################


y_mu=b[1]*rep(1,n)+a00*W0S

y = data_train$Conc

y

#z= b[1]*rep(1,n) + b[2]*X[,2]+ W0S+eps
#eps=rnorm(n,0,tau2)



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
THETA<-matrix(nrow=S,ncol=n+3)


# Initial point of the chain 
theta<-cbind(b[1],tau2,t(W0S),phi0) # Initial point of the chain 
theta

THETA[1,]=theta


X=rep(1,n)
sost1=t(X)%*%X
C=diag(a00*rep(1,n))
sost2=t(C)%*%C


################################
####### Gibbs sampling #########

##THETA (bo,tau2,wos1,wos2,..wosN,phi0)

r=2
for(r in 2:500) {
  
  # generate a new beta value from its full conditional
  if(r%%20==0){
    cat(paste(r,"\n"))
  }
  
  
  betahat= solve(sost1)%*%t(X)%*%(y-a00*THETA[r-1,3:(n+2)])#da mettere i valori dentro THETA
  sigmaBpost=  solve((sost1)/THETA[r-1,2] +sigmaB)
  muPost= sigmaBpost %*% ((( (sost1)/THETA[r-1,2])%*%betahat) + sigmaB%*%b_mu_prior)
  
  
  bpost= rmvnorm(1, muPost, sigmaBpost)
  bpost
  theta[1]=bpost[1]
  
  # generate a new tau^2 from its full conditional
  mu=X%*%t(bpost)+a00*THETA[r-1,3:(n+2)]
  mu
  an= a1+ n/2
  bn= b1+ sum((y - mu)^2)/2 
  
  theta[2]<- rinvgamma(1, an, bn)
  
  
  ## generate a new W0S from its full conditional
  that=solve(sost2)%*%t(C)%*%(y-X%*%t(bpost))
  that
  
  sigmaBpostW0=solve(((sost2)/ theta[2]+diag(n)*.00001)+solve(COV00+diag(n)*.00001))
  
  muPostW0= sigmaBpostW0%*% (((sost2)/ theta[2]+diag(n)*.00001)%*%that + solve(COV00+diag(n)*.00001)%*%b0s)
  
  
  bpostW0= rmvnorm(1, muPostW0, sigmaBpostW0+diag(n)*.00001)
  
  
  for (j in 1:n) {
    theta[j+2]=bpostW0[j]
  }
  
  
  ## generate a new phi0 from its full conditional
  
  phi0_post=rep(0,length(phi0_seq))
  for (i in 1:length(phi0_seq)) {
    Fernando=COV00_v[[i]]
    
    C=diag(a00*rep(1,n))
    that=solve(sost2)%*%t(C)%*%(y-X%*%t(bpost)) 
    
    sigmaBpostW0=solve(((sost2)/ theta[2]+diag(n)*.00001)+solve(Fernando+diag(n)*.00001))
    
    
    muPostW0= sigmaBpostW0%*% (((sost2)/ theta[2]+diag(n)*.00001)%*%that + solve(Fernando+diag(n)*.00001)%*%b0s)
    
    bpostW0= rmvnorm(1, muPostW0, sigmaBpostW0+diag(n)*.00001)
    
    phi0_post[i]=log(1/length(phi0_seq))-log(sqrt(det(Fernando+diag(n)*.00001)))-((bpostW0)%*%solve(Fernando+diag(n)*.00001)%*%t(bpostW0))/(2)
    
    
  }
  
  phi0_post=exp(phi0_post)/sum(exp(phi0_post))
  theta[n+3]=sample(phi0_seq,1,replace=FALSE,phi0_post)
  
  
  COV00=COV00_v[[which(phi0_seq==theta[n+3])]]
  

  theta
  THETA[r,]<-t(theta)
  THETA
}


plot(phi1_post)
plot(phi0_post)

THETA



##traceplots and histogram of the parameters of which we have assumed a prior

x11()
par(mfrow=c(1,3))
plot(ts(THETA[,1])) 
plot(ts(THETA[,2])) 
hist(THETA[,n+3],breaks = 100)




###we decided to take the last 200 iterations and to take the average of them as final value

THETA_last=THETA[(2):500,]

plot(ts(THETA_last[,1])) ###b0
plot(ts(THETA_last[,2])) ###tau2

hist(THETA_last[,n+3],breaks = 100) ###phi0



##this are our final result of the parameters
b0_full=mean(THETA_last[,1])
b0_full
tau_full=mean(THETA_last[,2])
tau_full

w0s_monitor=colMeans(THETA_last[,3:(n+2)])




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
COV00=kernel(D0,mean(THETA_last[,n+3]))

D01=matrix(NA,n1,n1)


for (i in 1:n1) {
  for (j in 1:n1) {
    D01[i,j] =Eucl(s2[i,], s2[j,])
  }
}
D01
COV01 = kernel(D01,mean(THETA_last[,n+3])) 
is_symmetric_matrix(COV01)

# compute the distances between new and old
D0mix=matrix(NA,n1,n)

for (i in 1:n1){
  for (j in 1:n) {
    D0mix[i,j]=Eucl(s2[i,], s1[j,])
    
  }
}
D0mix
COV0mix = kernel(D0mix,mean(THETA_last[,n+3]))


W0new=COV0mix%*%solve(COV00 + tau_full*diag(length(w0s_monitor)))%*%w0s_monitor




ynew=b0_full*rep(1,n1)+a00*W0new
#for(i in 1:n1){
# if(ynew[i]<(-10^3)){ynew[i]=0}}

#Estimated Gaussian process with 101 monitors with the new GPs with the previous y
ggplot()+
  geom_tile(aes(x=G[,1],y=G[,2],fill=ynew))+
  scale_fill_gradient2(low = "green",mid="white",high = "red")+
  geom_point(aes(x=s1[,1],y=s1[,2],col=y,size=y))


#3d plot made with points , in red y

rgl::plot3d(x=G[,1],y=G[,2],z=ynew,xlim=c(min(grid),max(grid)),ylim=c(min(grid2),max(grid2)))
rgl::plot3d(x=s1[,1],y=s1[,2],z=y,col=2)

#3d continuous plot, in red y
rgl::surface3d(x=grid,y=grid2,z=matrix(ynew,40,69),col="lightyellow",size=.1)
rgl::points3d(x=s1[,1],y=s1[,2],z=y,col="red",size=5)






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

