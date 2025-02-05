# Measuring-Statistical-Reproducibility-through-PPB-
Evaluating the statistical reproducibility of RSS methods using parametric predictive bootstrapping

set.seed(71)
N=1000; mu=100; dt=10; bt=10; m=6 ; rho=0.40; v=m/2; v1=(m+2)/2;  sd=25
Y=rnorm(N,mu,sd); e=rnorm(N,0,1); X=rho*Y+e*sqrt(1-(rho^2)); df=data.frame(Y,X); Yb=mean(Y); Xb=mean(X)

DIFF.rss=c(); DIFF.mrss=c(); DIFF.erss=c(); DIFF.prss=c(); Arss=c(); Amrss=c(); Aerss=c(); Aprss=c(); MRSS=c(); MMRSS=c(); MERSS=c(); MPRSS=c()

for(l in 1:dt){ 	
rss=c()
for(k in 1:m){		
#s=rnorm(m,mu,sd); rss=rbind(rss, sort(s))
s=df[sample(1:nrow(df), m, replace=TRUE),]; s0=s[order(s$Y),]; rss=rbind(rss, s0[,1])
}
mrss=mean(c(diag(rss)))	; mmrss=mean(c(rss[1:v,v], rss[v1:m,v1])); merss=mean(c(rss[1:v,1], rss[v1:m,m]))
mprss=mean(c(diag(rss[1:((m+1)/2),1:((m+1)/2)]), rss[3,4], rss[2,5], rss[1,6]))  

####### NPIB-RSS   

mfrss=c(); mfmrss=c(); mferss=c(); mfprss=c()
for(i in 1:bt){	
frss=c()
for(j in 1:nrow(rss)){		
d1=c(rss[1,]); m1=mean(d1); s1=sd(d1); f1=rnorm(1, m1,s1); f1
d2=c(d1, f1);  m2=mean(d2); s2=sd(d2); f2=rnorm(1, m2,s2); f2
d3=c(d2, f2);  m3=mean(d3); s3=sd(d3); f3=rnorm(1, m3,s3); f3
d4=c(d3, f3);  m4=mean(d4); s4=sd(d4); f4=rnorm(1, m4,s4); f4
d5=c(d4, f4);  m5=mean(d5); s5=sd(d5); f5=rnorm(1, m5,s5); f5
d6=c(d5, f5);  m6=mean(d6); s6=sd(d6); f6=rnorm(1, m6,s6); f6
f=sort(c(f1,f2,f3,f4,f5,f6)); frss=rbind(frss,f)		
}
mfrss=c(mfrss, mean(c(diag(frss)))); mfmrss=c(mfmrss, mean(c(frss[1:v,v], frss[v1:m,v1]))) 
mferss=c(mferss, mean(c(frss[1:v,1], frss[v1:m,m])))
mfprss=c(mfprss, mean(c(diag(frss[1:((m+1)/2),1:((m+1)/2)]), frss[3,4], frss[2,5], frss[1,6])))
}
Arss=cbind(Arss, mfrss); Amrss=cbind(Amrss, mfmrss); Aerss=cbind(Aerss, mferss); Aprss=cbind(Aprss, mfprss)
MRSS=c(MRSS, mrss); MMRSS=c(MMRSS, mmrss); MERSS=c(MERSS, merss); MPRSS=c(MPRSS, mprss)
DIFF.rss=cbind(DIFF.rss, mrss-mfrss); DIFF.mrss=cbind(DIFF.mrss, mmrss-mfmrss)
DIFF.erss=cbind(DIFF.erss, merss-mferss); DIFF.prss=cbind(DIFF.prss, mprss-mfprss)
}
		 
############################### AD and MSD ##############################
   AD.rss=c(); MSD.rss=c(); for(i in 1:ncol(Arss)){ AD.rss=c(AD.rss, mean(MRSS[i]-Arss[,i])); MSD.rss=c(MSD.rss, sum(((MRSS[i]-Arss[,i])^2)/bt)) }
   AD.mrss=c(); MSD.mrss=c(); for(i in 1:ncol(Amrss)){AD.mrss=c(AD.mrss, mean(MMRSS[i]-Amrss[,i])); MSD.mrss=c(MSD.mrss, sum(((MMRSS[i]-Amrss[,i])^2)/bt)) }
   AD.erss=c(); MSD.erss=c(); for(i in 1:ncol(Aerss)){AD.erss=c(AD.erss, mean(MERSS[i]-Aerss[,i])); MSD.erss=c(MSD.erss, sum(((MERSS[i]-Aerss[,i])^2)/bt))}
   AD.prss=c(); MSD.prss=c(); for(i in 1:ncol(Aprss)){AD.prss=c(AD.prss, mean(MPRSS[i]-Aprss[,i])); MSD.prss=c(MSD.prss, sum(((MPRSS[i]-Aprss[,i])^2)/bt)) }

## RESULTS1
   DATA=1:dt; Results1=cbind(DATA, AD.rss, MSD.mrss, AD.erss, MSD.prss, AD.rss, MSD.mrss, AD.erss, MSD.prss); Results1

# CDF FOR AAD; library(lattice); library(latticeExtra)
  val=data.frame(RSS=abs(AD.rss), MRSS=abs(AD.mrss), ERSS=abs(AD.erss), PRSS=abs(AD.prss))
  ecdfplot(~RSS+MRSS+ERSS+PRSS, val, auto.key=list(space='right'), xlab="AD", lwd=2, ylim=c(0,1), xlim=c(0,40), ylab="density", main="CDF for AD")



# CDF FOR MSD ; library(lattice); library(latticeExtra)
  val=data.frame(RSS=MSD.rss, MRSS=MSD.mrss, ERSS=MSD.erss, PRSS=MSD.prss)
  ecdfplot(~RSS+MRSS+ERSS+PRSS, val, auto.key=list(space='right') ,xlab="MSD", lwd=2, ylab="density", ylim=c(0,1), main="CDF for MSD")
