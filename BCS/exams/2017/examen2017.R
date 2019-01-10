#Préparation des données - inutile pour les étudiants
mut=barresMutationsEnseignants
bac=resultats_bac

d=merge(mut, bac, by.x="CodeEtab", by.y="code_etablissement", all.x=T)
#d2=merge(d, DNB, by.x="NomEtab", by.y="nom", all.x=T)
d=d[,-c(8:11)]

write.csv(d, file="examdata2017.csv")


#Q2
dpareto=function(x,m,alpha, log=F){
  if(log) return(log(alpha)+alpha*log(m)-(alpha+1)*log(x) +log(x>=m))
  return((x>=m)*alpha*m^alpha/(x^(alpha+1)))
}

rpareto=function(n,m,alpha){
  return(m*exp(rexp(n, alpha)))
}

#Q5
z=d$Barre
n=length(z)
m=21
a=2
b=2

apost=a+n
bpost=b+sum(log(z/m))
apost/bpost #posterior mean
apost/bpost^2 #posterior variance

#Q6
# Marginal lkd for a Gamma(a,b) prior
logmlpareto=function(x, m, a, b){
  n=length(x)
  return(a*log(b)-lgamma(a)-n*log(m)-sum(log(x))+lgamma(a+n)-(a+n)*log(b+sum(x/m)))
}

#Q7
zm=z[d$Matiere=="MATHS"]
za=z[d$Matiere=="ANGLAIS"]
(BF=logmlpareto(c(zm,za),m,a,b)-logmlpareto(zm,m,a,b)-logmlpareto(za,m,a,b))
# decisive evidence in favour of model 2

#Q9
alpha=.5
# posterior mean
min(z)*n*alpha/(n*alpha+1)
# posterior variance
min(z)^2*(n*alpha/(n*alpha+2)-(n*alpha/(n*alpha+1))^2)


#Q14
gibbsmix=function(z, niter, a=2, b=2){
  n=length(z)
  alpha1=rep(NA, niter)
  alpha2=rep(NA, niter)
  m1=rep(NA, niter)
  m2=rep(NA, niter)
  W=sample(2, n, rep=T)
  p=rep(NA, niter)
  
  p[1]=.5
  alpha1[1]=.5
  alpha2[1]=.5
  m1[1]=21
  m2[1]=800
  
  for(i in 2:niter){
    for(j in 1:n) W[j]=sample(2, 1, prob=c(p[i-1]*dpareto(z[j], m1[i-1], alpha1[i-1], log=F), (1-p[i-1])*dpareto(z[j], m2[i-1], alpha2[i-1], log=F)))
    n1=sum(W==1)
    n2=sum(W==2)
    alpha1[i]=rgamma(1, a+n1, b+sum(log(z[W==1]/m1[i-1])))
    alpha2[i]=rgamma(1, a+n2, b+sum(log(z[W==2]/m2[i-1])))
    m1[i]=runif(1)^(1/(n1*alpha1[i]))*min(z[W==1])
    m2[i]=runif(1)^(1/(n2*alpha2[i]))*min(z[W==2])
    p[i]=rbeta(1, n1+1, n2+1)
  }
  plot(m2, type="l")
  return(matrix(c(alpha1, alpha2, m1, m2, p), nrow=niter))
}

niter=1e3
gm=gibbsmix(z,niter)

ESS=niter/(2*sum(acf(gm[, 4], plot=F)$acf)-1)

lkdmix=function(x,alpha1, alpha2, m1, m2, p){
  return((p<1)*( p*dpareto(x, m1, alpha1) + (1-p)*dpareto(x, m2, alpha2)))
}

logpostmix=function(x,alpha1, alpha2, m1, m2, p, a=2, b=2){
  return(sum(log(lkdmix(x, alpha1, alpha2, m1, m2, p))) + dgamma(alpha1, a, b, log=T) + dgamma(alpha2, a, b, log=T) -log(m1) -log(m2))
}

#estimate log marginal lkd by importance sampling
require(mvtnorm)
postmean=colMeans(gm)
postvar=cov(gm)
NIS=1e2
param=rmvnorm(NIS, postmean, postvar) #we could make this more efficient by ensuring that m1<21 and p<1.
s=rep(0, NIS)
for( i in 1:NIS){
  s[i]=logpostmix(z, param[i,1], param[i,2], param[i,3], param[i,4], param[i,5])-dmvnorm(param[i,], postmean, postvar, log=T)
}

#Q19 Adapt code for mode with no mixture
gibbsnomix=function(z, niter, a=2, b=2){
  n=length(z)
  alpha1=rep(NA, niter)
  m1=rep(NA, niter)
  
  alpha1[1]=.5
  m1[1]=21

  for(i in 2:niter){
    alpha1[i]=rgamma(1, a+n, b+sum(log(z/m1[i-1])))
    m1[i]=runif(1)^(1/(n*alpha1[i]))*min(z)
  }
  return(matrix(c(alpha1, m1), nrow=niter))
}

niter=1e3
gm2=gibbsnomix(z,niter)

niter/(2*sum(acf(gm[, 2], plot=F)$acf)-1)

logpostnomix=function(x,alpha1,  m1, a=2, b=2){
  return(sum(dpareto(x, m1, alpha1, log=T)) + dgamma(alpha1, a, b, log=T) -log(m1))
}

#estimate log marginal lkd by importance sampling
postmean2=colMeans(gm2)
postvar2=cov(gm2)
NIS=1e2
param2=rmvnorm(NIS, postmean2, postvar2)
s2=rep(0, NIS)
for( i in 1:NIS){
  s2[i]=logpostnomix(z, param2[i,1], param2[i,2]-dmvnorm(param2[i,], postmean2, postvar2, log=T))
}

mean(exp(s))
maxs=max(max(s), max(s2))
s=s-maxs
s2=s2-maxs
mean(exp(s))/mean(exp(s2))


# Section on generalized Pareto
#Q21
dgpd=function(x, m, alpha, tau, log=F){
  logd=log(alpha)-log(m)-log(tau)-(alpha+1)*log(1+(x-m)/(tau*m))+log(x>=m)
  if(log) return(logd)
  return(exp(logd))
}

logpostgpd=function(x, m, alpha, tau, a=2, b=2){
  return(sum(dgpd(x, m, alpha, tau, log=T)) + dgamma(alpha, a, b, log=T) + dexp(tau, 1,  log=T))
}

MH=function(data, niter, sigmaalpha=1, sigmatau=1, alphastart=1, taustart=1){
  m=21
  alpha=rep(NA, niter)
  tau=rep(NA, niter)
  
  alpha[1]=alphastart
  tau[1]=taustart
  
  for(i in 2:niter){
    #our proposal kernel updates either alpha or tau
    if(runif(1)<.5){
      alphaprop=alpha[i-1]
      tauprop=rnorm(1, tau[i-1], sigmatau)
    }
    else{
      tauprop=tau[i-1]
      alphaprop=rnorm(1, alpha[i-1], sigmaalpha)
    }
    
    if(tauprop<0 | alphaprop <0){
      acc=0
    }
    else{
      acc=logpostgpd(data, m, alphaprop, tauprop) - logpostgpd(data, m, alpha[i-1], tau[i-1])
    }
    if(log(runif(1))<acc){
      alpha[i]=alphaprop
      tau[i]=tauprop
    }
    else{
      alpha[i]=alpha[i-1]
      tau[i]=tau[i-1]
    }
  }
  
  return(matrix(c(alpha, tau), nrow=niter))
}

niter=1e4
smh=MH(z, niter, .1, .1, .4,.2)

par(mfrow=c(2,1))
plot(smh[,1], type="l")
plot(smh[,2], type="l")

niter/(2*sum(acf(smh[,2], plot=F)$acf)-1)

quantile(smh[,2], c(.025, .975))

#Q25 Validation
v=rpareto(1000, 21, .5)
sval=MH(v,  niter, .01 , 0.01, 1, 1)
quantile(sval[,2], c(.025, .975))

#Q26 marginal lkd
postmean3=colMeans(smh)
postvar3=var(smh)
NIS=1e3
param3=rmvnorm(NIS, postmean3, postvar3)
s3=rep(0, NIS)
for( i in 1:NIS){
  s3[i]=logpostgpd(z, 21, param3[i,1], param3[i,2])-dmvnorm(param3[i,], postmean3, postvar3, log=T)
}
#Q27 BF
smh2=MH(z, niter, .1,0,.4,1) #set sigmatau=0 to sample from the model with fixed tau
postmean4=mean(smh2[,1])
postvar4=var(smh2[,1])
NIS=1e3
param4=rnorm(NIS, postmean4, sqrt(postvar4))
s4=rep(0, NIS)
for( i in 1:NIS){
  s4[i]=logpostgpd(z, 21, param4[i], 1)-dnorm(param4[i], postmean4, sqrt(postvar4), log=T)
}
