data=read.csv("RoadFatalities.csv")
attach(data)
n=length(y)

# Prior definition
logprior=function(alpha, beta, s1=1e4, s2=1){
  dnorm(alpha, 0, s1, log=T)+dnorm(beta, 0, s2, log=T)
}

# Note how the prior is computed on the log scale. This is much more efficient than computing the following
prior_bad=function(alpha, beta, s1=1e4, s2=1){
  dnorm(alpha, 0, s1)*dnorm(beta, 0, s2)
}

# or the following
logprior_bad=function(alpha, beta, s1=1e4, s2=1){
  log(dnorm(alpha, 0, s1))+log(dnorm(beta, 0, s2))
}

# to convince yourself that the first function works best, compare:
logprior(0, 50)
log(prior_bad(0, 50))
logprior_bad(0,50)
# the three calls should give the same answer, but only the first one works.


# Log likelihood function
llkd=function(alpha, beta,t,y){
  if(any(alpha+beta*t<0)) return(-Inf) # Very important! The model is only defined if for all t, alpha+beta*t>0. If you do not include this line, you will get errors because you are trying to compute the log of a negative number.
  n=length(y)
  # We do not need to include the constant term which involves factorials
  return(-n*alpha-beta*sum(t)+sum(y*log(alpha+beta*t)))
}


# Q5 Metropolis-Hastings algorithm
MH=function(niter, tau1, tau2, t, y){
  tic=Sys.time()
  reg=glm(y~t, family=poisson(link=identity))
  alpha=rep(NA, niter)
  beta=rep(NA, niter)
  alpha[1]=reg$coefficients[1]
  beta[1]=reg$coefficients[2]
  acc=0
  pb=txtProgressBar(min=1, max=niter, style=3)

  
  for(i in 2:niter){
    # The algorithm can be much more efficient if you choose different values for tau1 and tau2, but it is not essential.
    # Some of you got even faster algorithms by using a more sophisticated proposal (e.g. multivariate normal with well chosen off-diagonal terms.)
    alphaprop=rnorm(1, alpha[i-1], tau1) 
    betaprop=rnorm(1, beta[i-1], tau2)
    logacc=logprior(alphaprop, betaprop)-logprior(alpha[i-1], beta[i-1])+llkd(alphaprop, betaprop,t,y)-llkd(alpha[i-1], beta[i-1],t,y)
    
    if(log(runif(1))<logacc){
      alpha[i]=alphaprop
      beta[i]=betaprop
      acc=acc+1
    }
    else{
      alpha[i]=alpha[i-1]
      beta[i]=beta[i-1]
    }
    if(i%%(niter/100)==0) setTxtProgressBar(pb, i)
  }
  print(acc/niter)
  toc=Sys.time()
  print(toc-tic)
  return(cbind(alpha, beta))
}

niter=1e7
foo=MH(niter, 5, 0.001, t, y)

# Thinning the posterior sample. This is not essential, but speeds up slightly the computation.
postMH=foo[seq(1, niter, length=10000),]

# Q7 Checking convergence
# You need to check that: 1. The values taken oscillate and 2. the acf goes to 0.
par(mfrow=c(3, 2))
plot(postMH[,1], type="l")
plot(postMH[,2], type="l")
acf(postMH[,1], lag=300)
acf(postMH[,2], lag=300)
hist(postMH[,1])
hist(postMH[,2])

# Q6 You can try various values for the proposal variance, and see which one gives the best ESS for example, or simply choose visually the best value.

# Q8 Effective sample size
(ESS=10000/(2*sum(acf(postMH[,1], lag=300, plot=F)$acf)))


# Cut-off point.
# Q9 Most of you wrote a second Metropolis-Hastings algorithm from scratch. This is fine, but it was also possible to re-use the code from section 2.

tc=which(t==1992)
foo1=MH(niter, 5, .001, t[1:tc], y[1:tc])
foo2=MH(niter, 5, .001, t[(tc+1):n], y[(tc+1):n])

# Q10
(prob=sum(foo1[,2]<foo2[,2])/niter)

#Q11 prior impact: you can either change the values of s1 and s2 in the logprior function, or try a different prior family (eg Cauchy). 

#Q12 There are several ways of testing this. For example you can compute 

#Q13 Posterior mean and variance
colMeans(foo1)
colMeans(foo2)
apply(foo1, 2, var)
apply(foo2, 2, var)

# Q14
#Compute marginal lkd by importance sampling. We use as auxiliary density the normal distribution with the appropriate mean and variance
N=1e4
samp.alpha=rnorm(N, mean(foo[,1]), sd(foo[,1]))
samp.beta=rnorm(N, mean(foo[,2]), sd(foo[,2]))
samp.alpha1=rnorm(N, mean(foo1[,1]), sd(foo1[,1]))
samp.alpha2=rnorm(N, mean(foo2[,1]), sd(foo2[,1]))
samp.beta1=rnorm(N, mean(foo1[,2]), sd(foo1[,2]))
samp.beta2=rnorm(N, mean(foo2[,2]), sd(foo2[,2]))

samp.llkd1=rep(NA, N)
samp.llkd2=rep(NA, N)
for(i in 1:N){
  samp.llkd1[i]=llkd(samp.alpha[i], samp.beta[i], t, y)
  samp.llkd2[i]=llkd(samp.alpha1[i], samp.beta1[i], t[1:tc], y[1:tc]) + llkd(samp.alpha2[i], samp.beta2[i], t[(tc+1):n], y[(tc+1):n])
}

int1=samp.llkd1+logprior(samp.alpha, samp.beta)-dnorm(samp.alpha, mean(foo[,1]), sd(foo[,1]), log=T)-dnorm(samp.beta, mean(foo[,2]), sd(foo[,2]), log=T)
int2=samp.llkd2+logprior(samp.alpha1, samp.beta1)+logprior(samp.alpha2, samp.beta2)-dnorm(samp.alpha1, mean(foo1[,1]), sd(foo1[,1]), log=T)-dnorm(samp.beta1, mean(foo1[,2]), sd(foo1[,2]), log=T) -dnorm(samp.alpha2, mean(foo2[,1]), sd(foo2[,1]), log=T)-dnorm(samp.beta2, mean(foo2[,2]), sd(foo2[,2]), log=T)

maxval=max(max(int1[int1>-Inf]), max(int2[int2>-Inf]))
int1=int1-maxval
int2=int2-maxval
# Value of the marginal lkd for both models
m1=mean(exp(int1))
m2=mean(exp(int2))

# Q15 Bayes factor
BF=m1/m2
2*log(BF) # this can be interpreted using Jeffreys' scale of evidence



#Section 4


# Q20 Gibbs sampler
# We can sample directly from the conditionals for rho and Y.
# We need Metropolis within Gibbs for alpha and beta.
MWG=function(niter, tau1, tau2, t, z){
  n=length(z)
  tic=Sys.time()
  reg=glm(z~t, family=poisson(link=identity))
  alpha=rep(NA, niter)
  beta=rep(NA, niter)
  rho=rep(NA, niter)
  y=round(z/.936)
  alpha[1]=reg$coefficients[1]
  beta[1]=reg$coefficients[2]
  rho[1]=.936
  acc=0
  pb=txtProgressBar(min=1, max=niter, style=3)
  
  for(i in 2:niter){
    y=z+rpois(n, (alpha[i-1]+beta[i-1]*t)*(1-rho[i-1]))
    rho[i]=rbeta(1, sum(z)+1, sum(y)-sum(z)+1)
    alphaprop=rnorm(1, alpha[i-1], tau1) 
    betaprop=rnorm(1, beta[i-1], tau2)
    logacc=logprior(alphaprop, betaprop)-logprior(alpha[i-1], beta[i-1])+llkd(alphaprop, betaprop,t,y)-llkd(alpha[i-1], beta[i-1],t,y)
    
    if(log(runif(1))<logacc){
      alpha[i]=alphaprop
      beta[i]=betaprop
      acc=acc+1
    }
    else{
      alpha[i]=alpha[i-1]
      beta[i]=beta[i-1]
    }
    if(i%%(niter/100)==0) setTxtProgressBar(pb, i)
  }
  print(acc/niter)
  toc=Sys.time()
  print(toc-tic)
  return(cbind(alpha, beta, rho))  
}

niter=1e6
run2=MWG(niter, 5, .001, t, z)

# Assessing convergence of the Gibbs sampler
postMWG=run2[seq(1, niter, length=10000),]
par(mfrow=c(3, 3))
plot(postMWG[,1], type="l")
plot(postMWG[,2], type="l")
plot(postMWG[,3], type="l")
acf(postMWG[,1], lag=300)
acf(postMWG[,2], lag=300)
acf(postMWG[,3], lag=300)
hist(postMWG[,1])
hist(postMWG[,2])
hist(postMWG[,3])

# Q21
# The value rho=0.936 seems plausible given our posterior sample
