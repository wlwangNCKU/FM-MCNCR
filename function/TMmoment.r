library(mvtnorm)
library(tmvtnorm)
library(cubature)

# generate random sample from truncated multivariate normal independent distribution
rtmni = function(n, mu, Sigma, nu=NULL, rho=NULL, distr=c('MVN','MVT','MCN'), lower=rep(-Inf,length(mu)), upper=rep(Inf, length(mu)))
{
  begin = proc.time()[1]
  distr = distr[1]
  p = length(mu)
  if(p == 1){
   Sd = sqrt(Sigma)
  } else{
  Lam = diag(sqrt(diag(Sigma)))
  Lam.inv = diag(1/sqrt(diag(Sigma)))
  R = Lam.inv %*%Sigma%*% Lam.inv 
  }
#  cat('Generate random samples of size', n, 'from the', p,'-variate', distr, 'distribution', '\n')
  Y = matrix(NA, nrow=n, ncol=p)
  Tau = numeric(n)
  iter = 1
  while(iter <= n){
   if(distr=='MVN') tau = 1  
   if(distr=='MVT') tau = rgamma(1, shape=nu/2, rate=nu/2)
   if(distr=='MCN'){
     a = runif(1)
     tau = ifelse(a<=nu, rho, 1)
   }
   tau.sq = sqrt(tau)
   if(p==1){
     yi = mu + 1/tau.sq * Sd * rnorm(1, mean=0, sd=1) 
   } else  yi = mu + 1/tau.sq * Lam %*% c(rmvnorm(1, mean=rep(0,p), sigma=R))
   if((sum(yi>=lower)+sum(yi<=upper))==(2*p)){
     Y[iter, ] = yi
     Tau[iter] = tau
     iter = iter + 1
   } else next
  }
  end = proc.time()[1]
  run.sec = end - begin
  return(list(Y=Y, tau=Tau))
}

# PDF of multivariate normal-independent distribution
dmni = function(x, mu=rep(0,length(mu)), Sigma=diag(length(mu)), nu=NULL, rho=NULL, distr=c('MVN','MVT','MCN'))
{
   distr = distr[1]
   p = length(mu)
   if(distr=='MVN'){
     if(p==1){ 
           den = dnorm(x, mu, sd=sqrt(Sigma))
    } else den = dmvnorm(x, mean=mu, sigma=round(Sigma,5), log=F)
   } 
   if(distr=='MVT'){
     if(p==1){ 
           den = dt((x-mu)/sqrt(Sigma), df=nu)/sqrt(Sigma)
    } else den = dmvt(x, delta=mu, sigma=Sigma, df=nu, log=F)
   } 
   if(distr=='MCN'){
     if(p==1){
       Delta = c((x-mu)^2/c(Sigma))
       den = c(1/((2*pi)^(p/2)*sqrt(Sigma))) * (nu*rho^(p/2)*exp(-rho*Delta/2) + (1-nu)*exp(-Delta/2))
    } else{ 
       Delta = c(t(x-mu)%*%solve(Sigma)%*%(x-mu))
       den = 1/((2*pi)^(p/2)*sqrt(det(Sigma))) * (nu*rho^(p/2)*exp(-rho*Delta/2) + (1-nu)*exp(-Delta/2)) 
   }}
   return(den)
}

# CDF of multivariate normal-independent distribution
pmni = function(lower=rep(-Inf, length(mu)), upper=rep(Inf, length(mu)), mu=rep(0,length(mu)), Sigma=diag(length(mu)), nu=NULL, rho=NULL, distr=c('MVN','MVT','MCN'))
{
   require(mvtnorm)
   GB = GenzBretz(maxpts = 5e4, abseps = 1e-9, releps = 0)
   distr = distr[1]
   p = length(mu)
   if(distr=='MVN'){
    if(p==1){
      Sd = sqrt(Sigma)
      cdf = pnorm(upper, mean=mu, sd=Sd) - pnorm(lower, mean=mu, sd=Sd)
    } else cdf = pmvnorm(lower=lower, upper=upper, mean=mu, sigma=round(Sigma, 4), algorithm = GB)[1]
   }
   if(distr=='MVT'){
    if(p==1){
      Sd = sqrt(Sigma)
      cdf = pt((upper-mu)/Sd, df=nu) - pt((lower-mu)/Sd, df=nu)
    } else cdf = pmvt(lower=lower-mu, upper=upper-mu, df=nu, sigma=Sigma)[1]
   }
   if(distr=='MCN'){
      if(p==1){
       Sd = sqrt(Sigma)
       cdf = (nu * pnorm(upper, mean=mu, sd=Sd/sqrt(rho)) + (1-nu) * pnorm(upper, mean=mu, sd=Sd)) - (nu * pnorm(lower, mean=mu, sd=Sd/sqrt(rho)) + (1-nu)*pnorm(lower, mean=mu, sd=Sd)) 
      }
    cdf = nu * pmvnorm(lower=lower, upper=upper, mean=mu, sigma=round(Sigma/rho, 4), algorithm = GB)[1] + (1-nu) * pmvnorm(lower=lower, upper=upper, mean=mu, sigma=round(Sigma, 4), algorithm = GB)[1]
#     if(p==1) cdf1 = integrate(dmni, lower=lower, upper=upper, mu=mu, Sigma=Sigma, nu=nu, rho=rho, distr='MCN')$value
#     cdf1 = cubintegrate(dmni, lower=lower, upper=upper, mu=mu, Sigma=Sigma, nu=nu, rho=rho, distr='MCN')$integral
##    cdf = adaptIntegrate(dmni, lowerLimit=lower, upperLimit=upper, mu=mu, Sigma=Sigma, nu=nu, rho=rho, distr='MCN')$integral
#   cat('cdf=', cdf, '; cdf1=', cdf1, '\n')
   }
   return(cdf)  
}

# PDF of truncated multivariate normal-independent distribution
dtmni = function(x, mu, Sigma, nu=NULL, rho=NULL, distr=c('MVN','MVT','MCN'), a.low=rep(-Inf, length(mu)), a.upp=rep(Inf, length(mu)))
{
   distr = distr[1]
   p = length(mu)
   if((sum(x>=a.low)+sum(x<=a.upp))!=(2*p)){
    Tden=0
   } else {
   if(distr=='MVN'){
       cdf.inv = 1/pmni(lower=a.low, upper=a.upp, mu, Sigma, distr='MVN')
       Tden = cdf.inv * dmni(x, mu, Sigma, nu, distr='MVN') 
   }
   if(distr=='MVT'){
       cdf.inv = 1/pmni(lower=a.low, upper=a.upp, mu, Sigma, nu, distr='MVT')
       Tden = cdf.inv * dmni(x, mu, Sigma, nu, distr='MVT') 
   }
   if(distr=='MCN'){
     cdf.inv = 1/pmni(lower=a.low, upper=a.upp, mu, Sigma, nu=nu, rho=rho, distr='MCN')
     Tden = cdf.inv * dmni(x, mu, Sigma, nu=nu, rho=rho, distr='MCN') 
   }}
   return(Tden)
}

ptmni = function(lower=rep(-Inf, length(mu)), upper=rep(Inf, length(mu)), mu=rep(0,length(mu)), Sigma=diag(length(mu)), nu=NULL, rho=NULL, distr=c('MVN','MVT','MCN'), a.low=rep(-Inf, length(mu)), a.upp=rep(Inf, length(mu)))
{
   require(mvtnorm)
   GB = GenzBretz(maxpts = 5e4, abseps = 1e-9, releps = 0)
   distr = distr[1]
   p = length(mu)
   for(i in 1: p){
     if(lower[i]<a.low[i]) lower[i]=a.low[i]
     if(upper[i]>a.upp[i]) upper[i]=a.upp[i]
   }
   if(distr=='MVN'){
       cdf.inv = 1/pmni(lower=a.low, upper=a.upp, mu, Sigma, distr='MVN')
       if(p==1) Tcdf = cdf.inv * (pnorm(upper, mean=mu, sd=Sd) - pnorm(lower, mean=mu, sd=Sd))
       else Tcdf = cdf.inv * pmvnorm(lower=lower, upper=upper, mean=mu, sigma=round(Sigma, 4), algorithm = GB)[1]
   }
   if(distr=='MVT'){
       cdf.inv = 1/pmni(lower=a.low, upper=a.upp, mu, Sigma, nu, distr='MVT')
       if(p==1) Tcdf = cdf.inv * (pt((upper-mu)/sqrt(Sigma), df=nu) - pt((lower-mu)/sqrt(Sigma), df=nu))
       else Tcdf = cdf.inv * pmvt(lower=lower-mu, upper=upper-mu, df=nu, sigma=Sigma)[1]
   }
   if(distr=='MCN'){
     cdf.inv = 1/pmni(lower=a.low, upper=a.upp, mu, Sigma, nu=nu, rho=rho, distr='MCN')
     if(p==1){
       Sd = sqrt(Sigma)
       cdf = (nu * pnorm(upper, mean=mu, sd=Sd/sqrt(rho)) + (1-nu) * pnorm(upper, mean=mu, sd=Sd)) - (nu * pnorm(lower, mean=mu, sd=Sd/sqrt(rho)) + (1-nu)*pnorm(lower, mean=mu, sd=Sd)) 
     }
     cdf = nu * pmvnorm(lower=lower, upper=upper, mean=mu, sigma=round(Sigma/rho, 4), algorithm = GB)[1] + (1-nu) * pmvnorm(lower=lower, upper=upper, mean=mu, sigma=round(Sigma, 4), algorithm = GB)[1]
     Tcdf = cdf.inv * cdf
#     Tcdf1 = cdf.inv * cubintegrate(dmni, lower=lower, upper=upper, mu=mu, Sigma=Sigma, nu=nu, rho=rho, distr='MCN')$integral
#    cat('Tcdf=', Tcdf, '; Tcdf1=', Tcdf1, '\n')
   }
   return(Tcdf)
}

TMNI.moment = function(mu, Sigma, nu=NULL, rho=NULL, distr=c('MVN','MVT','MCN'), a.low=rep(-Inf, length(mu)), a.upp=rep(Inf, length(mu)))
{ 
  distr=distr[1]
  require(mvtnorm)
  GB = GenzBretz(maxpts = 5e4, abseps = 1e-9, releps = 0)
  p = length(mu)
  if(p == 1){
   Sd = sqrt(Sigma)
  } else{
   Lambda = diag(sqrt(diag(Sigma)))
   Lambda.inv = diag(1/sqrt(diag(Sigma)))
   R = Lambda.inv %*% Sigma %*% Lambda.inv
   if(det(R)<=0) stop("The R matrix must be inversible!")
  }
  if(distr=='MVN'){
    a.low = ifelse(a.low==-Inf,rep(-5e1+mu,p), a.low)
    a.upp = ifelse(a.upp==Inf,rep(5e1+mu,p), a.upp)
    if(p == 1){
     a = (a.low - mu)/Sd
     b = (a.upp - mu)/Sd
     EX = Sd*(dnorm(a)-dnorm(b))/ (pnorm(b)-pnorm(a))
     EY = mu + EX
     CovY = CovX = Sd^2*(1+(a*(dnorm(a))-b*dnorm(b))/(pnorm(b)-pnorm(a))-((dnorm(a)-dnorm(b))/(pnorm(b)-pnorm(a)))^2)
     EXX = CovX + EX^2
     EYY = CovY + EY^2
    } else{
    a = c(Lambda.inv %*% (a.low - mu))
    b = c(Lambda.inv %*% (a.upp - mu))
    al0 = pmvnorm(lower = a, upper = b, sigma = round(R,4), algorithm = GB)[1]
  ### pdf & cdf
    f1a = dnorm(a)
    f1b = dnorm(b)
    f2 = matrix(NA, p, p)
    G1a = G1b = rep(NA, p)
    G2 = matrix(NA, p, p)
    for(r in 1:p){
      temp = R[-r,r]
      S1 = R[-r,-r] - temp %*% t(R[r,-r])
      mua = temp * a[r]; low = a[-r]-mua; upp = b[-r]-mua
      G1a[r] = pmvnorm(lower = low, upper = upp, sigma = round(S1,4), algorithm = GB)[1]
      mub = temp * b[r]; low = a[-r]-mub; upp = b[-r]-mub
      G1b[r] = pmvnorm(lower = low, upper = upp, sigma = round(S1,4), algorithm = GB)[1]
    }
    qa = f1a*G1a; qb = f1b*G1b
    EX = c(R %*% (qa-qb)) / al0

    H = matrix(0,p,p)
    for(r in 1:(p-1)){
      for(s in (r+1):p){
        rs = c(r,s)
        pdf.aa = dmvnorm(t(c(a[r],a[s])),sigma=round(R[rs,rs],5), log =F)
        pdf.ab = dmvnorm(t(c(a[r],b[s])),sigma=round(R[rs,rs],5), log =F)
        pdf.ba = dmvnorm(t(c(b[r],a[s])),sigma=round(R[rs,rs],5), log =F)
        pdf.bb = dmvnorm(t(c(b[r],b[s])),sigma=round(R[rs,rs],5), log =F)
        if(p==2){cdf.aa=cdf.ab=cdf.ba=cdf.bb=1}
        if(p>2){
          tmp = R[-rs,rs]%*%solve(R[rs,rs])
          mu.aa = c(tmp%*%c(a[r],a[s]))
          mu.ab = c(tmp%*%c(a[r],b[s]))
          mu.ba = c(tmp%*%c(b[r],a[s]))
          mu.bb = c(tmp%*%c(b[r],b[s]))
          R21 = R[-rs,-rs] - R[-rs,rs]%*%solve(R[rs,rs]) %*% R[rs,-rs]
          cdf.aa = pmvnorm(lower = a[-rs], upper = b[-rs], mean=mu.aa, sigma = round(R21,5), algorithm = GB)[1]
          cdf.ab = pmvnorm(lower = a[-rs], upper = b[-rs], mean=mu.ab, sigma = round(R21,5), algorithm = GB)[1]
          cdf.ba = pmvnorm(lower = a[-rs], upper = b[-rs], mean=mu.ba, sigma = round(R21,5), algorithm = GB)[1]
          cdf.bb = pmvnorm(lower = a[-rs], upper = b[-rs], mean=mu.bb, sigma = round(R21,5), algorithm = GB)[1]
        }
        H[r,s] = H[s,r] = pdf.aa*cdf.aa - pdf.ab*cdf.ab - pdf.ba*cdf.ba + pdf.bb*cdf.bb
      }}
    D = matrix(0,p,p)
    diag(D) = a * qa - b * qb - diag(R%*%H)
    EXX = R + R %*% (H + D) %*% R / al0
  }}
  if(distr=='MVT'){
    if(nu <= 2) stop("The first moment exists only when the degree of freedom is larger than 2!")
    if(nu <= 4) stop("The theoretical second moment exists only when the degrees of freedom is larger than 4!")
    a.low = ifelse(a.low==-Inf,rep(-1e12,p), a.low)
    a.upp = ifelse(a.upp==Inf,rep(1e12,p), a.upp)
    if(p == 1){
     a = (a.low - mu)/Sd
     b = (a.upp - mu)/Sd
     ka = gamma((nu+1)/2)/((pt(b, df=nu)-pt(a, df=nu))*gamma(nu/2)*sqrt(nu*pi))
     EX = (ka*nu)/(nu-1) * ((1+a^2/nu)^(-(nu-1)/2) - (1+b^2/nu)^(-(nu-1)/2))
     tau.sq = sqrt((nu-2)/nu)
     EXX = (nu*(nu-1))/(nu-2) * ((pt(b*tau.sq, df=nu-2)-pt(a*tau.sq, df=nu-2)) / (pt(b, df=nu)-pt(a, df=nu))) - nu
     EY = mu + Sd * EX
     EYY = Sigma * EXX + mu^2 + 2 * mu * Sd * EX
     CovX = EXX - EX^2
     CovY = EYY - EY^2 
    } else{
    aL.ast = c(Lambda.inv %*% (a.low - mu))
    aU.ast = c(Lambda.inv %*% (a.upp - mu))
    tau.t = pmvt(lower = aL.ast, upper = aU.ast, sigma = R, df = round(nu), algorithm = GB)[1]
    la1 = (nu-2)/nu; la2 = (nu-4)/nu
    da = (nu-1)/(nu+aL.ast^2); db = (nu-1)/(nu+aU.ast^2)
    f1a = sqrt(la1)*dt(sqrt(la1)*aL.ast,df=nu-2)
    f1b = sqrt(la1)*dt(sqrt(la1)*aU.ast,df=nu-2)
    G1a = G1b = rep(NA, p)
    for(r in 1:p){
      temp = R[-r,r]
      S1 = R[-r,-r] - temp %*% t(R[r,-r])
      mua = temp * aL.ast[r]; low = aL.ast[-r]-mua; upp = aU.ast[-r]-mua
      G1a[r] = ifelse(p==2,pt(upp/sqrt(S1/da[r]),df=nu-1)-pt(low/sqrt(S1/da[r]),df=nu-1),pmvt(lower = low, upper = upp, sigma = S1/da[r], df = round(nu-1), algorithm = GB)[1])
      mub = temp * aU.ast[r]; low = aL.ast[-r]-mub; upp = aU.ast[-r]-mub
      G1b[r] = ifelse(p==2,pt(upp/sqrt(S1/db[r]),df=nu-1)-pt(low/sqrt(S1/db[r]),df=nu-1),pmvt(lower = low, upper = upp, sigma = S1/db[r], df = round(nu-1), algorithm = GB)[1])
    }
    qa = f1a*G1a; qb = f1b*G1b
    EX = c(R %*% (qa-qb)) / tau.t / la1
    H = matrix(0,p,p)
    for(r in 1:(p-1)){
     for(s in (r+1):p){
       rs = c(r,s)
       pdf.aa = dmvt(c(aL.ast[r],aL.ast[s]),sigma=R[rs,rs]/la2,df=nu-4, log =F)
       pdf.ab = dmvt(c(aL.ast[r],aU.ast[s]),sigma=R[rs,rs]/la2,df=nu-4, log =F)
       pdf.ba = dmvt(c(aU.ast[r],aL.ast[s]),sigma=R[rs,rs]/la2,df=nu-4, log =F)
       pdf.bb = dmvt(c(aU.ast[r],aU.ast[s]),sigma=R[rs,rs]/la2,df=nu-4, log =F)
       if(p==2){cdf.aa=cdf.ab=cdf.ba=cdf.bb=1}
       if(p>2){
         tmp = R[-rs,rs]%*%solve(R[rs,rs])
         mu.aa = c(tmp%*%c(aL.ast[r],aL.ast[s]))
         mu.ab = c(tmp%*%c(aL.ast[r],aU.ast[s]))
         mu.ba = c(tmp%*%c(aU.ast[r],aL.ast[s]))
         mu.bb = c(tmp%*%c(aU.ast[r],aU.ast[s]))
         daa = (nu-2)/(nu+(aL.ast[r]^2-2*R[r,s]*aL.ast[r]*aL.ast[s]+aL.ast[s]^2)/(1-R[r,s]^2))
         dab = (nu-2)/(nu+(aL.ast[r]^2-2*R[r,s]*aL.ast[r]*aU.ast[s]+aU.ast[s]^2)/(1-R[r,s]^2))
         dba = (nu-2)/(nu+(aU.ast[r]^2-2*R[r,s]*aU.ast[r]*aL.ast[s]+aL.ast[s]^2)/(1-R[r,s]^2))
         dbb = (nu-2)/(nu+(aU.ast[r]^2-2*R[r,s]*aU.ast[r]*aU.ast[s]+aU.ast[s]^2)/(1-R[r,s]^2))
         R21 = R[-rs,-rs] - R[-rs,rs]%*%solve(R[rs,rs]) %*% R[rs,-rs]
         cdf.aa = ifelse(p==3,pt((aU.ast[-rs]-mu.aa)/sqrt(R21/daa),df=nu-2)-pt((aL.ast[-rs]-mu.aa)/sqrt(R21/daa),df=nu-2),pmvt(lower = aL.ast[-rs]-mu.aa, upper = aU.ast[-rs]-mu.aa, sigma = R21/daa, df=round(nu-2), algorithm = GB)[1])
         cdf.ab = ifelse(p==3,pt((aU.ast[-rs]-mu.ab)/sqrt(R21/dab),df=nu-2)-pt((aL.ast[-rs]-mu.ab)/sqrt(R21/dab),df=nu-2),pmvt(lower = aL.ast[-rs]-mu.ab, upper = aU.ast[-rs]-mu.ab, sigma = R21/dab, df=round(nu-2), algorithm = GB)[1])
         cdf.ba = ifelse(p==3,pt((aU.ast[-rs]-mu.ba)/sqrt(R21/dba),df=nu-2)-pt((aL.ast[-rs]-mu.ba)/sqrt(R21/dba),df=nu-2),pmvt(lower = aL.ast[-rs]-mu.ba, upper = aU.ast[-rs]-mu.ba, sigma = R21/dba, df=round(nu-2), algorithm = GB)[1])
         cdf.bb = ifelse(p==3,pt((aU.ast[-rs]-mu.bb)/sqrt(R21/dbb),df=nu-2)-pt((aL.ast[-rs]-mu.bb)/sqrt(R21/dbb),df=nu-2),pmvt(lower = aL.ast[-rs]-mu.bb, upper = aU.ast[-rs]-mu.bb, sigma = R21/dbb, df=round(nu-2), algorithm = GB)[1])
       }
       H[r,s] = H[s,r] = pdf.aa*cdf.aa - pdf.ab*cdf.ab - pdf.ba*cdf.ba + pdf.bb*cdf.bb
     }}
     H = H / la2
     D = matrix(0,p,p)
     diag(D) = aL.ast * qa - aU.ast * qb - diag(R%*%H)
     al1 = pmvt(lower = aL.ast, upper = aU.ast, sigma = R/la1, df=round(nu-2), algorithm = GB)[1]
     EXX = (al1 * R + R %*% (H + D) %*% R) / tau.t / la1
  }}
  if(distr=='MCN'){
    if(nu <= 0 | nu>=1) stop("The mixing probabiity must be between 0 to 1!")
    a.low = ifelse(a.low==-Inf,rep(-3e1+mu,p), a.low)
    a.upp = ifelse(a.upp==Inf,rep(3e1+mu,p), a.upp)
    if(p == 1){
     a = (a.low - mu)/Sd
     b = (a.upp - mu)/Sd
     Fa = pmni(lower=-Inf, upper=a, mu=0, Sigma=1, nu=nu, rho=rho, distr='MCN')
     Fb = pmni(lower=-Inf, upper=b, mu=0, Sigma=1, nu=nu, rho=rho, distr='MCN')
     sqrho = sqrt(rho)
     EX = (nu/sqrho*(dnorm(sqrho*a)-dnorm(sqrho*b)) + (1-nu)*(dnorm(a)-dnorm(b)))/(Fb-Fa)
     EXX = 1/rho + (nu/sqrho*(a*dnorm(sqrho*a)-b*dnorm(sqrho*b)) + (1-nu)*((a*dnorm(a)-b*dnorm(b))+(1/rho-1)*(pnorm(a)-pnorm(b))))/(Fb-Fa)
     EY = mu + Sd * EX
     EYY = Sigma * EXX + mu^2 + 2 * mu * Sd * EX
     CovX = EXX - EX^2
     CovY = EYY - EY^2 
    } else{    
    aL.ast = c(Lambda.inv %*% (a.low - mu))
    aU.ast = c(Lambda.inv %*% (a.upp - mu))
    tau.c = pmni(lower = aL.ast, upper = aU.ast, mu=rep(0,p), Sigma=R, nu=nu, rho=rho, distr='MCN')
    G1a = G1b = rep(NA, p)
    irho = 1/rho
    irho2 = 1/rho^2
    for(r in 1:p){
      temp = R[-r,r]
      S1 = R[-r,-r] - temp %*% t(R[r,-r])
      sd1 = sqrt(S1)
      mua = temp * aL.ast[r]
      f1a = dnorm(aL.ast[r], mean=0, sd=sqrt(irho))/rho * pmni(lower = aL.ast[-r], upper = aU.ast[-r], mu=mua, Sigma=S1, nu=nu, rho=rho, distr='MCN')
      la1 = (dmni(aL.ast[r], mu=0, Sigma=1, nu=nu, rho=rho, distr='MCN')/rho + (1-nu)*(1-irho)*dnorm(aL.ast[r]) - dnorm(aL.ast[r],mean=0,sd=sqrt(irho))/rho)        
      G1a[r] = ifelse(p==2, f1a+la1*(pnorm(aU.ast[-r], mean=mua, sd=sd1)-pnorm(aL.ast[-r], mean=mua, sd=sd1)), f1a+la1*pmvnorm(lower=aL.ast[-r], upper=aU.ast[-r], mean=mua, sigma=round(S1,5), algorithm = GB))
      mub = temp * aU.ast[r]
      f1b = dnorm(aU.ast[r], mean=0, sd=sqrt(irho))/rho * pmni(lower = aL.ast[-r], upper = aU.ast[-r], mu=mub, Sigma=S1, nu=nu, rho=rho, distr='MCN') 
      la2 = (dmni(aU.ast[r], mu=0, Sigma=1, nu=nu, rho=rho, distr='MCN')/rho + (1-nu)*(1-irho)*dnorm(aU.ast[r]) - dnorm(aU.ast[r],mean=0,sd=sqrt(irho))/rho)        
      G1b[r] = ifelse(p==2, f1b+la2*(pnorm(aU.ast[-r], mean=mub, sd=sd1)-pnorm(aL.ast[-r], mean=mub, sd=sd1)), f1b+la2*pmvnorm(lower=aL.ast[-r], upper=aU.ast[-r], mean=mub, sigma=round(S1,5), algorithm = GB))
    }
    qa = G1a; qb = G1b
    EX = c(R %*% (qa-qb)) / tau.c

    H = matrix(0,p,p)
    for(r in 1:(p-1)){
     for(s in (r+1):p){
       rs = c(r,s)
       pdf1.aa = dmvnorm(t(c(aL.ast[r],aL.ast[s])), mean=rep(0,2), sigma=round(irho*R[rs, rs], 5))*irho2
       pdf1.ab = dmvnorm(t(c(aL.ast[r],aU.ast[s])), mean=rep(0,2), sigma=round(irho*R[rs, rs], 5))*irho2
       pdf1.ba = dmvnorm(t(c(aU.ast[r],aL.ast[s])), mean=rep(0,2), sigma=round(irho*R[rs, rs], 5))*irho2
       pdf1.bb = dmvnorm(t(c(aU.ast[r],aU.ast[s])), mean=rep(0,2), sigma=round(irho*R[rs, rs], 5))*irho2
       pdf2.aa = dmni(c(aL.ast[r],aL.ast[s]), mu=rep(0,2), Sigma=R[rs, rs], nu=nu, rho=rho, distr='MCN')*irho2 + (1-nu)*(1-irho2)*dmvnorm(t(c(aL.ast[r],aL.ast[s])), mean=rep(0,2), sigma=round(R[rs, rs], 5)) - pdf1.aa
       pdf2.ab = dmni(c(aL.ast[r],aU.ast[s]), mu=rep(0,2), Sigma=R[rs, rs], nu=nu, rho=rho, distr='MCN')*irho2 + (1-nu)*(1-irho2)*dmvnorm(t(c(aL.ast[r],aU.ast[s])), mean=rep(0,2), sigma=round(R[rs, rs], 5)) - pdf1.ab
       pdf2.ba = dmni(c(aU.ast[r],aL.ast[s]), mu=rep(0,2), Sigma=R[rs, rs], nu=nu, rho=rho, distr='MCN')*irho2 + (1-nu)*(1-irho2)*dmvnorm(t(c(aU.ast[r],aL.ast[s])), mean=rep(0,2), sigma=round(R[rs, rs], 5)) - pdf1.ba
       pdf2.bb = dmni(c(aU.ast[r],aU.ast[s]), mu=rep(0,2), Sigma=R[rs, rs], nu=nu, rho=rho, distr='MCN')*irho2 + (1-nu)*(1-irho2)*dmvnorm(t(c(aU.ast[r],aU.ast[s])), mean=rep(0,2), sigma=round(R[rs, rs], 5)) - pdf1.bb 
       if(p==2){
        cdf1.aa = cdf1.ab = cdf1.ba= cdf1.bb = cdf2.aa = cdf2.ab = cdf2.ba= cdf2.bb = 1
       }
       if(p>2){
         tmp = R[-rs,rs]%*%solve(R[rs,rs])
         mu.aa = c(tmp%*%c(aL.ast[r],aL.ast[s]))
         mu.ab = c(tmp%*%c(aL.ast[r],aU.ast[s]))
         mu.ba = c(tmp%*%c(aU.ast[r],aL.ast[s]))
         mu.bb = c(tmp%*%c(aU.ast[r],aU.ast[s]))
         R21 = R[-rs,-rs] - R[-rs,rs] %*% solve(R[rs,rs]) %*% R[rs,-rs]
         if(p==3) sd21 = sqrt(c(R21))
         cdf1.aa = pmni(lower = aL.ast[-rs], upper = aU.ast[-rs], mu=mu.aa, Sigma=R21, nu=nu, rho=rho, distr='MCN')
         cdf1.ab = pmni(lower = aL.ast[-rs], upper = aU.ast[-rs], mu=mu.ab, Sigma=R21, nu=nu, rho=rho, distr='MCN')
         cdf1.ba = pmni(lower = aL.ast[-rs], upper = aU.ast[-rs], mu=mu.ba, Sigma=R21, nu=nu, rho=rho, distr='MCN')
         cdf1.bb = pmni(lower = aL.ast[-rs], upper = aU.ast[-rs], mu=mu.bb, Sigma=R21, nu=nu, rho=rho, distr='MCN')
         cdf2.aa = ifelse(p==3, (pnorm(aU.ast[-rs], mean=mu.aa, sd=sd21) - pnorm(aL.ast[-rs], mean=mu.aa, sd=sd21)), pmvnorm(lower = aL.ast[-rs], upper = aU.ast[-rs], mean=mu.aa, sigma=round(R21,5)))
         cdf2.ab = ifelse(p==3, (pnorm(aU.ast[-rs], mean=mu.ab, sd=sd21) - pnorm(aL.ast[-rs], mean=mu.ab, sd=sd21)), pmvnorm(lower = aL.ast[-rs], upper = aU.ast[-rs], mean=mu.ab, sigma=round(R21,5)))
         cdf2.ba = ifelse(p==3, (pnorm(aU.ast[-rs], mean=mu.ba, sd=sd21) - pnorm(aL.ast[-rs], mean=mu.ba, sd=sd21)), pmvnorm(lower = aL.ast[-rs], upper = aU.ast[-rs], mean=mu.ba, sigma=round(R21,5)))
         cdf2.bb = ifelse(p==3, (pnorm(aU.ast[-rs], mean=mu.bb, sd=sd21) - pnorm(aL.ast[-rs], mean=mu.bb, sd=sd21)), pmvnorm(lower = aL.ast[-rs], upper = aU.ast[-rs], mean=mu.bb, sigma=round(R21,5)))
       }
       H[r,s] = H[s,r] = ((pdf1.aa*cdf1.aa+pdf2.aa*cdf2.aa)- (pdf1.ab*cdf1.ab+pdf2.ab*cdf2.ab) - (pdf1.ba*cdf1.ba+pdf2.ba*cdf2.ba) + (pdf1.bb*cdf1.bb+pdf2.bb*cdf2.bb))
     }}
     K = matrix(0,p,p)
     diag(K) = aL.ast * qa - aU.ast * qb - diag(R%*%H)
     al1 = tau.c/rho + (1-nu)*(1-irho)*pmvnorm(lower = aL.ast, upper = aU.ast, mean=rep(0,p), sigma=round(R,5))
     EXX = (al1*R + R %*% (H + K) %*% R) / tau.c
  }}
  if(p >= 2){
   EY = c(mu + Lambda %*% EX)
   EYY = mu%*%t(mu) + Lambda%*%EX%*%t(mu) + mu%*%t(EX)%*%Lambda + Lambda%*%EXX%*%Lambda
   CovY = EYY-(EY)%*%t(EY)
  }
  return(list(EY=EY, EYY=EYY, CovY=CovY, EX=EX, EXX=EXX))
}

TMCN.moment = function(mu, Sigma=diag(length(mu)), nu=NULL, rho=NULL, 
                      a.low=rep(-Inf, length(mu)), a.upp=rep(Inf, length(mu)))
{
  GB = GenzBretz(maxpts = 5e4, abseps = 1e-9, releps = 0)
  Phi1 = pmvnorm(lower=a.low, upper=a.upp, mean=mu, sigma=Sigma/rho, algorithm = GB)[1]
  Phi2 = pmvnorm(lower=a.low, upper=a.upp, mean=mu, sigma=Sigma, algorithm = GB)[1]
  MCNcdf = nu * Phi1 + (1-nu) * Phi2
  ans1 = TMNI.moment(mu, Sigma/rho, distr='MVN', a.low=a.low, a.upp=a.upp)
  ans2 = TMNI.moment(mu, Sigma, distr='MVN', a.low=a.low, a.upp=a.upp)
  EU1 = ans1$EY
  EUU1 = ans1$EYY
  EU2 = ans2$EY
  EUU2 = ans2$EYY
  EY = (nu * Phi1 * EU1 + (1-nu) * Phi2 * EU2) / MCNcdf
  EYY = (nu * Phi1 * EUU1 + (1-nu) * Phi2 * EUU2) / MCNcdf
  return(list(EY=EY, EYY=EYY))
}
 