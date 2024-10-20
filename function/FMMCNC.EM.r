library(combinat)

# Basic functions:
tr = function(M)  sum(diag(M))
vech.posi=function(dim) cbind(rep(1:dim, 1:dim), unlist(mapply(':', 1, 1:dim)))
as.vech = function(M)
{
    Dim = dim(M)[1]
    posi = vech.posi(Dim)
    temp = paste('M', posi[,1], posi[,2], sep = '')
    vech.M = matrix(M[posi], nrow=1)
    colnames(vech.M) = temp
    return(vech.M)
}

mse = function(v1, v2) mean((v1 - v2)^2)

#Generate FM-MCN data
rfmmcn = function(n, para, distr=c('MCN','MVN'))
{
  distr = distr[1]
  w = para$w
  mu = para$mu 
  Sigma = para$Sigma 
  nu = para$nu 
  rho = para$rho
  g = length(w); p = nrow(mu)
  N = rmultinom(n=1, size=n, prob=w)
  cls = rep((1:g), N)
  Tau = c()
  Y = NULL
  for(i in 1: g){
   if(distr == 'MVN'){
    tau = rep(1, N[i])
   } 
   if(distr == 'MCN'){
     a = runif(N[i])
     tau = ifelse(a<=nu[i], rho[i], 1)
     tau.sq = sqrt(tau)
   }
   Y = rbind(Y, matrix(rep(mu[,i], N[i]), nrow=N[i], byrow=T) + 1/sqrt(tau) * rmvnorm(N[i], mean=rep(0,p), sigma=Sigma[,,i]))
   Tau = c(Tau, tau)
  } 
  return(list(Y = Y, cls = cls, tau = Tau))
}

#Create censored data
rcen = function(Y, cen.rate)
{
  n = nrow(Y); p = ncol(Y)
  CEN = matrix(0, n, p)
  ubd = apply(Y, 2, quantile, prob=cen.rate)
  for(k in 1: p) Y[,k] = ifelse(Y[,k] <= ubd[k], ubd[k], Y[,k])
  for(k in 1: p) CEN[,k] = ifelse(Y[,k] == ubd[k], 1, 0)
  return(list(Y=Y, CEN=CEN))
}

# Density for FM-MCN distribution
dfmmcn = function(x, g=ncol(mu), w=rep(1/g, g), mu=matrix(rep(0,nrow(mu)*g),ncol=g),
                  Sigma=array(diag(nrow(mu)), dim=c(nrow(mu),nrow(mu),g)),
                  nu=NULL, rho=NULL, distr=c('MVN','MVT','MCN'))
{
   distr = distr[1]
   den = 0
   for(i in 1:g) den = den + w[i] * dmni(x, mu=mu[,i], Sigma=Sigma[,,i], nu=nu[i], rho=rho[i], distr=distr)
   return(den)
}

# ECM Algorithm for FMMCN-CR and FMMN-CR Models
FMMCNC.EM = function(Yc, cen, g, distr=c('MVN','MCN'), true.clus=NULL, tol=1e-6, max.iter = 1000, per=10, mu.true=NULL, SE=FALSE)
{
  require(mvtnorm)
  GB = GenzBretz(maxpts = 5e4, abseps = 1e-9, releps = 0)
  begin = proc.time()[1]
  distr = distr[1]
 # initial values:
  n = nrow(Yc)
  p = ncol(Yc)
  vechS = vech.posi(p)

# initial values
  if(g == 1){ 
    clus = rep(1, n)
  } else{
    if(length(true.clus) == 0 | length(unique(true.clus)) != g){
     kmY = kmeans(Yc, g, nstart = 25)
     clus = kmY$cluster
    } else{ clus = true.clus}
  }
  w = rep(1/g, g)
  mu = matrix(NA, p, g)
  Sigma = Sig.inv = array(NA, dim = c(p,p,g))
  for(i in 1:g)
  {
    w[i] = nrow(Yc[clus == i,])/n
    mu[,i] = colMeans(Yc[clus == i,])
    Sigma[,,i] = cov(Yc[clus == i,])
  }
  if(distr == 'MVN'){ rho = rep(1, g)}
  if(distr == 'MCN'){ 
    nu = rep(0.25, g) 
    rho = rep(0.2, g)
  }
# fix label switch
  if(length(mu.true) != 0){
     allper = permn(g) 
     pk = factorial(g)
     mu.norm = length(pk)
     for(k in 1: pk) mu.norm[k] = sum((mu[,allper[[k]]] - mu.true)^2)
     cho = order(mu.norm)[1]
     w = w[allper[[cho]]]
     mu = mu[,allper[[cho]]]
     Sigma = Sigma[,,allper[[cho]]]
  }   
  for(i in 1: g) Sig.inv[,,i] = solve(Sigma[,,i])

  Ip = diag(p)
  po = p - rowSums(cen)
  cen.subj = which(rowSums(cen) != 0)
  Nc = length(cen.subj)
  obs.subj = which(rowSums(cen) == 0)
  No = n - Nc
  ind.cen = colSums(t(cen) * 2 ^ (1:p - 1))
  num.class.cen = length(unique(ind.cen))
  row.posi = O.list = C.list = as.list(numeric(num.class.cen))
  uni.ind = unique(ind.cen)
  
  for(j in 1: num.class.cen){
    row.posi[[j]] = which(ind.cen == uni.ind[j])
    O.list[[j]] = matrix(Ip[!cen[row.posi[[j]][1],],], ncol = p)
    C.list[[j]] = matrix(Ip[which(cen[row.posi[[j]][1], ]==1), ], ncol=p)
  }
  
  wden = matrix(NA, n, g)
  det.o = matrix(NA, n, g)
  delta.o = matrix(NA, n, g)
  pdf = cdf = w.cdf.pdf = matrix(1, n, g)
  cent = array(NA, dim = c(p, n, g))

# old observed-data log-likelihood
  for(i in 1: g){
   cent[,,i] = t(Yc) - mu[,i]
  
   for(j in 1:num.class.cen){
    Cj = C.list[[j]]
    ind = row.posi[[j]]
    no.ind = length(ind)
    pjc = nrow(Cj)
    if(pjc == p){
      det.o[ind, i] = 1
      delta.o[ind, i] = 0
      mu.co = matrix(rep(mu[, i], no.ind), ncol=no.ind)
      Scc.o = Sigma[,,i]
    } else{
      O = O.list[[j]]
      OSO = O %*% Sigma[,,i] %*% t(O)
      cent.ind = cent[, ind, i]
      det.o[ind, i] = det(OSO)
      Soo = t(O) %*% solve(OSO) %*% O
      delta.o[ind, i] = colSums(cent.ind * (Soo %*% cent.ind))
      mu.co = Cj %*% matrix(rep(mu[, i], no.ind), ncol=no.ind) + Cj %*% Sigma[,,i] %*% Soo %*% cent.ind
      Scc.o = Cj %*% (Ip - Sigma[,,i] %*% Soo) %*% Sigma[,,i] %*% t(Cj)
    }
    if(distr=='MVN'){
      if(pjc == 1){
        for(s in 1: no.ind) cdf[ind[s], i] = pnorm(Cj%*%Yc[ind[s], ], mean=mu.co[,s], sd=sqrt(Scc.o))
      }
      if(pjc >= 2){
        for(s in 1: no.ind) cdf[ind[s], i] = ptmvnorm(lowerx=rep(-Inf,pjc), upperx=c(Cj%*%Yc[ind[s], ]), mean=c(mu.co[,s]), sigma=round(Scc.o, 6))
    }}
    if(distr=='MCN'){
      if(pjc == 0){
        for(s in 1: no.ind) w.cdf.pdf[ind[s], i] = nu[i] * dmvnorm(t(Yc[ind[s], ]), mean=mu[, i], sigma=round(Sigma[,,i]/rho[i],6)) + (1-nu[i]) * dmvnorm(t(Yc[ind[s], ]), mean=mu[, i], sigma=round(Sigma[,,i], 6)) 
      }
      if(pjc == 1){
        for(s in 1: no.ind){
        fyo = dmni(c(O%*%Yc[ind[s], ]), mu=c(O%*%mu[, i]), Sigma=round(OSO, 6), nu=nu[i], rho=rho[i], distr='MCN') 
        nu2.1 = nu[i] * dmvnorm(t(c(O%*%Yc[ind[s], ])), mean=c(O%*%mu[, i]), sigma=round(OSO/rho[i], 6)) / fyo
        w.cdf.pdf[ind[s], i] = (nu2.1 * pnorm(Cj%*%Yc[ind[s], ], mean=mu.co[,s], sd=sqrt(Scc.o/rho[i])) + (1-nu2.1) * pnorm(Cj%*%Yc[ind[s], ], mean=mu.co[,s], sd=sqrt(Scc.o))) * fyo
      }}                                   
      if(pjc >= 2 & pjc < (p-1)){
        for(s in 1: no.ind){ 
        fyo = dmni(c(O%*%Yc[ind[s], ]), mu=c(O%*%mu[, i]), Sigma=round(OSO,6), nu=nu[i], rho=rho[i], distr='MCN')
        nu2.1 = nu[i] * dmvnorm(t(c(O%*%Yc[ind[s], ])), mean=c(O%*%mu[, i]), sigma=round(OSO/rho[i], 6)) / fyo
        w.cdf.pdf[ind[s], i] = (nu2.1 * ptmvnorm(lowerx=rep(-Inf,pjc), upperx=c(Cj%*%Yc[ind[s], ]), mean=c(mu.co[,s]), sigma=round(Scc.o/rho[i],6)) + (1-nu2.1) * ptmvnorm(lowerx=rep(-Inf,pjc), upperx=c(Cj%*%Yc[ind[s], ]), mean=c(mu.co[,s]), sigma=round(Scc.o, 6))) * fyo  
      }} 
      if(pjc == (p-1)){
        for(s in 1: no.ind){ 
        fyo = dmni(c(O%*%Yc[ind[s], ]), mu=c(O%*%mu[, i]), Sigma=round(OSO, 6), nu=nu[i], rho=rho[i], distr='MCN')
        nu2.1 = nu[i] * dnorm(c(O%*%Yc[ind[s], ]), mean=c(O%*%mu[, i]), sd=sqrt(OSO/rho[i])) / fyo
        w.cdf.pdf[ind[s], i] = (nu2.1 * ptmvnorm(lowerx=rep(-Inf,pjc), upperx=c(Cj%*%Yc[ind[s], ]), mean=c(mu.co[,s]), sigma=round(Scc.o/rho[i],6)) + (1-nu2.1) * ptmvnorm(lowerx=rep(-Inf,pjc), upperx=c(Cj%*%Yc[ind[s], ]), mean=c(mu.co[,s]), sigma=round(Scc.o, 6))) * fyo  
      }} 
      if(pjc == p){
        for(s in 1: no.ind) w.cdf.pdf[ind[s], i] = nu[i] * ptmvnorm(lowerx=rep(-Inf,pjc), upperx=c(Cj%*%Yc[ind[s], ]), mean=c(mu.co[,s]), sigma=round(Scc.o/rho[i],6)) + (1-nu[i]) * ptmvnorm(lowerx=rep(-Inf,pjc), upperx=c(Cj%*%Yc[ind[s], ]), mean=c(mu.co[,s]), sigma=round(Scc.o,6))
    }}
    }
    if(distr=='MVN') wden[, i] = log(w[i]) -.5*po*log(2*pi) -.5*log(det.o[,i]) -.5*delta.o[,i] + log(cdf[,i])
    if(distr=='MCN') wden[, i] = log(w[i]) + log(w.cdf.pdf[, i])
    }
    max.wden = apply(wden, 1, max)
    wden = exp(wden - max.wden)
    indv.den = rowSums(wden)
    log.indv.den = log(indv.den) + max.wden
    iter.lnL = loglik.old = sum(log.indv.den)    

    iter = 0
    old.par = NULL
    if(distr == 'MVN'){
      for(i in 1: g) old.par = c(old.par, c(w[i], mu[,i], Sigma[,,i][vechS]))
    }
    if(distr == 'MCN'){
      for(i in 1: g) old.par = c(old.par, c(w[i], mu[,i], Sigma[,,i][vechS], nu[i], rho[i]))
    }
    iter.EST = c(iter, old.par)
    cat(paste(rep("=", 50), collapse = ""), "\n")
    cat("EM is running for the ", g, "-component FM censored model with ", distr, "distribution: censoring proportion ", mean(cen)*100, '%.', "\n")
    cat("Initial log-likelihood = ", loglik.old, "\n")

    xi = matrix(1, n, g)
    xiyhat = yhat = array(NA, dim=c(n, p, g))
    xiy2hat = y2hat = array(NA, dim=c(p, p, n, g))
    repeat
    {
      iter = iter + 1 
      a.hat = wden / indv.den
      Ni = colSums(a.hat)
      w = Ni / n

#E-step
      for(i in 1: g){
       for(j in 1:num.class.cen){
        Cj = C.list[[j]]
        ind = row.posi[[j]]               
        no.ind = length(ind)
        pjc = nrow(Cj)     #p-pjo
        
        if(pjc == 0){ 
          cent.ind = cent[, ind, i]
          d = diag(t(cent.ind) %*% Sig.inv[,,i] %*% cent.ind)
          if(distr=='MCN'){
            w.rho = w1 = numeric(no.ind) 
            for(s in 1: no.ind){
              dcn = dmni(c(Yc[ind[s], ]), mu=mu[, i], Sigma=round(Sigma[,,i], 6), nu=nu[i], rho=rho[i], distr='MCN')
              w.rho[s] = nu[i] * dmvnorm(t(c(Yc[ind[s], ])), mean=mu[, i], sigma=round(Sigma[,,i]/rho[i], 6)) / dcn 
              w1[s] = (1-nu[i]) * dmvnorm(t(c(Yc[ind[s], ])), mean=mu[, i], sigma=round(Sigma[,,i],6)) / dcn
            }
            xi[ind, i] = w.rho 
          }
          yhat[ind, , i] = Yc[ind, ]
          xiyhat[ind, ,i] = xi[ind, i] * Yc[ind, ] 
          for(s in 1: no.ind){
           y2hat[,,ind[s], i] = Yc[ind[s], ] %*% t(Yc[ind[s], ])
           xiy2hat[,,ind[s], i] = xi[ind[s], i]*Yc[ind[s], ] %*% t(Yc[ind[s], ])
          }} else{

          yc.hat = xi.yc.hat = array(NA, dim = c(pjc, no.ind, g))
          yc2.hat = xi.yc2.hat = array(NA, dim=c(pjc, pjc, no.ind, g))          
          if(pjc == p){
            mu.co = matrix(rep(mu[, i], no.ind), ncol=no.ind)
            Scc.o = Sigma[,,i]

            if(distr=='MVN'){
              for(s in 1: no.ind){
               EX = TMNI.moment(mu=mu.co[,s], Sigma=round(Scc.o, 6), distr='MVN', a.low=rep(-Inf, pjc), a.upp=c(Cj%*%Yc[ind[s], ])) 
               yc.hat[, s, i] = EX$EY
               xiy2hat[,, ind[s], i] = y2hat[,, ind[s], i] = EX$EYY  
             }
             xiyhat[ind, , i] = yhat[ind, , i] = t(t(Cj) %*% yc.hat[,,i])
            }
            if(distr=='MCN'){
              for(s in 1: no.ind){
               yjc = c(Cj%*%Yc[ind[s], ])
               pcn = 1 / pmni(lower=rep(-Inf, pjc), upper=yjc, mu=mu.co[,s], Sigma=round(Sigma[,,i], 6), nu=nu[i], rho=rho[i], distr='MCN')
               pn.rho = pmvnorm(lower = rep(-Inf, pjc), upper = yjc, mean=mu.co[,s], sigma = round(Sigma[,,i]/rho[i], 6), algorithm = GB)[1]
               xi[ind[s], i] = nu[i] * pn.rho * pcn 
               EX.rho = TMNI.moment(mu=mu.co[,s], Sigma=round(Scc.o, 6), nu=nu[i], rho=rho[i], distr='MCN', a.low=rep(-Inf, pjc), a.upp=yjc)
               EX1 = TMNI.moment(mu=mu.co[,s], Sigma=round(Scc.o/rho[i], 6), distr='MVN', a.low=rep(-Inf, pjc), a.upp=yjc) 
               xi.yc.hat[, s, i] = xi[ind[s], i] * EX1$EY
               xiy2hat[,, ind[s], i] = xi[ind[s], i] * EX1$EYY 
               yc.hat[, s, i] = EX.rho$EY
               y2hat[,, ind[s], i] = EX.rho$EYY
              }
              xiyhat[ind, , i] = t(t(Cj) %*% xi.yc.hat[,, i]) 
              yhat[ind, , i] = t(t(Cj) %*% yc.hat[,, i])
            }  
            }else{
            O = O.list[[j]]
            OSO = O %*% Sigma[,,i] %*% t(O)
            Soo = t(O) %*% solve(OSO) %*% O
            cent.ind = cent[, ind, i]
            mu.co = Cj %*% matrix(rep(mu[, i], no.ind), ncol=no.ind) + Cj %*% Sigma[,,i] %*% Soo %*% cent.ind
            Scc.o = Cj %*% (Ip - Sigma[,,i] %*% Soo) %*% Sigma[,,i] %*% t(Cj)
            d = diag(t(cent.ind) %*% Soo %*% cent.ind)
            if(distr=='MVN'){
              for(s in 1: no.ind){
               EX = TMNI.moment(mu=mu.co[,s], Sigma=round(Scc.o, 6), distr='MVN', a.low=rep(-Inf, pjc), a.upp=c(Cj%*%Yc[ind[s], ])) 
               yc.hat[, s, i] = EX$EY
               yc2.hat[,, s, i] = EX$EYY
               yjo = O %*% Yc[ind[s],]
               xiy2hat[,, ind[s], i] = y2hat[,, ind[s], i] = xi[ind[s], i]*(t(O)%*% yjo %*% t(yjo) %*% O + t(O) %*% yjo %*% t(yc.hat[,s,i]) %*% Cj + t(Cj) %*% yc.hat[,s,i] %*% t(yjo) %*% O + t(Cj) %*% yc2.hat[,,s,i] %*% Cj)
              }
            xiyhat[ind, , i] = yhat[ind, , i] = xi[ind, i]*(t(t(O) %*% O %*% t(matrix(Yc[ind, ], ncol=p)) + t(Cj) %*% yc.hat[,,i])) 
            } 
            if(distr=='MCN'){
              w.rho = w1 = pn.rho = pn1 = w.cdf = numeric(no.ind)
              mu.o = O %*% mu[, i] 
              pjo = nrow(O)
              for(s in 1: no.ind){
                yjo = c(O %*% Yc[ind[s],])
                yjc = c(Cj %*% Yc[ind[s], ])
                dcn = dmni(yjo, mu=mu.o, Sigma=round(OSO, 6), nu=nu[i], rho=rho[i], distr='MCN')
                if(pjo == 1){
                  w.rho[s] = nu[i] * dnorm(yjo, mean=c(mu.o), sd=sqrt(OSO/rho[i])) / dcn
                  w1[s] = (1-nu[i]) * dnorm(yjo, mean=c(mu.o), sd=sqrt(OSO)) / dcn
                } else{
                  w.rho[s] = nu[i] * dmvnorm(t(yjo), mean=mu.o, sigma=round(OSO/rho[i],6)) / dcn
                  w1[s] = (1-nu[i]) * dmvnorm(t(yjo), mean=mu.o, sigma=OSO) / dcn
                }
                if(pjc == 1){
                  pn.rho[s] = pnorm(yjc, mean=c(mu.co[,s]), sd=sqrt(Scc.o/rho[i]))
                  pn1[s] = pnorm(yjc, mean=c(mu.co[,s]), sd=sqrt(Scc.o))
                } else{
                  pn.rho[s] = pmvnorm(lower = rep(-Inf, pjc), upper = yjc, mean=mu.co[,s], sigma = round(Scc.o/rho[i],6), algorithm = GB)[1]
                  pn1[s] = pmvnorm(lower = rep(-Inf, pjc), upper = yjc, mean=mu.co[,s], sigma = round(Scc.o, 6), algorithm = GB)[1]
                } 
                w.cdf[s] = w.rho[s]*pn.rho[s] + w1[s]*pn1[s] 
                EX.rho = TMNI.moment(mu=mu.co[,s], Sigma=round(Scc.o/rho[i], 6), distr='MVN', a.low=rep(-Inf, pjc), a.upp=yjc)
                EX1 = TMNI.moment(mu=mu.co[,s], Sigma=round(Scc.o, 6), distr='MVN', a.low=rep(-Inf, pjc), a.upp=yjc) 
                yc.hat[,s,i] = (w.rho[s]*pn.rho[s]*EX.rho$EY + w1[s]*pn1[s]*EX1$EY) / w.cdf[s]
                yc2.hat[,,s,i] = (w.rho[s]*pn.rho[s]*EX.rho$EYY + w1[s]*pn1[s]*EX1$EYY) / w.cdf[s]   
                xi[ind[s], i] = (w.rho[s]*pn.rho[s]) / w.cdf[s] 
                xi.yc.hat[, s, i] = xi[ind[s], i] * EX.rho$EY  
                xi.yc2.hat[,, s, i] = xi[ind[s], i] * EX.rho$EYY
                y2hat[,, ind[s], i] = (t(O)%*% yjo %*% t(yjo) %*% O) + (t(O) %*% yjo %*% t(yc.hat[,s, i]) %*% Cj) + (t(Cj) %*% yc.hat[, s, i] %*% t(yjo) %*% O) + (t(Cj) %*% yc2.hat[,, s, i] %*% Cj) 
                xiy2hat[,, ind[s], i] = xi[ind[s], i]*(t(O)%*% yjo %*% t(yjo) %*% O) + (t(O) %*% yjo %*% t(xi.yc.hat[,s,i]) %*% Cj) + (t(Cj) %*% xi.yc.hat[,s,i] %*% t(yjo) %*% O) + (t(Cj) %*% xi.yc2.hat[,,s,i] %*% Cj) 
              }
              yhat[ind, , i] = t(t(O) %*% O %*% t(matrix(Yc[ind, ], ncol=p))) + t(t(Cj) %*% yc.hat[,,i])
              xiyhat[ind, , i] = xi[ind, i]*t(t(O) %*% O %*% t(matrix(Yc[ind, ], ncol=p))) + t(t(Cj) %*% xi.yc.hat[,,i]) 
            } 
          }
      }}}    

# M-step
     for(i in 1: g){
      mu[,i] = colSums(a.hat[,i]*((rho[i]-1)*xiyhat[,,i] + yhat[,,i])) / sum(a.hat[,i]*((rho[i]-1)*xi[,i]+1))
      E1 = E2 = array(NA, dim=c(p, p, n))
      for(j in 1: n){
       E1[,,j] = a.hat[j, i] * (xiy2hat[,,j,i] - xiyhat[j,,i]%*%t(mu[,i]) - mu[,i]%*%t(xiyhat[j,,i]) + xi[j,i]*mu[,i]%*%t(mu[,i]))
       E2[,,j] = a.hat[j, i] * (y2hat[,,j,i] - yhat[j,,i]%*%t(mu[,i]) - mu[,i]%*%t(yhat[j,,i]) + mu[,i]%*%t(mu[,i]))
      }
      Sigma[,,i] = ((rho[i]-1)*apply(E1, 1:2, sum) + apply(E2, 1:2, sum)) / Ni[i]
      Sigma[,,i] = (Sigma[,,i] + t(Sigma[,,i]))/2
      Sig.inv[,,i] = solve(Sigma[,,i])
      if(distr == 'MCN'){
        nu[i] = sum(a.hat[,i] * xi[,i]) / Ni[i]
        deo = 0
        for(j in 1: n) deo = deo + tr((Sig.inv[,,i] %*% E1[,,j]))
        rho[i] = min(1, ((p*sum(a.hat[,i] * xi[,i])) / deo))
     }}

# new observed-data log-likelihood
  for(i in 1: g){
   cent[,,i] = t(Yc) - mu[,i]
  
   for(j in 1:num.class.cen){
    Cj = C.list[[j]]
    ind = row.posi[[j]]
    no.ind = length(ind)
    pjc = nrow(Cj)
    if(pjc == p){
      det.o[ind, i] = 1
      delta.o[ind, i] = 0
      mu.co = matrix(rep(mu[, i], no.ind), ncol=no.ind)
      Scc.o = Sigma[,,i]
    } else{
      O = O.list[[j]]
      OSO = O %*% Sigma[,,i] %*% t(O)
      cent.ind = cent[, ind, i]
      det.o[ind, i] = det(OSO)
      Soo = t(O) %*% solve(OSO) %*% O
      delta.o[ind, i] = colSums(cent.ind * (Soo %*% cent.ind))
      mu.co = Cj %*% matrix(rep(mu[, i], no.ind), ncol=no.ind) + Cj %*% Sigma[,,i] %*% Soo %*% cent.ind
      Scc.o = Cj %*% (Ip - Sigma[,,i] %*% Soo) %*% Sigma[,,i] %*% t(Cj)
    }
    if(distr=='MVN'){
      if(pjc == 1){
        for(s in 1: no.ind) cdf[ind[s], i] = pnorm(Cj%*%Yc[ind[s], ], mean=mu.co[,s], sd=sqrt(Scc.o))
      }
      if(pjc >= 2){
        for(s in 1: no.ind) cdf[ind[s], i] = ptmvnorm(lowerx=rep(-Inf,pjc), upperx=c(Cj%*%Yc[ind[s], ]), mean=c(mu.co[,s]), sigma=round(Scc.o, 6))
    }}
    if(distr=='MCN'){
      if(pjc == 0){
        for(s in 1: no.ind) w.cdf.pdf[ind[s], i] = nu[i] * dmvnorm(t(Yc[ind[s], ]), mean=mu[, i], sigma=round(Sigma[,,i]/rho[i],6)) + (1-nu[i]) * dmvnorm(t(Yc[ind[s], ]), mean=mu[, i], sigma=round(Sigma[,,i], 6)) 
      }
      if(pjc == 1){
        for(s in 1: no.ind){
        fyo = dmni(c(O%*%Yc[ind[s], ]), mu=c(O%*%mu[, i]), Sigma=round(OSO, 6), nu=nu[i], rho=rho[i], distr='MCN') 
        nu2.1 = nu[i] * dmvnorm(t(c(O%*%Yc[ind[s], ])), mean=c(O%*%mu[, i]), sigma=round(OSO/rho[i], 6)) / fyo
        w.cdf.pdf[ind[s], i] = (nu2.1 * pnorm(Cj%*%Yc[ind[s], ], mean=mu.co[,s], sd=sqrt(Scc.o/rho[i])) + (1-nu2.1) * pnorm(Cj%*%Yc[ind[s], ], mean=mu.co[,s], sd=sqrt(Scc.o))) * fyo
      }}                                   
      if(pjc >= 2 & pjc < (p-1)){
        for(s in 1: no.ind){ 
        fyo = dmni(c(O%*%Yc[ind[s], ]), mu=c(O%*%mu[, i]), Sigma=round(OSO,6), nu=nu[i], rho=rho[i], distr='MCN')
        nu2.1 = nu[i] * dmvnorm(t(c(O%*%Yc[ind[s], ])), mean=c(O%*%mu[, i]), sigma=round(OSO/rho[i], 6)) / fyo
        w.cdf.pdf[ind[s], i] = (nu2.1 * ptmvnorm(lowerx=rep(-Inf,pjc), upperx=c(Cj%*%Yc[ind[s], ]), mean=c(mu.co[,s]), sigma=round(Scc.o/rho[i],6)) + (1-nu2.1) * ptmvnorm(lowerx=rep(-Inf,pjc), upperx=c(Cj%*%Yc[ind[s], ]), mean=c(mu.co[,s]), sigma=round(Scc.o, 6))) * fyo  
      }} 
      if(pjc == (p-1)){
        for(s in 1: no.ind){ 
        fyo = dmni(c(O%*%Yc[ind[s], ]), mu=c(O%*%mu[, i]), Sigma=round(OSO, 6), nu=nu[i], rho=rho[i], distr='MCN')
        nu2.1 = nu[i] * dnorm(c(O%*%Yc[ind[s], ]), mean=c(O%*%mu[, i]), sd=sqrt(OSO/rho[i])) / fyo
        w.cdf.pdf[ind[s], i] = (nu2.1 * ptmvnorm(lowerx=rep(-Inf,pjc), upperx=c(Cj%*%Yc[ind[s], ]), mean=c(mu.co[,s]), sigma=round(Scc.o/rho[i],6)) + (1-nu2.1) * ptmvnorm(lowerx=rep(-Inf,pjc), upperx=c(Cj%*%Yc[ind[s], ]), mean=c(mu.co[,s]), sigma=round(Scc.o, 6))) * fyo  
      }} 
      if(pjc == p){
        for(s in 1: no.ind) w.cdf.pdf[ind[s], i] = nu[i] * ptmvnorm(lowerx=rep(-Inf,pjc), upperx=c(Cj%*%Yc[ind[s], ]), mean=c(mu.co[,s]), sigma=round(Scc.o/rho[i],6)) + (1-nu[i]) * ptmvnorm(lowerx=rep(-Inf,pjc), upperx=c(Cj%*%Yc[ind[s], ]), mean=c(mu.co[,s]), sigma=round(Scc.o,6))
    }}
    }
    if(distr=='MVN') wden[, i] = log(w[i]) -.5*po*log(2*pi) -.5*log(det.o[,i]) -.5*delta.o[,i] + log(cdf[,i])
    if(distr=='MCN') wden[, i] = log(w[i]) + log(w.cdf.pdf[, i])
    }
    max.wden = apply(wden, 1, max)
    wden = exp(wden - max.wden)
    indv.den = rowSums(wden)
    log.indv.den = log(indv.den) + max.wden
    loglik.new = sum(log.indv.den)    
    new.par = NULL
    for(i in 1: g){
     if(distr == 'MVN') new.par = c(new.par, c(w[i], mu[,i], Sigma[,,i][vechS]))
     if(distr == 'MCN') new.par = c(new.par, c(w[i], mu[,i], Sigma[,,i][vechS], nu[i], rho[i]))
    }
    iter.EST = rbind(iter.EST, c(iter, new.par))
    diff = loglik.new - loglik.old
    diff.par = mean(((old.par-new.par)/old.par)^2)
    iter.lnL = c(iter.lnL, loglik.new)

    if (iter%%per == 0){
    if(distr=='MVN') cat("iter", iter, ":", "log-likelihood =", loglik.new, "\t diff =", diff, "\t diff.para =", diff.par, '\t w=', w, '\t mu=', as.vector(mu), "\n")
    if(distr=='MCN') cat("iter", iter, ":", "log-likelihood =", loglik.new, "\t diff =", diff, "\t diff.para =", diff.par, '\t w=', w, '\t nu=', nu, '\t rho=', rho, '\t mu=', as.vector(mu), "\n")
    }
    if(diff < tol | diff.par < tol | iter > max.iter) break
    loglik.old = loglik.new
    old.par = new.par
    }
### Summarize results ###
    end = proc.time()[1]
    cat(paste(rep("-", 50), collapse = ""), "\n")
    run.sec = end - begin
# estimation
    EST = NULL
    if(distr=='MVN'){
     para = list(w=w, mu=mu, Sigma=Sigma)
     for(i in 1: g) EST = c(EST, c(w[i], mu[,i], Sigma[,,i][vechS]))
    }
    if(distr=='MCN'){
     para = list(w=w, mu=mu, Sigma=Sigma, nu=nu, rho=rho)
     for(i in 1: g) EST = c(EST, c(w[i], mu[,i], Sigma[,,i][vechS], nu[i], rho[i]))
    }
# model selection
    m = length(EST) - 1
    aic = 2 * m - 2 * loglik.new
    bic = m * log(n) - 2 * loglik.new
    model.inf = list(m = m, iter = iter, loglik = loglik.new, aic = aic, bic = bic, iter.lnL = iter.lnL, iter.EST = iter.EST, run=run.sec)
    cat("The CPU time takes", run.sec, "seconds.\n")
    cat('iter = ', iter, ',\t obs.loglik = ', loglik.new, '\t aic = ', aic, '\t bic = ', bic, '\t diff = ', diff, '\t diff.para =', diff.par, sep = '', '\n')
#    cat('w =', w, '\n')
#    cat('mu =', '\n')
#    print(mu)
#    cat('Sigma =', '\n')
#    print(Sigma)
#    if(distr == 'MCN') cat('nu =', nu, ',\t rho =', rho, '\n')
# clustering
    post.cls = matrix(apply(a.hat, 1, order), nrow = g)[g,]
    if(length(true.clus) != 0){
       CCR = 1 - classError(true.clus,post.cls)$errorRate
       ARI = adjustedRandIndex(true.clus, post.cls)
    } else{
      CCR = ARI = NULL
    } 
    pre.cls = list(wden = wden, a.hat = a.hat, post.cls = post.cls, CCR = CCR, ARI = ARI)
    cat('CCR =', CCR, ',\t ARI =', ARI, '\n')
    cat(paste(rep('=', 50), sep = '', collapse = ''), '\n')   
 # imputation
    z = matrix(0, n, g)
    for(i in 1: g) z[which(post.cls==i), i] = 1
    xi.hat1 = rowSums(a.hat * xi)
    xi.hat2 = rowSums(z * xi)
    xi = list(xi1 = xi.hat1, xi2 = xi.hat2)
    Ahat = Zhat = array(0, dim=c(n, p, g))
    for(i in 1: g){
     Ahat[,,i] = rep(a.hat[,i], p)
     Zhat[,,i] = rep(z[,i], p)
    }
    Yhat1 = apply(Ahat * yhat, 1:2, sum) 
    Yhat2 = apply(Zhat * yhat, 1:2, sum) 
    yhat = list(yhat1 = Yhat1, yhat2 = Yhat2)
    ychat1 = as.vector(t(Yhat1))[which(as.vector(t(cen)) == 1)]
    ychat2 = as.vector(t(Yhat2))[which(as.vector(t(cen)) == 1)]
    ychat = list(ychat1 = ychat1, ychat2 = ychat2)
    cat(paste(rep("-", 50), collapse = ""), "\n")
    if(SE == TRUE){
      IM = I.FMMCNC(para, Yc=Yc, cen=cen, distr=distr, EST=EST, wden=wden)
      return(list(model.inf = model.inf, para = para, EST = EST, IM = IM, pre.cls = pre.cls, xi = xi, yhat = yhat, ychat = ychat))
    } else{
      return(list(model.inf = model.inf, para = para, EST = EST, pre.cls = pre.cls, xi = xi, yhat = yhat, ychat = ychat))  
    }
}

# Louis' information matrix
I.FMMCNC = function(para.est, Yc, cen, distr=c('MVN','MCN'), EST, wden)
{
  GB = GenzBretz(maxpts = 5e4, abseps = 1e-9, releps = 0)
  n = nrow(Yc)
  p = ncol(Yc)
  vechS = vech.posi(p)
# parameter estimates
  w = para.est$w
  g = length(w)
  mu = para.est$mu
  Sigma = para.est$Sigma
  Sig.inv = array(NA, dim=c(p,p,g))
  for(i in 1: g) Sig.inv[,,i] = solve(Sigma[,,i])
  if(distr=='MVN'){ nu = rep(0, g); rho = rep(1, g)} 
  if(distr=='MCN'){ nu = para.est$nu; rho = para.est$rho} 

  Ip = diag(p)
  po = p - rowSums(cen)
  cen.subj = which(rowSums(cen) != 0)
  Nc = length(cen.subj)
  obs.subj = which(rowSums(cen) == 0)
  No = n - Nc
  ind.cen = colSums(t(cen) * 2 ^ (1:p - 1))
  num.class.cen = length(unique(ind.cen))
  row.posi = O.list = C.list = as.list(numeric(num.class.cen))
  uni.ind = unique(ind.cen)

  for(j in 1: num.class.cen){
    row.posi[[j]] = which(ind.cen == uni.ind[j])
    O.list[[j]] = matrix(Ip[!cen[row.posi[[j]][1],],], ncol = p)
    C.list[[j]] = matrix(Ip[which(cen[row.posi[[j]][1], ]==1), ], ncol=p)
  }

  indv.den = rowSums(wden)
  a.hat = wden / indv.den
  Ni = colSums(a.hat)

  cent = array(NA, dim = c(p, n, g))
  xi = matrix(1, n, g)
  yhat = xiyhat = array(NA, dim=c(n, p, g))
  y2hat = xiy2hat = array(NA, dim=c(p, p, n, g))
  for(i in 1: g){
   cent[,,i] = t(Yc) - mu[,i]

   for(j in 1:num.class.cen){
    Cj = C.list[[j]]
    ind = row.posi[[j]]               
    no.ind = length(ind)
    pjc = nrow(Cj)     #p-pjo
    
    if(pjc == 0){ 
      cent.ind = cent[, ind, i]
      d = diag(t(cent.ind) %*% Sig.inv[,,i] %*% cent.ind)
      if(distr=='MCN'){
        w.rho = w1 = numeric(no.ind) 
        for(s in 1: no.ind){
          dcn = dmni(c(Yc[ind[s], ]), mu=mu[, i], Sigma=round(Sigma[,,i], 6), nu=nu[i], rho=rho[i], distr='MCN')
          w.rho[s] = nu[i] * dmvnorm(t(c(Yc[ind[s], ])), mean=mu[, i], sigma=round(Sigma[,,i]/rho[i], 6)) / dcn 
          w1[s] = (1-nu[i]) * dmvnorm(t(c(Yc[ind[s], ])), mean=mu[, i], sigma=round(Sigma[,,i],6)) / dcn
        }
        xi[ind, i] = w.rho 
      }
      yhat[ind, , i] = Yc[ind, ]
      xiyhat[ind, ,i] = xi[ind, i] * Yc[ind, ] 
      for(s in 1: no.ind){
       y2hat[,,ind[s], i] = Yc[ind[s], ] %*% t(Yc[ind[s], ])
       xiy2hat[,,ind[s], i] = xi[ind[s], i]*Yc[ind[s], ] %*% t(Yc[ind[s], ])
      }} else{

      yc.hat = xi.yc.hat = array(NA, dim = c(pjc, no.ind, g))
      yc2.hat = xi.yc2.hat = array(NA, dim=c(pjc, pjc, no.ind, g))          
      if(pjc == p){
        mu.co = matrix(rep(mu[, i], no.ind), ncol=no.ind)
        Scc.o = Sigma[,,i]
        if(distr=='MVN'){
          for(s in 1: no.ind){
           EX = TMNI.moment(mu=mu.co[,s], Sigma=round(Scc.o, 6), distr='MVN', a.low=rep(-Inf, pjc), a.upp=c(Cj%*%Yc[ind[s], ])) 
           yc.hat[, s, i] = EX$EY
           xiy2hat[,, ind[s], i] = y2hat[,, ind[s], i] = EX$EYY  
          }
          xiyhat[ind, , i] = yhat[ind, , i] = t(t(Cj) %*% yc.hat[,,i])
        }
        if(distr=='MCN'){
          for(s in 1: no.ind){
           yjc = c(Cj%*%Yc[ind[s], ])
           pcn = 1 / pmni(lower=rep(-Inf, pjc), upper=yjc, mu=mu.co[,s], Sigma=round(Sigma[,,i], 6), nu=nu[i], rho=rho[i], distr='MCN')
           pn.rho = pmvnorm(lower = rep(-Inf, pjc), upper = yjc, mean=mu.co[,s], sigma = round(Sigma[,,i]/rho[i], 6), algorithm = GB)[1]
           xi[ind[s], i] = nu[i] * pn.rho * pcn 
           EX.rho = TMNI.moment(mu=mu.co[,s], Sigma=round(Scc.o, 6), nu=nu[i], rho=rho[i], distr='MCN', a.low=rep(-Inf, pjc), a.upp=yjc)
           EX1 = TMNI.moment(mu=mu.co[,s], Sigma=round(Scc.o/rho[i], 6), distr='MVN', a.low=rep(-Inf, pjc), a.upp=yjc) 
           xi.yc.hat[, s, i] = xi[ind[s], i] * EX1$EY
           xiy2hat[,, ind[s], i] = xi[ind[s], i] * EX1$EYY 
           yc.hat[, s, i] = EX.rho$EY
           y2hat[,, ind[s], i] = EX.rho$EYY
          }
          xiyhat[ind, , i] = t(t(Cj) %*% xi.yc.hat[,, i]) 
          yhat[ind, , i] = t(t(Cj) %*% yc.hat[,, i])
        }} else{
          O = O.list[[j]]
          OSO = O %*% Sigma[,,i] %*% t(O)
          Soo = t(O) %*% solve(OSO) %*% O
          cent.ind = cent[, ind, i]
          mu.co = Cj %*% matrix(rep(mu[, i], no.ind), ncol=no.ind) + Cj %*% Sigma[,,i] %*% Soo %*% cent.ind
          Scc.o = Cj %*% (Ip - Sigma[,,i] %*% Soo) %*% Sigma[,,i] %*% t(Cj)
          d = diag(t(cent.ind) %*% Soo %*% cent.ind)
          if(distr=='MVN'){
            for(s in 1: no.ind){
             EX = TMNI.moment(mu=mu.co[,s], Sigma=round(Scc.o, 6), distr='MVN', a.low=rep(-Inf, pjc), a.upp=c(Cj%*%Yc[ind[s], ])) 
             yc.hat[, s, i] = EX$EY
             yc2.hat[,, s, i] = EX$EYY
             yjo = O %*% Yc[ind[s],]
             xiy2hat[,, ind[s], i] = y2hat[,, ind[s], i] = xi[ind[s], i]*(t(O)%*% yjo %*% t(yjo) %*% O + t(O) %*% yjo %*% t(yc.hat[,s,i]) %*% Cj + t(Cj) %*% yc.hat[,s,i] %*% t(yjo) %*% O + t(Cj) %*% yc2.hat[,,s,i] %*% Cj)
            }
          xiyhat[ind, , i] = yhat[ind, , i] = xi[ind, i]*(t(t(O) %*% O %*% t(matrix(Yc[ind, ], ncol=p)) + t(Cj) %*% yc.hat[,,i])) 
          } 
          if(distr=='MCN'){
            w.rho = w1 = pn.rho = pn1 = w.cdf = numeric(no.ind)
            mu.o = O %*% mu[, i] 
            pjo = nrow(O)
            for(s in 1: no.ind){
              yjo = c(O %*% Yc[ind[s],])
              yjc = c(Cj %*% Yc[ind[s], ])
              dcn = dmni(yjo, mu=mu.o, Sigma=round(OSO, 6), nu=nu[i], rho=rho[i], distr='MCN')
              if(pjo == 1){
                w.rho[s] = nu[i] * dnorm(yjo, mean=c(mu.o), sd=sqrt(OSO/rho[i])) / dcn
                w1[s] = (1-nu[i]) * dnorm(yjo, mean=c(mu.o), sd=sqrt(OSO)) / dcn
              } else{
                w.rho[s] = nu[i] * dmvnorm(t(yjo), mean=mu.o, sigma=round(OSO/rho[i],6)) / dcn
                w1[s] = (1-nu[i]) * dmvnorm(t(yjo), mean=mu.o, sigma=OSO) / dcn
              }
              if(pjc == 1){
                pn.rho[s] = pnorm(yjc, mean=c(mu.co[,s]), sd=sqrt(Scc.o/rho[i]))
                pn1[s] = pnorm(yjc, mean=c(mu.co[,s]), sd=sqrt(Scc.o))
              } else{
                pn.rho[s] = pmvnorm(lower = rep(-Inf, pjc), upper = yjc, mean=mu.co[,s], sigma = round(Scc.o/rho[i],6), algorithm = GB)[1]
                pn1[s] = pmvnorm(lower = rep(-Inf, pjc), upper = yjc, mean=mu.co[,s], sigma = round(Scc.o, 6), algorithm = GB)[1]
              } 
              w.cdf[s] = w.rho[s]*pn.rho[s] + w1[s]*pn1[s] 
              EX.rho = TMNI.moment(mu=mu.co[,s], Sigma=round(Scc.o/rho[i], 6), distr='MVN', a.low=rep(-Inf, pjc), a.upp=yjc)
              EX1 = TMNI.moment(mu=mu.co[,s], Sigma=round(Scc.o, 6), distr='MVN', a.low=rep(-Inf, pjc), a.upp=yjc) 
              yc.hat[,s,i] = (w.rho[s]*pn.rho[s]*EX.rho$EY + w1[s]*pn1[s]*EX1$EY) / w.cdf[s]
              yc2.hat[,,s,i] = (w.rho[s]*pn.rho[s]*EX.rho$EYY + w1[s]*pn1[s]*EX1$EYY) / w.cdf[s]   
              xi[ind[s], i] = (w.rho[s]*pn.rho[s]) / w.cdf[s] 
              xi.yc.hat[, s, i] = xi[ind[s], i] * EX.rho$EY  
              xi.yc2.hat[,, s, i] = xi[ind[s], i] * EX.rho$EYY
              y2hat[,, ind[s], i] = (t(O)%*% yjo %*% t(yjo) %*% O) + (t(O) %*% yjo %*% t(yc.hat[,s, i]) %*% Cj) + (t(Cj) %*% yc.hat[, s, i] %*% t(yjo) %*% O) + (t(Cj) %*% yc2.hat[,, s, i] %*% Cj) 
              xiy2hat[,, ind[s], i] = xi[ind[s], i]*(t(O)%*% yjo %*% t(yjo) %*% O) + (t(O) %*% yjo %*% t(xi.yc.hat[,s,i]) %*% Cj) + (t(Cj) %*% xi.yc.hat[,s,i] %*% t(yjo) %*% O) + (t(Cj) %*% xi.yc2.hat[,,s,i] %*% Cj) 
            }
            yhat[ind, , i] = t(t(O) %*% O %*% t(matrix(Yc[ind, ], ncol=p))) + t(t(Cj) %*% yc.hat[,,i])
            xiyhat[ind, , i] = xi[ind, i]*t(t(O) %*% O %*% t(matrix(Yc[ind, ], ncol=p))) + t(t(Cj) %*% xi.yc.hat[,,i]) 
          }}
   }}}

# individual expected score vector
   q1 = p*(p+1)/2 
   m = g * (1 + p + q1)
   if(distr == 'MCN') m = m + 2 * g 

   dot.Sig = array(0, dim=c(p, p, q1))
   for(l in 1: q1){
    dotS = matrix(0, p, p)
    dotS[matrix(vechS[l, ], 1)]=dotS[matrix(rev(vechS[l, ]), 1)] = 1
    dot.Sig[,,l] = dotS
   }

   m1 = m / g
   aa = seq(0, m, m1)
   EIs = matrix(NA, nrow=m, ncol=n)
   for(i in 1: g){
      EIs[(aa[i]+1): aa[(i+1)], ][1, ] = a.hat[, i]/w[i]
      for(j in 1: n){
        EIs[(aa[i]+1): aa[(i+1)], ][2:(1+p), j] = (a.hat[j, i] * Sig.inv[,,i]) %*% (yhat[j, ,i] - mu[,i] + (rho[i]-1)*(xiyhat[j,,i] - xi[j,i]*mu[,i])) 
        A1 = y2hat[,,j,i] - yhat[j,,i] %*% t(mu[,i]) - mu[,i] %*% t(yhat[j,,i]) + mu[,i] %*% t(mu[,i])
        A2 = xiy2hat[,,j,i] - xiyhat[j,,i] %*% t(mu[,i]) - mu[,i] %*% t(xiyhat[j,,i]) + (xi[j,i] * mu[,i]) %*% t(mu[,i])
       for(k in 1: q1){
        SSdot = Sig.inv[,,i] %*% dot.Sig[,,k]
        EIs[(aa[i]+1): aa[(i+1)], ][(1+p+k), j] = 0.5 * a.hat[j,i] * (tr((A1 + (rho[i]-1) * A2) %*% SSdot %*% Sig.inv[,,i]) - tr(SSdot))
       }
       if(distr == 'MCN') EIs[(aa[i]+1): aa[(i+1)], ][(3+p+q1), j] = 0.5 * a.hat[j,i] * (p*xi[j,i]/rho[i] - tr(A2 %*% Sig.inv[,,i]))
      }
      if(distr == 'MCN') EIs[(aa[i]+1): aa[(i+1)], ][(2+p+q1), ] = a.hat[,i] * (xi[,i] - nu[i]) / (nu[i] * (1-nu[i]))
   }
   
# Meilijson (1989) formula
   I.theta = EIs[,1] %*% t(EIs[,1])
   for(j in 2: n) I.theta = I.theta + EIs[,j] %*% t(EIs[,j])
   sd.theta = rep(NA, m)
   V.theta = try(solve(I.theta), silent=F)
   if(class(V.theta)[1] != "try-error") sd.theta = sqrt(diag(V.theta))
   if(distr == 'MVN'){
     EST = sd.mvn = NULL
     for(i in 1: g){
      EST = c(EST, c(w[i], mu[,i], Sigma[,,i][vechS], nu[i], rho[i]))
      sd.mvn = c(sd.mvn, c(sd.theta[(aa[i]+1): aa[(i+1)]], 0, 0))
     }
     sd.theta = sd.mvn
   }
   theta.hat = rbind(c(EST), sd.theta)
   if(g==1) theta.hat = theta.hat[,-1]
   return(list(theta.hat=theta.hat, I.theta = I.theta, V.theta = V.theta, SD=sd.theta))
}
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    