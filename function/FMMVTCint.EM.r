# ECM Algorithm for FM-MVT and FM-MVN Models
FMMVTCint.EM = function(Yc1, Yc2, cen, g, distr=c('MVN','MVT'), true.clus=NULL, tol=1e-6, max.iter = 1000, per=10, SE=FALSE)
{
  require(mvtnorm)
  GB = GenzBretz(maxpts = 5e4, abseps = 1e-9, releps = 0)
  begin = proc.time()[1]
  distr = distr[1]
 # initial values:
  n = nrow(Yc2)
  p = ncol(Yc2)
  vechS = vech.posi(p)

  if(g == 1){ 
    clus = rep(1, n)
  } else{
    if(length(true.clus) == 0 | length(unique(true.clus)) != g){
     kmY = kmeans(Yc2, g, nstart = 25)
     clus = kmY$cluster
    } else{ clus = true.clus}
  }
  rk1 = rev(rank(table(clus)/n))
  w = rep(1/g, g)
  mu = matrix(NA, p, g)
  Sigma = Sig.inv = array(NA, dim = c(p,p,g))
  for(i in 1:g)
  {
    w[i] = nrow(Yc2[clus == rk1[i],])/n
    mu[,i] = colMeans(Yc2[clus == rk1[i],])
    Sigma[,,i] = cov(Yc2[clus == rk1[i],])
  }
  if(distr == 'MVT') nu = rep(5, g)
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
  
# old observed-data log-likelihood
  wden = matrix(NA, n, g)
  det.o = matrix(NA, n, g)
  delta.o = matrix(NA, n, g)
  pdf = cdf = w.cdf.pdf = matrix(1, n, g)
  cent = array(NA, dim = c(p, n, g))

  for(i in 1: g){
   cent[,,i] = t(Yc2) - mu[,i]

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
        for(s in 1: no.ind) cdf[ind[s], i] = pnorm(Cj%*%Yc2[ind[s], ], mean=mu.co[,s], sd=sqrt(Scc.o)) - pnorm(Cj%*%Yc1[ind[s], ], mean=mu.co[,s], sd=sqrt(Scc.o))
      }
      if(pjc >= 2){
        for(s in 1: no.ind) cdf[ind[s], i] = ptmvnorm(lowerx=c(Cj%*%Yc1[ind[s], ]), upperx=c(Cj%*%Yc2[ind[s], ]), mean=c(mu.co[,s]), sigma=round(Scc.o, 6))
    }}
    if(distr=='MVT'){
      DF = round(nu[i]+po[ind][1])
      if(pjc == 1){
        for(s in 1: no.ind) cdf[ind[s], i] = pt(((Cj%*%Yc2[ind[s], ])-mu.co[,s])/sqrt(Scc.o),df=DF) - pt(((Cj%*%Yc1[ind[s], ])-mu.co[,s])/sqrt(Scc.o),df=DF)
      }
      if(pjc >= 2){
        for(s in 1: no.ind) cdf[ind[s], i] = ptmvt(lowerx=c(Cj%*%Yc1[ind[s], ]), upperx=c(Cj%*%Yc2[ind[s], ]), mean=c(mu.co[,s]), sigma=Scc.o,df=DF)
    }}
    }                        
    if(distr=='MVN') wden[, i] = log(w[i]) -.5*po*log(2*pi) -.5*log(det.o[,i]) -.5*delta.o[,i] + log(cdf[,i])
    if(distr=='MVT') wden[, i] = log(w[i]) + lgamma((nu[i]+po)/2) - lgamma(nu[i]/2) - .5*po*log(pi*nu[i])- .5*log(det.o[,i])  - .5*(nu[i]+po)*log(1+delta.o[,i]/nu[i]) +log(cdf[,i])
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
    if(distr == 'MVT'){
      for(i in 1: g) old.par = c(old.par, c(w[i], mu[,i], Sigma[,,i][vechS], nu[i]))
    }
    iter.EST = c(iter, old.par)
    cat(paste(rep("=", 50), collapse = ""), "\n")
    cat("EM is running for the", g, "-component FM censored model with", distr, "distribution: censoring proportion", mean(cen)*100, '%.', "\n")
    cat("Initial log-likelihood = ", loglik.old, "\n")

    tau = matrix(1, n, g)
    tauyhat = array(NA, dim=c(n, p, g))
    tauy2hat = array(NA, dim=c(p, p, n, g))
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
          d = diag(t(cent.ind)%*%solve(Sigma[,,i])%*%cent.ind)
          if(distr=='MVT') tau[ind, i] = (nu[i] + p)/(nu[i] + d)
          tauyhat[ind, ,i] = tau[ind, i] * Yc2[ind, ]
          for(s in 1: no.ind) tauy2hat[,,ind[s],i] = tau[ind[s],i]*Yc2[ind[s], ] %*% t(Yc2[ind[s], ])
        } else{

          yc.hat = tau.yc.hat =  array(NA, dim = c(pjc, no.ind, g))
          yc2.hat = tau.yc2.hat = array(NA, dim=c(pjc, pjc, no.ind, g))          
          if(pjc == p){
            mu.co = matrix(rep(mu[, i], no.ind), ncol=no.ind)
            Scc.o = Sigma[,,i] 

            if(distr=='MVN'){
              for(s in 1: no.ind){
               EX = TMNI.moment(mu=mu.co[,s], Sigma=Scc.o, distr='MVN', a.low=c(Cj%*%Yc1[ind[s], ]), a.upp=c(Cj%*%Yc2[ind[s], ])) 
               yc.hat[,s,i] = EX$EY
               tauy2hat[,,ind[s],i] = EX$EYY  
             }
             tauyhat[ind, ,i] = t(t(Cj) %*% yc.hat[,,i])
            }
            if(distr=='MVT'){
              for(s in 1: no.ind){
               yjc1 = c(Cj%*%Yc1[ind[s], ])
               yjc2 = c(Cj%*%Yc2[ind[s], ])
               tau[ind[s], i] = ptmvt(lowerx=yjc1, upperx=yjc2, mean=mu[,i], sigma=nu[i]/(nu[i]+2)*Sigma[,,i], df=round(nu[i]+2))/ptmvt(lowerx=yjc1, upperx=yjc2, mean=mu[,i], sigma=Sigma[,,i], df=round(nu[i]))
               EX = TMNI.moment(mu=mu.co[,s], Sigma=nu[i]/(nu[i]+2)*Sigma[,,i], nu=nu[i]+2, distr='MVT', a.low=yjc1, a.upp=yjc2)
               yc.hat[,s,i] = EX$EY
               tauy2hat[,,ind[s],i] = tau[ind[s], i] * EX$EYY  
              }
              tauyhat[ind,,i] = tau[ind,i]*(t(t(Cj) %*% yc.hat[,,i])) 
            }}else{
            O = O.list[[j]]
            OSO = O %*% Sigma[,,i] %*% t(O)
            Soo = t(O) %*% solve(OSO) %*% O
            cent.ind = cent[, ind, i]
            mu.co = Cj %*% matrix(rep(mu[, i], no.ind), ncol=no.ind) + Cj %*% Sigma[,,i] %*% Soo %*% cent.ind
            Scc.o = Cj %*% (Ip - Sigma[,,i] %*% Soo) %*% Sigma[,,i] %*% t(Cj)
            d=diag(t(cent.ind)%*%Soo%*%cent.ind)
            if(distr=='MVN'){
              for(s in 1: no.ind){
               EX = TMNI.moment(mu=mu.co[,s], Sigma=Scc.o, distr='MVN', a.low=c(Cj%*%Yc1[ind[s], ]), a.upp=c(Cj%*%Yc2[ind[s], ])) 
               yc.hat[,s,i] = EX$EY
               yc2.hat[,,s,i] = EX$EYY
               yjo = O %*% Yc2[ind[s],]
               tauy2hat[,,ind[s],i] = tau[ind[s],i]*(t(O)%*% yjo %*% t(yjo) %*% O + t(O) %*% yjo %*% t(yc.hat[,s,i]) %*% Cj + t(Cj) %*% yc.hat[,s,i] %*% t(yjo) %*% O + t(Cj) %*% yc2.hat[,,s,i] %*% Cj)
              }
            tauyhat[ind,,i] = tau[ind,i]*(t(t(O) %*% O %*% t(matrix(Yc2[ind, ], ncol=p)) + t(Cj) %*% yc.hat[,,i])) 
            } 
            if(distr=='MVT'){
              tmp = (nu[i] + po[ind]) / (nu[i] + d)
              DF = round(nu[i]+po[ind][1])
              for(s in 1: no.ind){
               yjc1 = c(Cj%*%Yc1[ind[s], ])
               yjc2 = c(Cj%*%Yc2[ind[s], ])
               newScc.o = rep(1/tmp[s], pjc) * Scc.o 
               if(pjc == 1){ cdf.rate = (pt((yjc2-c(mu.co[,s]))/sqrt((DF/(DF+2))*newScc.o), df=DF+2)-pt((yjc1-c(mu.co[,s]))/sqrt((DF/(DF+2))*newScc.o), df=DF+2)) / (pt((yjc2-c(mu.co[,s]))/sqrt(newScc.o), df=DF)-pt((yjc1-c(mu.co[,s]))/sqrt(newScc.o), df=DF))
               } else cdf.rate = ptmvt(lowerx=yjc1, upperx=yjc2, mean=c(mu.co[,s]), sigma=rep((DF/(DF+2)),pjc)*newScc.o, df=DF+2) / ptmvt(lowerx=yjc1, upperx=yjc2, mean=c(mu.co[,s]), sigma=newScc.o, df=DF)
               tau[ind[s],i] = tmp[s] * cdf.rate
               EX = TMNI.moment(mu=mu.co[,s], Sigma=rep(DF/(DF+2), pjc)*newScc.o, nu=DF+2, distr='MVT', a.low=yjc1, a.upp=yjc2)
               yc.hat[,s,i] = EX$EY
               yc2.hat[,,s,i] = EX$EYY
               yjo = O %*% Yc2[ind[s],]
               tauy2hat[,,ind[s],i] = tau[ind[s],i]*(t(O)%*% yjo %*% t(yjo) %*% O + t(O) %*% yjo %*% t(yc.hat[,s,i]) %*% Cj + t(Cj) %*% yc.hat[,s,i] %*% t(yjo) %*% O + t(Cj) %*% yc2.hat[,,s,i] %*% Cj)
              }
            tauyhat[ind,,i] = tau[ind,i]*(t(t(O) %*% O %*% t(matrix(Yc2[ind, ], ncol=p)) + t(Cj) %*% yc.hat[,,i])) 
            }
          }
      }}}    
# M-step
     for(i in 1: g){ 
#      tauyhat = matrix(as.numeric(tauyhat),n,p)
      mu[,i] = colSums(a.hat[,i] * tauyhat[,,i])/sum(a.hat[,i] * tau[,i])
      E1 = array(NA, dim=c(p, p, n))
      for(j in 1: n) E1[,,j] = a.hat[j, i] * tauy2hat[,,j,i]
      Ta = matrix(rep(a.hat[,i], p), ncol=p)
      Sigma[,,i] = (apply(E1, 1:2, sum) - colSums(Ta * tauyhat[,,i])%*%t(mu[,i]) - mu[,i]%*%t(colSums(Ta * tauyhat[,,i])) + sum(a.hat[,i]*tau[,i])*mu[,i]%*%t(mu[,i])) / Ni[i]
      Sigma[,,i] = (Sigma[,,i] + t(Sigma[,,i]))/2
      if(distr=='MVT'){
        nu[i] = nlminb(start = nu[i], objective = FMMVTCint.nu.fn, lower = 2+1e-6, upper =Inf, gup=i, w=w, mu=mu, Sigma=Sigma, nu=nu, Yc1=Yc1, Yc2=Yc2, cen=cen)$par     
      }}

# new observed-data log-likelihood
  for(i in 1: g){
   cent[,,i] = t(Yc2) - mu[,i]

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
        for(s in 1: no.ind) cdf[ind[s], i] = pnorm(Cj%*%Yc2[ind[s], ], mean=mu.co[,s], sd=sqrt(Scc.o)) - pnorm(Cj%*%Yc1[ind[s], ], mean=mu.co[,s], sd=sqrt(Scc.o))
      }
      if(pjc >= 2){
        for(s in 1: no.ind) cdf[ind[s], i] = ptmvnorm(lowerx=c(Cj%*%Yc1[ind[s], ]), upperx=c(Cj%*%Yc2[ind[s], ]), mean=c(mu.co[,s]), sigma=round(Scc.o, 6))
       }}
    if(distr=='MVT'){
      DF = round(nu[i]+po[ind][1])
      if(pjc == 1){
        for(s in 1: no.ind) cdf[ind[s], i] = pt(((Cj%*%Yc2[ind[s], ])-mu.co[,s])/sqrt(Scc.o),df=DF) - pt(((Cj%*%Yc1[ind[s], ])-mu.co[,s])/sqrt(Scc.o),df=DF)
      }
      if(pjc >= 2){
        for(s in 1: no.ind) cdf[ind[s], i] = ptmvt(lowerx=c(Cj%*%Yc1[ind[s], ]), upperx=c(Cj%*%Yc2[ind[s], ]), mean=c(mu.co[,s]), sigma=Scc.o,df=DF)
    }}
    }                        
    if(distr=='MVN') wden[, i] = log(w[i]) -.5*po*log(2*pi) -.5*log(det.o[,i]) -.5*delta.o[,i] + log(cdf[,i])
    if(distr=='MVT') wden[, i] = log(w[i]) + lgamma((nu[i]+po)/2) - lgamma(nu[i]/2) - .5*po*log(pi*nu[i])- .5*log(det.o[,i])  - .5*(nu[i]+po)*log(1+delta.o[,i]/nu[i]) +log(cdf[,i])
    }
    max.wden = apply(wden, 1, max)
    wden = exp(wden - max.wden)
    indv.den = rowSums(wden)
    log.indv.den = log(indv.den) + max.wden
    loglik.new = sum(log.indv.den)    
    new.par = NULL
    for(i in 1: g){
     if(distr == 'MVN') new.par = c(new.par, c(w[i], mu[,i], Sigma[,,i][vechS]))
     if(distr == 'MVT') new.par = c(new.par, c(w[i], mu[,i], Sigma[,,i][vechS], nu[i]))
    }
    iter.EST = rbind(iter.EST, c(iter, new.par))
    diff = loglik.new - loglik.old
    diff.par = mean(((old.par-new.par)/old.par)^2)
    iter.lnL = c(iter.lnL, loglik.new)

    if (iter%%per == 0){
    if(distr=='MVN') cat("iter", iter, ":", "log-likelihood =", loglik.new, "\t diff =", diff, "\t diff.para =", diff.par, '\t w=', w, '\t mu=', as.vector(mu), "\n")
    if(distr=='MVT') cat("iter", iter, ":", "log-likelihood =", loglik.new, "\t diff =", diff, "\t diff.para =", diff.par, '\t w=', w, '\t nu=', nu, '\t mu=', as.vector(mu), "\n")
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
    if(g == 1){
     rk = 1
    } else{
    rk = order(w, decreasing = TRUE)
    w = w[rk]
    mu = mu[, rk]
    Sigma = Sigma[,,rk]
    }
    EST = NULL
    if(distr=='MVN'){
     para.est = list(w=w, mu=mu, Sigma=Sigma)
     for(i in 1: g) EST = c(EST, c(w[i], mu[,i], Sigma[,,i][vechS]))
    }
    if(distr=='MVT'){
     nu = nu[rk]
     para.est = list(w=w, mu=mu, Sigma=Sigma, nu=nu)
     for(i in 1: g) EST = c(EST, c(w[i], mu[,i], Sigma[,,i][vechS], nu[i]))
    }
# model selection
    m = length(EST) - 1
    aic = 2 * m - 2 * loglik.new
    bic = m * log(n) - 2 * loglik.new
    model.inf = list(m = m, iter = iter, loglik = loglik.new, aic = aic, bic = bic, iter.lnL = iter.lnL, iter.EST = iter.EST, run=run.sec)
    cat("The CPU time takes", run.sec, "seconds.\n")
    cat('iter = ', iter, ',\t obs.loglik = ', loglik.new, '\t aic = ', aic, '\t bic = ', bic, sep = '', '\n')
#    cat('w =', w, '\n')
#    cat('mu =', '\n')
#    print(mu)
#    cat('Sigma =', '\n')
#    print(Sigma)
#    if(distr == 'MVT') cat('nu =', nu, '\n')
# clustering
    if(g == 1) post.cls = matrix(apply(a.hat, 1, order), nrow = g)[g,]
    if(g != 1) post.cls = matrix(apply(a.hat[,rk], 1, order), nrow = g)[g,]
    if(length(true.clus) != 0){
       CCR = 1 - classError(true.clus,post.cls)$errorRate
       ARI = adjustedRandIndex(true.clus, post.cls)
    } else{
      CCR = ARI = NULL
    } 
    if(g == 1) pre.cls = list(wden = wden, a.hat = a.hat, post.cls = post.cls, CCR = CCR, ARI = ARI)
    if(g != 1) pre.cls = list(wden = wden[,rk], a.hat = a.hat[,rk], post.cls = post.cls, CCR = CCR, ARI = ARI)
    cat(paste(rep('=', 50), sep = '', collapse = ''), '\n')   
 # imputation
    z = matrix(0, n, g)
    for(i in 1: g) z[which(post.cls==i), i] = 1
    Ahat = Zhat = array(0, dim=c(n, p, g))
    if(g == 1){
    tau.hat1 = rowSums(a.hat * tau)
    tau.hat2 = rowSums(z * tau)
    for(i in 1: g){
     Ahat[,,i] = rep(a.hat[,i]/tau[,i], p)
     Zhat[,,i] = rep(z[,i]/tau[,i], p)
    }
    Yhat1 = apply(Ahat * tauyhat, 1:2, sum) 
    Yhat2 = apply(Zhat * tauyhat, 1:2, sum) 
    } else{
    tau.hat1 = rowSums(a.hat[,rk] * tau[,rk])
    tau.hat2 = rowSums(z * tau[,rk])
    for(i in 1: g){
     Ahat[,,i] = rep(a.hat[,rk][,i]/tau[,rk][,i], p)
     Zhat[,,i] = rep(z[,i]/tau[,rk][,i], p)
    }
    Yhat1 = apply(Ahat * tauyhat[,,rk], 1:2, sum) 
    Yhat2 = apply(Zhat * tauyhat[,,rk], 1:2, sum) 
    }
    tau.hat = list(tau1 = tau.hat1, tau2 = tau.hat2)
    yhat = list(yhat1 = Yhat1, yhat2 = Yhat2)
    ychat1 = as.vector(t(Yhat1))[which(as.vector(t(cen)) == 1)]
    ychat2 = as.vector(t(Yhat2))[which(as.vector(t(cen)) == 1)]
    ychat = list(ychat1 = ychat1, ychat2 = ychat2)
    if(SE == TRUE){
     if(g == 1) IM = I.FMMVTCint(para.est, Yc1=Yc1, Yc2=Yc2, cen=cen, distr=distr, EST=EST, wden=wden)
     if(g != 1) IM = I.FMMVTCint(para.est, Yc1=Yc1, Yc2=Yc2, cen=cen, distr=distr, EST=EST, wden=wden[, rk])
     return(list(iter=iter, model.inf = model.inf, para = para.est, EST = EST, pre.cls=pre.cls, IM=IM, yhat = yhat, ychat = ychat, tau=tau.hat))
    } else{
     return(list(iter=iter, model.inf = model.inf, para = para.est, EST = EST, pre.cls=pre.cls, yhat = yhat, ychat = ychat, tau=tau.hat))
    }
}

# MVT log-likelihood function for nu
FMMVTCint.nu.fn = function(par, gup, w, mu, Sigma, nu, Yc1, Yc2, cen)
{
  nu[gup] = par
  n = nrow(Yc2)
  p = ncol(Yc2)
  g = length(w)
  vechS = vech.posi(p)

  Ip = diag(p)
  po = p - rowSums(cen)
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

  for(i in 1: g){
   cent[,,i] = t(Yc2) - mu[,i]

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
    DF = round(nu[i]+po[ind][1])
    if(pjc == 1){
      for(s in 1: no.ind) cdf[ind[s], i] = pt(((Cj%*%Yc2[ind[s], ])-mu.co[,s])/sqrt(Scc.o),df=DF) - pt(((Cj%*%Yc1[ind[s], ])-mu.co[,s])/sqrt(Scc.o),df=DF)
    }
    if(pjc >= 2){
      for(s in 1: no.ind) cdf[ind[s], i] = ptmvt(lowerx=c(Cj%*%Yc1[ind[s], ]), upperx=c(Cj%*%Yc2[ind[s], ]), mean=c(mu.co[,s]), sigma=Scc.o,df=DF)
    }}                        
    wden[, i] = log(w[i]) + lgamma((nu[i]+po)/2) - lgamma(nu[i]/2) - .5*po*log(pi*nu[i])- .5*log(det.o[,i])  - .5*(nu[i]+po)*log(1+delta.o[,i]/nu[i]) +log(cdf[,i])
    }
    max.wden = apply(wden, 1, max)
    wden = exp(wden - max.wden)
    indv.den = rowSums(wden)
    log.indv.den = log(indv.den) + max.wden
    return(-sum(log.indv.den))
}

# Louis' information matrix
I.FMMVTCint = function(para.est, Yc1, Yc2, cen, distr=c('MVN','MVT'), EST, wden)
{
  GB = GenzBretz(maxpts = 5e4, abseps = 1e-9, releps = 0)
  n = nrow(Yc2)
  p = ncol(Yc2)
  vechS = vech.posi(p)
# parameter estimates
  w = para.est$w
  g = length(w)
  mu = para.est$mu
  Sigma = para.est$Sigma
  Sig.inv = array(NA, dim = c(p,p,g))
  for(i in 1:g) Sig.inv[,,i] = solve(Sigma[,,i])
  if(distr=='MVN') nu = rep(100, g)
  if(distr=='MVT') nu = para.est$nu 

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

  tau = ka = matrix(1, n, g)
  tauyhat = array(NA, dim=c(n, p, g))
  tauy2hat = array(NA, dim=c(p, p, n, g))
  cent = array(NA, dim = c(p, n, g))

  for(i in 1: g){
   cent[,,i] = t(Yc2) - mu[,i]

   for(j in 1:num.class.cen){
     Cj = C.list[[j]]
     ind = row.posi[[j]]               
     no.ind = length(ind)
     pjc = nrow(Cj)     #p-pjo
        
     if(pjc == 0){ 
        cent.ind = cent[, ind, i]
        d = diag(t(cent.ind)%*%solve(Sigma[,,i])%*%cent.ind)
        if(distr=='MVT'){
         tau[ind, i] = (nu[i] + p)/(nu[i] + d)
         ka[ind, i] = digamma((nu[i] + p)/2) - log((nu[i] + d)/2)
        }
        tauyhat[ind, ,i] = tau[ind, i] * Yc2[ind, ]
        for(s in 1: no.ind) tauy2hat[,,ind[s],i] = tau[ind[s],i]*Yc2[ind[s], ] %*% t(Yc2[ind[s], ])
     } else{

     yc.hat = tau.yc.hat =  array(NA, dim = c(pjc, no.ind, g))
     yc2.hat = tau.yc2.hat = array(NA, dim=c(pjc, pjc, no.ind, g))          
     if(pjc == p){
        mu.co = matrix(rep(mu[, i], no.ind), ncol=no.ind)
        Scc.o = Sigma[,,i] 

        if(distr=='MVN'){
          for(s in 1: no.ind){
           EX = TMNI.moment(mu=mu.co[,s], Sigma=Scc.o, distr='MVN', a.low=c(Cj%*%Yc1[ind[s], ]), a.upp=c(Cj%*%Yc2[ind[s], ])) 
           yc.hat[,s,i] = EX$EY
           tauy2hat[,,ind[s],i] = EX$EYY  
        }
        tauyhat[ind, ,i] = t(t(Cj) %*% yc.hat[,,i])
        }
        if(distr=='MVT'){
          for(s in 1: no.ind){
           yjc1 = c(Cj%*%Yc1[ind[s], ])
           yjc2 = c(Cj%*%Yc2[ind[s], ])
           tau[ind[s], i] = ptmvt(lowerx=yjc1, upperx=yjc2, mean=mu[,i], sigma=nu[i]/(nu[i]+2)*Sigma[,,i], df=round(nu[i]+2))/ptmvt(lowerx=yjc1, upperx=yjc2, mean=mu[,i], sigma=Sigma[,,i], df=round(nu[i]))
           EX = TMNI.moment(mu=mu.co[,s], Sigma=nu[i]/(nu[i]+2)*Sigma[,,i], nu=nu[i]+2, distr='MVT', a.low=yjc1, a.upp=yjc2)
           yc.hat[,s,i] = EX$EY
           tauy2hat[,,ind[s],i] = tau[ind[s], i] * EX$EYY  
           yjc.gen = t(rtmni(10, mu=mu.co[,s], Sigma=nu[i]/(nu[i]+2)*Sigma[,,i], nu=round(nu[i]+2), distr='MVT', lower=yjc1, upper=yjc2)$Y)
           cent.yc = matrix((yjc.gen - mu[,i]), nrow=pjc)
           dc = diag(t(cent.yc)%*%Sig.inv[,,i]%*%cent.yc)
           ka[ind[s],i] = digamma((nu[i]+p)/2) - mean(log((nu[i]+dc)/2))
          }
          tauyhat[ind,,i] = tau[ind,i]*(t(t(Cj) %*% yc.hat[,,i])) 
        }}else{
        O = O.list[[j]]
        OSO = O %*% Sigma[,,i] %*% t(O)
        Soo = t(O) %*% solve(OSO) %*% O
        cent.ind = cent[, ind, i]
        mu.co = Cj %*% matrix(rep(mu[, i], no.ind), ncol=no.ind) + Cj %*% Sigma[,,i] %*% Soo %*% cent.ind
        Scc.o = Cj %*% (Ip - Sigma[,,i] %*% Soo) %*% Sigma[,,i] %*% t(Cj)
        d=diag(t(cent.ind)%*%Soo%*%cent.ind)
        if(distr=='MVN'){
          for(s in 1: no.ind){
          EX = TMNI.moment(mu=mu.co[,s], Sigma=Scc.o, distr='MVN', a.low=c(Cj%*%Yc1[ind[s], ]), a.upp=c(Cj%*%Yc2[ind[s], ])) 
          yc.hat[,s,i] = EX$EY
          yc2.hat[,,s,i] = EX$EYY
          yjo = O %*% Yc2[ind[s],]
          tauy2hat[,,ind[s],i] = tau[ind[s],i]*(t(O)%*% yjo %*% t(yjo) %*% O + t(O) %*% yjo %*% t(yc.hat[,s,i]) %*% Cj + t(Cj) %*% yc.hat[,s,i] %*% t(yjo) %*% O + t(Cj) %*% yc2.hat[,,s,i] %*% Cj)
        }
        tauyhat[ind,,i] = tau[ind,i]*(t(t(O) %*% O %*% t(matrix(Yc2[ind, ], ncol=p)) + t(Cj) %*% yc.hat[,,i])) 
        } 
        if(distr=='MVT'){
          tmp = (nu[i] + po[ind]) / (nu[i] + d)
          DF = round(nu[i]+po[ind][1])
          for(s in 1: no.ind){
           yjc1 = c(Cj%*%Yc1[ind[s], ])
           yjc2 = c(Cj%*%Yc2[ind[s], ])
           newScc.o = rep(1/tmp[s], pjc) * Scc.o 
           if(pjc == 1){ cdf.rate = (pt((yjc2-c(mu.co[,s]))/sqrt((DF/(DF+2))*newScc.o), df=DF+2)-pt((yjc1-c(mu.co[,s]))/sqrt((DF/(DF+2))*newScc.o), df=DF+2)) / (pt((yjc2-c(mu.co[,s]))/sqrt(newScc.o), df=DF) - pt((yjc1-c(mu.co[,s]))/sqrt(newScc.o), df=DF))
           } else cdf.rate = ptmvt(lowerx=yjc1, upperx=yjc2, mean=c(mu.co[,s]), sigma=rep((DF/(DF+2)),pjc)*newScc.o, df=DF+2) / ptmvt(lowerx=yjc1, upperx=yjc2, mean=c(mu.co[,s]), sigma=newScc.o, df=DF)
           tau[ind[s],i] = tmp[s] * cdf.rate
           EX = TMNI.moment(mu=mu.co[,s], Sigma=rep(DF/(DF+2), pjc)*newScc.o, nu=DF+2, distr='MVT', a.low=yjc1, a.upp=yjc2)
           yc.hat[,s,i] = EX$EY
           yc2.hat[,,s,i] = EX$EYY
           yjo = O %*% Yc2[ind[s],]
           tauy2hat[,,ind[s],i] = tau[ind[s],i]*(t(O)%*% yjo %*% t(yjo) %*% O + t(O) %*% yjo %*% t(yc.hat[,s,i]) %*% Cj + t(Cj) %*% yc.hat[,s,i] %*% t(yjo) %*% O + t(Cj) %*% yc2.hat[,,s,i] %*% Cj)
           yjc.gen = t(rtmni(10, mu=mu.co[,s], Sigma=rep(DF/(DF+2), pjc)*newScc.o, nu=DF+2, distr='MVT', lower=yjc1, upper=yjc2)$Y)
           cent.yc = matrix((yjc.gen - mu.co[,s]), nrow=pjc)
           d.co = diag(t(cent.yc)%*%solve(Scc.o)%*%cent.yc)
           ka[ind[s],i] = digamma((nu[i]+p)/2) - mean(log((nu[i]+d[s]+d.co)/2))  
           }
         tauyhat[ind,,i] = tau[ind,i]*(t(t(O) %*% O %*% t(matrix(Yc2[ind, ], ncol=p)) + t(Cj) %*% yc.hat[,,i])) 
         }
      }}}}    

# individual expected score vector
   q1 = p*(p+1)/2 
   m = g * (1 + p + q1)
   if(distr == 'MVT') m = m + g 

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
      EIs[(aa[i]+1): aa[(i+1)], ][1, ] = a.hat[, i]/w[i] - a.hat[,g]/w[g]
      for(j in 1: n){
        EIs[(aa[i]+1): aa[(i+1)], ][2:(1+p), j] = (a.hat[j, i] * Sig.inv[,,i]) %*% (tauyhat[j, ,i] - tau[j,i]*mu[,i]) 
        A2 = tauy2hat[,,j,i] - tauyhat[j,,i] %*% t(mu[,i]) - mu[,i] %*% t(tauyhat[j,,i]) + (tau[j,i] * mu[,i]) %*% t(mu[,i])
       for(k in 1: q1){
        SSdot = Sig.inv[,,i] %*% dot.Sig[,,k]
        EIs[(aa[i]+1): aa[(i+1)], ][(1+p+k), j] = .5 * a.hat[j,i] * (tr(A2 %*% SSdot %*% Sig.inv[,,i]) - tr(SSdot))
       }}
      if(distr == 'MVT') EIs[(aa[i]+1): aa[(i+1)], ][(2+p+q1), ] = .5 * a.hat[,i] * (ka[,i] - tau[,i] - digamma(nu[i]/2)/gamma(nu[i]/2) - log(nu[i]/2) -1)
   }
   
# Meilijson (1989) formula
   I.theta = EIs[,1] %*% t(EIs[,1])
   for(j in 2: n) I.theta = I.theta + EIs[,j] %*% t(EIs[,j])
   sd.theta = rep(NA, m)
   Det = (g-1)*m1+1 
   V.theta = try(solve(I.theta[-Det, -Det]), silent=F)
   if(class(V.theta)[1] != "try-error") sd.theta = sqrt(diag(V.theta))
   if(distr == 'MVN'){
     EST = sd.mvn = NULL
     for(i in 1: g){
      EST = c(EST, c(w[i], mu[,i], Sigma[,,i][vechS], nu[i]))
      sd.mvn = c(sd.mvn, c(sd.theta[(aa[i]+1): aa[(i+1)]], 0))
     }
     sd.theta = sd.mvn[-c(m+(g-1))]
   }
   if(distr == 'MVN'){
     theta.hat = rbind(c(EST[-c((g-1)*(1+p+q1+1)+1)]), sd.theta)
   } else theta.hat = rbind(c(EST[-Det]), sd.theta)
   if(g==1) theta.hat = theta.hat[,-1]
   return(list(theta.hat=theta.hat, I.theta = I.theta, V.theta = V.theta, SD=sd.theta))
}
