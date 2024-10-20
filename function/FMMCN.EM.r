library(combinat)

# ECM Algorithm for FMMCN-R and FMMN-R Models 
FMMCN.EM = function(Y, g, distr=c('MVN','MCN'), true.clus=NULL, tol=1e-5, max.iter=1000, per=10, mu.true=NULL, SE=FALSE)
{
    require(mvtnorm)
    begin = proc.time()[1]
    n = nrow(Y)
    p = ncol(Y)
    vechS = vech.posi(p)
# initial values
    if(g == 1){ 
      clus = rep(1, n)
    } else{
      if(length(true.clus) == 0 | length(unique(true.clus)) != g){
       kmY = kmeans(Y, g, nstart = 25)
       clus = kmY$cluster
      } else{ clus = true.clus}
    }
    w = rep(1/g, g)
    mu = matrix(NA, p, g)
    Sigma = Sig.inv = array(NA, dim = c(p,p,g))
    for(i in 1:g)
    {
      w[i] = nrow(Y[clus == i,])/n
      mu[,i] = colMeans(Y[clus == i,])
      Sigma[,,i] = cov(Y[clus == i,])
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

    wden = matrix(NA, n, g)
    delta = pdf = w.pdf = matrix(1, n, g)
    cent = array(NA, dim = c(p, n, g))
# old observed-data log-likelihood
    for(i in 1: g){
     cent[,,i] = t(Y) - mu[,i]
     delta[,i] = diag(t(cent[,,i]) %*% Sig.inv[,,i] %*% cent[,,i])
     if(distr=='MCN'){
         for(j in 1: n) w.pdf[j, i] = nu[i] * dmvnorm(t(Y[j, ]), mean=mu[,i], sigma=round(Sigma[,,i]/rho[i],6)) + (1-nu[i]) * dmvnorm(t(Y[j, ]), mean=mu[,i], sigma=round(Sigma[,,i],6)) 
     }
     if(distr=='MVN') wden[, i] = log(w[i]) -.5*p*log(2*pi) -.5*det(Sigma[,,i]) -.5*delta[,i]
     if(distr=='MCN') wden[, i] = log(w[i]) + log(w.pdf[, i])
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
    cat("EM is running for the ", g, "-component FM-MR model with ", distr, "distribution.", "\n")
    cat("Initial log-likelihood = ", loglik.old, "\n")

    w.rho = w1 = tau = matrix(1, n, g)
    tauycent2 = array(NA, dim=c(p, p, n))
    repeat
    {
      iter = iter + 1 
      a.hat = wden / indv.den
      Ni = colSums(a.hat)
      w = Ni / n

#E-step
      for(i in 1: g){
        if(distr=='MCN'){
          for(j in 1: n){
            dcn = dmni(c(Y[j, ]), mu=mu[,i], Sigma=round(Sigma[,,i],6), nu=nu[i], rho=rho[i], distr='MCN')
            w.rho[j, i] = nu[i] * dmvnorm(t(c(Y[j, ])), mean=mu[,i], sigma=round(Sigma[,,i]/rho[i], 6)) / dcn 
            w1[j, i] = (1-nu[i]) * dmvnorm(t(c(Y[j, ])), mean=mu[,i], sigma=round(Sigma[,,i],6)) / dcn
          }
          tau[,i] = rho[i] * w.rho[,i] + w1[,i]
      }} 

# M-step
      y = as.vector(t(Y))
      for(i in 1: g){
        TSig = kronecker(diag(n), Sigma[,,i])
        TSig.inv = kronecker(diag(n), Sig.inv[,,i])
        mu[,i] = colSums(a.hat[,i] * tau[,i] * Y) / sum(a.hat[,i] * tau[,i])
        cent[,,i] = t(Y) - mu[,i]
        for(j in 1: n) tauycent2[,,j] = (a.hat[j,i] * tau[j,i] * cent[, j,i]) %*% t(cent[, j,i])
        Sigma[,,i] = apply(tauycent2, 1:2, sum)/ Ni[i]
        Sigma[,,i] = (Sigma[,,i] + t(Sigma[,,i]))/2
        Sig.inv[,,i] = solve(Sigma[,,i])
        if(distr == 'MCN'){
          nu[i] = sum(a.hat[,i] * w.rho[,i]) / Ni[i]
          d = diag(t(cent[,,i]) %*% Sig.inv[,,i] %*% cent[,,i])
          rho[i] = min(1, ((p*sum(a.hat[,i] * w.rho[,i])) / sum(a.hat[,i] * w.rho[,i] * d))) 
      }}

# new observed-data log-likelihood
     for(i in 1: g){
     cent[,,i] = t(Y) - mu[,i]
     delta[,i] = diag(t(cent[,,i]) %*% Sig.inv[,,i] %*% cent[,,i])
     if(distr=='MCN'){
         for(j in 1: n) w.pdf[j, i] = nu[i] * dmvnorm(t(Y[j, ]), mean=mu[,i], sigma=round(Sigma[,,i]/rho[i],6)) + (1-nu[i]) * dmvnorm(t(Y[j, ]), mean=mu[,i], sigma=round(Sigma[,,i],6)) 
     }
     if(distr=='MVN') wden[, i] = log(w[i]) -.5*p*log(2*pi) -.5*det(Sigma[,,i]) -.5*delta[,i]
     if(distr=='MCN') wden[, i] = log(w[i]) + log(w.pdf[, i])
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
    if(iter%%per == 0){
    if(distr=='MVN') cat("iter", iter, ":", "log-likelihood =", loglik.new, "\t diff =", diff, "\t diff.para =", diff.par, '\t w=', w, '\t mu=', as.vector(mu), "\n")
    if(distr=='MCN') cat("iter", iter, ":", "log-likelihood =", loglik.new, "\t diff =", diff, "\t diff.para =", diff.par, '\t w=', w, '\t nu=', nu, '\t rho=', rho, '\t mu=', as.vector(mu), "\n")
    }
    if(diff < tol | diff.par < tol | iter > max.iter) break
    loglik.old = loglik.new
    old.par = new.par
    }
    end = proc.time()[1]
    cat(paste(rep("-", 50), collapse = ""), "\n")
    run.sec = end - begin
# estimation
    EST = NULL
    if(distr=='MVN'){
     para = list(w=w, mu = mu, Sigma = Sigma)
     for(i in 1: g) EST = c(EST, c(w[i], mu[,i], Sigma[,,i][vechS]))
    }
    if(distr=='MCN'){
     para = list(w=w, mu = mu, Sigma = Sigma, nu = nu, rho = rho)
     for(i in 1: g) EST = c(EST, c(w[i], mu[,i], Sigma[,,i][vechS], nu[i], rho[i]))
    }
# model selection
    m = length(EST) - 1
    aic = 2 * m - 2 * loglik.new
    bic = m * log(n) - 2 * loglik.new
    model.inf = list(m=m, iter = iter, loglik = loglik.new, aic = aic, bic = bic, iter.lnL = iter.lnL, iter.EST = iter.EST, run=run.sec)
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
    xi.hat1 = rowSums(a.hat * w.rho)
    xi.hat2 = rowSums(z * w.rho)
    xi = list(xi1 = xi.hat1, xi2 = xi.hat2)
    cat(paste(rep("-", 50), collapse = ""), "\n")
    if(SE == TRUE){
      IM = I.FMMCN(para, Y=Y, distr=distr, EST=EST, wden=wden)
      return(list(model.inf = model.inf, para = para, EST = EST, IM = IM, pre.cls = pre.cls, xi = xi))
    } else{
      return(list(model.inf = model.inf, para = para, EST = EST, pre.cls = pre.cls, xi = xi))
    }
}

# Louis' information matrix
I.FMMCN = function(para.est, Y, distr=c('MVN','MCN'), EST, wden)
{
  GB = GenzBretz(maxpts = 5e4, abseps = 1e-9, releps = 0)
  n = nrow(Y)
  p = ncol(Y)
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

  indv.den = rowSums(wden)
  a.hat = wden / indv.den
  Ni = colSums(a.hat)
  w.rho = w1 = tau = matrix(1, n, g)
  cent = array(NA, dim = c(p, n, g))
  for(i in 1: g){
    cent[,,i] = t(Y) - mu[,i]
    if(distr=='MCN'){
      for(j in 1: n){
        dcn = dmni(c(Y[j, ]), mu=mu[,i], Sigma=round(Sigma[,,i],6), nu=nu[i], rho=rho[i], distr='MCN')
        w.rho[j, i] = nu[i] * dmvnorm(t(c(Y[j, ])), mean=mu[,i], sigma=round(Sigma[,,i]/rho[i], 6)) / dcn 
        w1[j, i] = (1-nu[i]) * dmvnorm(t(c(Y[j, ])), mean=mu[,i], sigma=round(Sigma[,,i],6)) / dcn
      }
      tau[,i] = rho[i] * w.rho[,i] + w1[,i]
  }} 
# individual expected score vector
  q1 = p*(p+1)/2 
  m = g * (1 + p + q1)
  if(distr == 'MCN') m = m + 2 * g 

  dot.Sig = array(0, dim=c(p, p, q1))
  for(l in 1: q1){
   dotS = matrix(0, p, p)
   dotS[matrix(vechS[l, ], 1)] = dotS[matrix(rev(vechS[l, ]), 1)] = 1
   dot.Sig[,,l] = dotS
  }

  m1 = m / g
  aa = seq(0, m, m1)
  EIs = matrix(NA, nrow=m, ncol=n)
  for(i in 1: g){
     EIs[(aa[i]+1): aa[(i+1)], ][1, ] = a.hat[, i]/w[i]
     for(j in 1: n){
      EIs[(aa[i]+1): aa[(i+1)], ][2:(1+p), j] = (a.hat[j, i] * tau[j, i]) * Sig.inv[,,i] %*% (Y[j,] - mu[,i]) 
      A1 = (tau[j,i] * cent[, j,i]) %*% t(cent[, j,i])
      A2 = (w.rho[j,i] * cent[, j,i]) %*% t(cent[, j,i])
     for(k in 1: q1){
      SSdot = Sig.inv[,,i] %*% dot.Sig[,,k]
      EIs[(aa[i]+1): aa[(i+1)], ][(1+p+k), j] = 0.5 * a.hat[j,i] * (tr(A1 %*% SSdot %*% Sig.inv[,,i]) - tr(SSdot))
     }
     if(distr == 'MCN') EIs[(aa[i]+1): aa[(i+1)], ][(3+p+q1), j] = 0.5 * a.hat[j,i] * (p*w.rho[j,i]/rho[i] - tr(A2 %*% Sig.inv[,,i]))
    }
    if(distr == 'MCN') EIs[(aa[i]+1): aa[(i+1)], ][(2+p+q1), ] = a.hat[,i] * (w.rho[,i] - nu[i]) / (nu[i] * (1-nu[i]))
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

