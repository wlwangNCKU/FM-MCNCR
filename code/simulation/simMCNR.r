################################################################################
#
#   Filename: simMCNR.r
#   Purpose: produce the intermediate results for Experiment 2
#   Input data files: function/TMmoment.r;
#                     function/FMMCNCR.EM.r;
#                     function/FMMCNR.EM.r;
#                     function/mclust.fn.r
#   Output data files: results/Experiment2/.../modelfit.txt;
#                      results/Experiment2/.../cls.txt;
#                      results/Experiment2/.../estSD.txt;
#                      results/Experiment2/.../MSEest.txt;
#                      results/Experiment2/.../MSEyfit.txt;
#                      results/Experiment2/.../MSEyc.txt
#   Required R packages: sampleSelection, maxLik, miscTools, combinat, Matrix,       
#                        cubature, mvtnorm, tmvtnorm, gmm, sandwich        
#
################################################################################

source(paste(PATH, 'function/TMmoment.r', sep=''))
source(paste(PATH, 'function/FMMCNCR.EM.r', sep=''))
source(paste(PATH, 'function/FMMCNR.EM.r', sep=''))
source(paste(PATH, 'function/mclust.fn.r', sep=''))

Rep = 1
total.rep = 100
ti = c(1, 2)
p = length(ti)
Xi = cbind(1, ti)
X = cbind(1, rep(ti, n))
g = 3
w = rep(1/g, g)
Beta = cbind(c(2, -3), c(2, 2), c(-2, 3))
Sigma = array(NA, dim=c(p, p, g))
Cor.g = c(0.1, 0.3, 0.5)
for(i in 1: g) diag(Sigma[,,i]) = rep(1,p)
for(i in 1: g){
 for(j in 1:(p-1))for(k in (j+1):p) Sigma[j,k,i] = Sigma[k,j,i] = Cor.g[i] * sqrt(Sigma[j,j,i] * Sigma[k,k,i]) 
} 
nu = rep(0.25, g)
rho = rep(0.2, g)
m = g + length(Beta) + g*(p*(p+1)/2 + 2)

para = list(w=w, Beta=Beta, Sigma=Sigma, nu=nu, rho=rho)
par.true = NULL
for(i in 1: g) par.true = c(par.true, c(w[i], Beta[,i], Sigma[,,i][vech.posi(p)], nu[i], rho[i]))
while(Rep <= total.rep)
{
cat(paste(c(rep('=', 25), rep(' ', 5), 'The ', Rep, 'th replication: n = ', n, rep(' ', 5), rep('=', 25)), sep = '', collapse = ''), '\n')
ygen = rfmmcn.reg(n=n, para=para, Xi=Xi, distr='MCN')
cls = ygen$cls

# censoring rate = 20% # 
crate2 = 0.2                               
cat(paste(c(rep('-', 20), rep(' ', 5), Rep, ': n = ', n, ', the censoring rate = ', crate2*100, '%', rep(' ', 5), rep('-', 20)), sep = '', collapse = ''), '\n')  
Cdata2 = rcen(ygen$Y, cen.rate=crate2)
Yc2 = Cdata2$Y; cen2 = Cdata2$CEN; 
estCC23 = try(FMMCNCR.EM(Yc=Yc2, X=X, cen=cen2, g=3, distr='MCN', true.clus=cls, tol=1e-5, Beta.true=para$Beta, SE=TRUE), silent=F)
          if(class(estCC23) == "try-error") next;
estC23 =  try(FMMCNR.EM(Y=Yc2, X=X, g=3, distr='MCN', true.clus=cls, tol=1e-5, Beta.true=para$Beta, SE=TRUE), silent=F)
          if(class(estC23) == "try-error") next;

estCC24 = try(FMMCNCR.EM(Yc=Yc2, X=X, cen=cen2, g=4, distr='MCN', tol=1e-5), silent=F)
          if(class(estCC24) == "try-error") next;
estC24 =  try(FMMCNR.EM(Y=Yc2, X=X, g=4, distr='MCN', tol=1e-5), silent=F)
          if(class(estC24) == "try-error") next;

estCC22 = try(FMMCNCR.EM(Yc=Yc2, X=X, cen=cen2, g=2, distr='MCN', tol=1e-5), silent=F)
          if(class(estCC22) == "try-error") next;
estC22 =  try(FMMCNR.EM(Y=Yc2, X=X, g=2, distr='MCN', tol=1e-5), silent=F)
          if(class(estC22) == "try-error") next;

estCC21 = try(FMMCNCR.EM(Yc=Yc2, X=X, cen=cen2, g=1, distr='MCN', tol=1e-5), silent=F)
          if(class(estCC21) == "try-error") next;
estC21 =  try(FMMCNR.EM(Y=Yc2, X=X, g=1, distr='MCN', tol=1e-5), silent=F)
          if(class(estC21) == "try-error") next;

# censoring rate = 10% # 
crate1 = 0.1
cat(paste(c(rep('-', 20), rep(' ', 5), Rep, ': n = ', n, ', the censoring rate = ', crate1*100, '%', rep(' ', 5), rep('-', 20)), sep = '', collapse = ''), '\n')    
Cdata1 = rcen(ygen$Y, cen.rate=crate1)
Yc1 = Cdata1$Y; cen1 = Cdata1$CEN; 
estCC13 = try(FMMCNCR.EM(Yc=Yc1, X=X, cen=cen1, g=3, distr='MCN', true.clus=cls, tol=1e-5, Beta.true=para$Beta, SE=TRUE), silent=F)
          if(class(estCC13) == "try-error") next;
estC13 =  try(FMMCNR.EM(Y=Yc1, X=X, g=3, distr='MCN', true.clus=cls, tol=1e-5, Beta.true=para$Beta, SE=TRUE), silent=F)
          if(class(estC13) == "try-error") next;

estCC14 = try(FMMCNCR.EM(Yc=Yc1, X=X, cen=cen1, g=4, distr='MCN', tol=1e-5), silent=F)
          if(class(estCC14) == "try-error") next;
estC14 =  try(FMMCNR.EM(Y=Yc1, X=X, g=4, distr='MCN', tol=1e-5), silent=F)
          if(class(estC14) == "try-error") next;

estCC12 = try(FMMCNCR.EM(Yc=Yc1, X=X, cen=cen1, g=2, distr='MCN', tol=1e-5), silent=F)
          if(class(estCC12) == "try-error") next;
estC12 =  try(FMMCNR.EM(Y=Yc1, X=X, g=2, distr='MCN', tol=1e-5), silent=F)
          if(class(estC12) == "try-error") next;

estCC11 = try(FMMCNCR.EM(Yc=Yc1, X=X, cen=cen1, g=1, distr='MCN', tol=1e-5), silent=F)
          if(class(estCC11) == "try-error") next;
estC11 =  try(FMMCNR.EM(Y=Yc1, X=X, g=1, distr='MCN', tol=1e-5), silent=F)
          if(class(estC11) == "try-error") next;

# model selection
Tab1 = rbind(c(Rep, crate1, 'loglik', estC11$model.inf$loglik, estC12$model.inf$loglik, estC13$model.inf$loglik, estC14$model.inf$loglik, estCC11$model.inf$loglik, estCC12$model.inf$loglik, estCC13$model.inf$loglik, estCC14$model.inf$loglik),
             c(Rep, crate2, 'loglik', estC21$model.inf$loglik, estC22$model.inf$loglik, estC23$model.inf$loglik, estC24$model.inf$loglik, estCC21$model.inf$loglik, estCC22$model.inf$loglik, estCC23$model.inf$loglik, estCC24$model.inf$loglik),
             c(Rep, crate1, 'aic', estC11$model.inf$aic, estC12$model.inf$aic, estC13$model.inf$aic, estC14$model.inf$aic, estCC11$model.inf$aic, estCC12$model.inf$aic, estCC13$model.inf$aic, estCC14$model.inf$aic),
             c(Rep, crate2, 'aic', estC21$model.inf$aic, estC22$model.inf$aic, estC23$model.inf$aic, estC24$model.inf$aic, estCC21$model.inf$aic, estCC22$model.inf$aic, estCC23$model.inf$aic, estCC24$model.inf$aic),
             c(Rep, crate1, 'bic', estC11$model.inf$bic, estC12$model.inf$bic, estC13$model.inf$bic, estC14$model.inf$bic, estCC11$model.inf$bic, estCC12$model.inf$bic, estCC13$model.inf$bic, estCC14$model.inf$bic),
             c(Rep, crate2, 'bic', estC21$model.inf$bic, estC22$model.inf$bic, estC23$model.inf$bic, estC24$model.inf$bic, estCC21$model.inf$bic, estCC22$model.inf$bic, estCC23$model.inf$bic, estCC24$model.inf$bic))

# classification
Tab2 = rbind(c(Rep, crate1, 'CCR', estC13$pre.cls$CCR, estCC13$pre.cls$CCR),
             c(Rep, crate2, 'CCR', estC23$pre.cls$CCR, estCC23$pre.cls$CCR),
             c(Rep, crate1, 'ARI', estC13$pre.cls$ARI, estCC13$pre.cls$ARI),
             c(Rep, crate2, 'ARI', estC23$pre.cls$ARI, estCC23$pre.cls$ARI))

# parameter estimation
Tab3 = rbind(c(Rep, crate1, 'MCN',  as.vector(t(estC13$IM$theta))), 
             c(Rep, crate1, 'MCNC', as.vector(t(estCC13$IM$theta))), 
             c(Rep, crate2, 'MCN',  as.vector(t(estC23$IM$theta))), 
             c(Rep, crate2, 'MCNC', as.vector(t(estCC23$IM$theta))))
              
Tab4 = rbind(c(Rep, crate1, mse(estC13$IM$theta[1,], par.true), mse(estCC13$IM$theta[1,], par.true)), 
             c(Rep, crate2, mse(estC23$IM$theta[1,], par.true), mse(estCC23$IM$theta[1,], par.true)))

# fitted responses 
y = as.vector(t(ygen$Y))
Tab5 = rbind(c(Rep, crate1, mse(estC13$yfit$yfit1, y), mse(estCC13$yfit$yfit1, y), 
                            mse(estC13$yfit$yfit2, y), mse(estCC13$yfit$yfit2, y)), 
             c(Rep, crate2, mse(estC23$yfit$yfit1, y), mse(estCC23$yfit$yfit1, y),
                            mse(estC23$yfit$yfit2, y), mse(estCC23$yfit$yfit2, y)))

# recovery of censored responses 
yc1 = as.vector(t(ygen$Y))[which(as.vector(t(cen1)) == 1)]
yc2 = as.vector(t(ygen$Y))[which(as.vector(t(cen2)) == 1)]
Tab6 = rbind(c(Rep, crate1, mse(estCC13$ychat$ychat1, yc1), mse(estCC13$ychat$ychat2, yc1)),
             c(Rep, crate2, mse(estCC23$ychat$ychat1, yc2), mse(estCC23$ychat$ychat2, yc2)))

# Explore outputs to folder
write.table(Tab1, paste(WD.PATH, 'modelfit.txt',sep=""), append=T, row.names = F, col.names = F)
write.table(Tab2, paste(WD.PATH, 'cls.txt',sep=""), append=T, row.names = F, col.names = F)
write.table(Tab3, paste(WD.PATH, 'estSD.txt',sep=""), append=T, row.names = F, col.names = F)
write.table(Tab4, paste(WD.PATH, 'MSEest.txt',sep=""), append=T, row.names = F, col.names = F)
write.table(Tab5, paste(WD.PATH, 'MSEyfit.txt',sep=""), append=T, row.names = F, col.names = F)
write.table(Tab6, paste(WD.PATH, 'MSEyc.txt',sep=""), append=T, row.names = F, col.names = F)

Rep = Rep + 1
}
 