################################################################################
#
#   Filename: simMCN.r
#   Purpose: produce the intermediate results for Experiment 1
#   Input data files: function/TMmoment.r;
#                     function/FMMCNC.EM.r;
#                     function/FMMCN.EM.r;
#                     function/mclust.fn.r
#   Output data files: results/Experiment1/.../modelfit.txt;
#                      results/Experiment1/.../cls.txt;
#                      results/Experiment1/.../estSD.txt;
#                      results/Experiment1/.../MSEest.txt;
#                      results/Experiment1/.../MSEyc.txt
#   Required R packages: sampleSelection, maxLik, miscTools, combinat, Matrix,       
#                        cubature, mvtnorm, tmvtnorm, gmm, sandwich        
#
################################################################################

source(paste(PATH, 'function/TMmoment.r', sep=''))
source(paste(PATH, 'function/FMMCNC.EM.r', sep=''))
source(paste(PATH, 'function/FMMCN.EM.r', sep=''))
source(paste(PATH, 'function/mclust.fn.r', sep=''))
 
Rep = 1 
total.rep = 100

g = 2; p = 3
w = rep(1/g, g)
mu = cbind(c(3, 2, 1), c(6, 1, 3))
Sigma = array(NA, dim=c(p, p, g))
Cor.g = c(0.25, 0.5)
for(i in 1: g) diag(Sigma[,,i]) = rep(1, p)
for(i in 1: g){
 for(j in 1:(p-1))for(k in (j+1):p) Sigma[j,k,i] = Sigma[k,j,i] = Cor.g[i] * sqrt(Sigma[j,j,i] * Sigma[k,k,i]) 
} 
nu = rep(0.25, g)
rho = c(0.1, 0.2)
m = g*(1 + p + p*(p+1)/2 + 2)

para = list(w=w, mu=mu, Sigma=Sigma, nu=nu, rho=rho)
par.true = NULL
for(i in 1: g) par.true = c(par.true, c(w[i], mu[,i], Sigma[,,i][vech.posi(p)], nu[i], rho[i]))
while(Rep <= total.rep)
{
cat(paste(c(rep('=', 25), rep(' ', 5), 'The ', Rep, 'th replication: n = ', n, ' sample from ', DIS, ' distribution.', rep(' ', 5), rep('=', 25)), sep = '', collapse = ''), '\n')
ygen = rfmmcn(n=n, para=para, distr=DIS)
cls = ygen$cls

# censoring rate = 20% # 
crate2 = 0.2                               
cat(paste(c(rep('-', 20), rep(' ', 5), Rep, ': n = ', n, ', the censoring rate = ', crate2*100, '%', rep(' ', 5), rep('-', 20)), sep = '', collapse = ''), '\n')  
Cdata2 = rcen(ygen$Y, cen.rate=crate2)
Yc2 = Cdata2$Y; cen2 = Cdata2$CEN; 

estC2 = try(FMMCN.EM(Y=Yc2, g=2, distr='MCN', true.clus=cls, tol=1e-5, mu.true=para$mu, SE=TRUE), silent=F)
        if(class(estC2) == "try-error") next;
estCC2 = try(FMMCNC.EM(Yc=Yc2, cen=cen2, g=2, distr='MCN', true.clus=cls, tol=1e-5, mu.true=para$mu, SE=TRUE), silent=F)
        if(class(estCC2) == "try-error") next;
estNC2 = try(FMMCNC.EM(Yc=Yc2, cen=cen2, g=2, distr='MVN', true.clus=cls, tol=1e-5, mu.true=para$mu, SE=TRUE), silent=F)
        if(class(estNC2) == "try-error") next;

# censoring rate = 10% # 
crate1 = 0.1
cat(paste(c(rep('-', 20), rep(' ', 5), Rep, ': n = ', n, ', the censoring rate = ', crate1*100, '%', rep(' ', 5), rep('-', 20)), sep = '', collapse = ''), '\n')    
Cdata1 = rcen(ygen$Y, cen.rate=crate1)
Yc1 = Cdata1$Y; cen1 = Cdata1$CEN;

estC1 = try(FMMCN.EM(Y=Yc1, g=2, distr='MCN', true.clus=cls, tol=1e-5, mu.true=para$mu, SE=TRUE), silent=F)
        if(class(estC1) == "try-error") next;
estCC1 = try(FMMCNC.EM(Yc=Yc1, cen=cen1, g=2, distr='MCN', true.clus=cls, tol=1e-5, mu.true=para$mu, SE=TRUE), silent=F)
         if(class(estCC1) == "try-error") next;
estNC1 = try(FMMCNC.EM(Yc=Yc1, cen=cen1, g=2, distr='MVN', true.clus=cls, tol=1e-5, mu.true=para$mu, SE=TRUE), silent=F)
        if(class(estNC1) == "try-error") next;

# model selection
Tab1 = rbind(c(Rep, crate1, 'loglik', estC1$model.inf$loglik, estNC1$model.inf$loglik, estCC1$model.inf$loglik),
             c(Rep, crate2, 'loglik', estC2$model.inf$loglik, estNC2$model.inf$loglik, estCC2$model.inf$loglik),
             c(Rep, crate1, 'aic', estC1$model.inf$aic, estNC1$model.inf$aic, estCC1$model.inf$aic),
             c(Rep, crate2, 'aic', estC2$model.inf$aic, estNC2$model.inf$aic, estCC2$model.inf$aic),
             c(Rep, crate1, 'bic', estC1$model.inf$bic, estNC1$model.inf$bic, estCC1$model.inf$bic),
             c(Rep, crate2, 'bic', estC2$model.inf$bic, estNC2$model.inf$bic, estCC2$model.inf$bic))

# classification
Tab2 = rbind(c(Rep, crate1, 'CCR', estC1$pre.cls$CCR, estNC1$pre.cls$CCR, estCC1$pre.cls$CCR),
             c(Rep, crate2, 'CCR', estC2$pre.cls$CCR, estNC2$pre.cls$CCR, estCC2$pre.cls$CCR),
             c(Rep, crate1, 'ARI', estC1$pre.cls$ARI, estNC1$pre.cls$ARI, estCC1$pre.cls$ARI),
             c(Rep, crate2, 'ARI', estC2$pre.cls$ARI, estNC2$pre.cls$ARI, estCC2$pre.cls$ARI))

# parameter estimation
Tab3 = rbind(c(Rep, crate1, 'MCN',  as.vector(t(estC1$IM$theta))), 
             c(Rep, crate1, 'MVNC', as.vector(t(estNC1$IM$theta))), 
             c(Rep, crate1, 'MCNC', as.vector(t(estCC1$IM$theta))), 
             c(Rep, crate2, 'MCN',  as.vector(t(estC2$IM$theta))), 
             c(Rep, crate2, 'MVNC', as.vector(t(estNC2$IM$theta))), 
             c(Rep, crate2, 'MCNC', as.vector(t(estCC2$IM$theta))))
              
Tab4 = rbind(c(Rep, crate1, mse(estC1$IM$theta[1,], par.true), mse(estNC1$IM$theta[1,], par.true), mse(estCC1$IM$theta[1,], par.true)), 
             c(Rep, crate2, mse(estC2$IM$theta[1,], par.true), mse(estNC2$IM$theta[1,], par.true), mse(estCC2$IM$theta[1,], par.true)))

# recovery of censored responses 
yc1 = as.vector(t(ygen$Y))[which(as.vector(t(cen1)) == 1)]
yc2 = as.vector(t(ygen$Y))[which(as.vector(t(cen2)) == 1)]
Tab5 = rbind(c(Rep, crate1, mse(estNC1$ychat$ychat1, yc1), mse(estCC1$ychat$ychat1, yc1), mse(estNC1$ychat$ychat2, yc1), mse(estCC1$ychat$ychat2, yc1)),
             c(Rep, crate2, mse(estNC2$ychat$ychat1, yc2), mse(estCC2$ychat$ychat1, yc2), mse(estNC2$ychat$ychat2, yc2), mse(estCC2$ychat$ychat2, yc2)))

# Explore outputs to folder
write.table(Tab1, paste(WD.PATH, 'modelfit.txt',sep=""), append=T, row.names = F, col.names = F)
write.table(Tab2, paste(WD.PATH, 'cls.txt',sep=""), append=T, row.names = F, col.names = F)
write.table(Tab3, paste(WD.PATH, 'estSD.txt',sep=""), append=T, row.names = F, col.names = F)
write.table(Tab4, paste(WD.PATH, 'MSEest.txt',sep=""), append=T, row.names = F, col.names = F)
write.table(Tab5, paste(WD.PATH, 'MSEyc.txt',sep=""), append=T, row.names = F, col.names = F)

Rep = Rep + 1
}
 