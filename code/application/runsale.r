################################################################################
#
#   Filename: runsale.r
#   Purpose: produce the fitting results for the Wholesale data
#   Input data files: data/Wholesale.RData;
#                     function/TMmoment.r;
#                     function/FMMVTCint.EM.r;
#                     function/FMMCNCint.EM.r;
#                     function/mclust.fn.r
#   Output data files: data/SaleFit.RData
#   Required R packages: mvtnorm, tmvtnorm, cubature, 
#                        gmm, combinat, sandwich, Matrix
#
################################################################################

source(paste(PATH, 'function/TMmoment.r', sep=''))
source(paste(PATH, 'function/FMMVTCint.EM.r', sep=''))
source(paste(PATH, 'function/FMMCNCint.EM.r', sep=''))
source(paste(PATH, 'function/mclust.fn.r', sep=''))

load(paste(PATH, 'data/Wholesale.RData', sep=''))
dataset = Wholesale[,4:7]
n = nrow(dataset); n
p = ncol(dataset); p
cls = Wholesale$Channel
Y = scale(dataset)

u1.std = (rep(0, p) - colMeans(dataset))/apply(dataset, 2, sd)
u2.std = (rep(100, p) - colMeans(dataset))/apply(dataset, 2, sd)
Yc1 = Yc2 = Y
cen = matrix(0, nrow=n, ncol=p)
for(k in 1: p){
 Yc1[which(Y[,k] <= u2.std[k]), k] = u1.std[k] 
 Yc2[which(Y[,k] <= u2.std[k]), k] = u2.std[k]
 cen[,k] = as.numeric(Yc2[,k] == u2.std[k])
}
colMeans(cen)
cls2 = Wholesale$Channel
set.seed(7)

# Model fitting #
repeat{ 
estNC2 = try(FMMCNCint.EM(Yc1=Yc1, Yc2=Yc2, cen=cen, g=2, distr='MVN', tol=1e-5, per=10, SE=TRUE), silent=F)
          if(class(estNC2) != "try-error") break;
}
CCRnc2 = 1 - classError(cls2, estNC2$pre.cls$post.cls)$errorRate
ARInc2 = adjustedRandIndex(cls2, estNC2$pre.cls$post.cls)
cat('CCR =', CCRnc2, ',\t ARI =', ARInc2, '\n')

repeat{
estCC2 = try(FMMCNCint.EM(Yc1=Yc1, Yc2=Yc2, cen=cen, g=2, distr='MCN', tol=1e-5, per=10, SE=TRUE), silent=F)
          if(class(estCC2) != "try-error") break;
}
CCRcc2 = 1 - classError(cls2, estCC2$pre.cls$post.cls)$errorRate
ARIcc2 = adjustedRandIndex(cls2, estCC2$pre.cls$post.cls)
cat('CCR =', CCRcc2, ',\t ARI =', ARIcc2, '\n')

repeat{
estTC2 = try(FMMVTCint.EM(Yc1=Yc1, Yc2=Yc2, cen=cen, g=2, distr='MVT', tol=1e-5, per=10, SE=TRUE), silent=F)
          if(class(estTC2) != "try-error") break;
}
CCRtc2 = 1 - classError(cls2, estTC2$pre.cls$post.cls)$errorRate
ARItc2 = adjustedRandIndex(cls2, estTC2$pre.cls$post.cls)
cat('CCR =', CCRtc2, ',\t ARI =', ARItc2, '\n')

save.image(paste(PATH, 'data/SaleFit.RData', sep=""))
