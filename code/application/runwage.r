#install.packages("sampleSelection")

################################################################################
#
#   Filename: runwage.r
#   Purpose: produce the fitting results for the Mroz87 data
#   Input data files: function/TMmoment.r;
#                     function/FMMVTCRint.EM.r;
#                     function/FMMCNCRint.EM.r;
#                     function/mclust.fn.r
#   Output data files: data/StdWageFit.RData
#   Required R packages: sampleSelection, maxLik, miscTools, combinat, 
#                        cubature, tmvtnorm, gmm, sandwich, Matrix, mvtnorm 
#
################################################################################

source(paste(PATH, 'function/TMmoment.r', sep=''))
source(paste(PATH, 'function/FMMVTCRint.EM.r', sep=''))
source(paste(PATH, 'function/FMMCNCRint.EM.r', sep=''))
source(paste(PATH, 'function/mclust.fn.r', sep=''))

library("sampleSelection")
data(Mroz87)
attach(Mroz87)

Y = cbind(hours, hushrs, wage, huswage)
n = nrow(Y); p = ncol(Y)
u1.std = (rep(0, 4) - colMeans(Y))/apply(Y, 2, sd)
u2.std = (c(40, 40, 2.1, 2.1) - colMeans(Y))/apply(Y, 2, sd)
Yc1 = Yc2 = scale(Y)

cen = matrix(0, nrow=n, ncol=p)
for(k in 1: p){
 Yc1[which(Y[,k] == 0), k] = u1.std[k] 
 Yc2[which(Y[,k] == 0), k] = u2.std[k]
 cen[,k] = as.numeric(Yc2[,k] == u2.std[k])
}
colMeans(cen)

x1 = cbind(kids5+kids618, age, educ)
x2 = cbind(kids5+kids618, husage, huseduc)
X = NULL
for(j in 1:n) X = rbind(X, kronecker(diag(2), rbind(c(t(x1[j,]), 
                        rep(0,ncol(x2))),c(rep(0,ncol(x1)), t(x2[j,])))))
set.seed(28) 

# Fit FM-MCNCR #
repeat{
estCC1 = try(FMMCNCRint.EM(Yc1=Yc1, Yc2=Yc2, X=X, cen=cen, g=1, distr='MCN', tol=1e-3, per=10, SE=TRUE), silent=F)
          if(class(estCC1) != "try-error") break;
}
repeat{
estCC2 = try(FMMCNCRint.EM(Yc1=Yc1, Yc2=Yc2, X=X, cen=cen, g=2, distr='MCN', tol=1e-3, per=10, SE=TRUE), silent=F)
          if(class(estCC2) != "try-error") break;
}
repeat{
estCC3 = try(FMMCNCRint.EM(Yc1=Yc1, Yc2=Yc2, X=X, cen=cen, g=3, distr='MCN', tol=1e-3, per=10, SE=TRUE), silent=F)
          if(class(estCC3) != "try-error") break;
}

# Fit FM-MVNCR #
repeat{
estNC1 = try(FMMCNCRint.EM(Yc1=Yc1, Yc2=Yc2, X=X, cen=cen, g=1, distr='MVN', tol=1e-3, per=10, SE=TRUE), silent=F)
          if(class(estNC1) != "try-error") break;
}
repeat{
estNC2 = try(FMMCNCRint.EM(Yc1=Yc1, Yc2=Yc2, X=X, cen=cen, g=2, distr='MVN', tol=1e-3, per=10, SE=TRUE), silent=F)
          if(class(estNC2) != "try-error") break;
}
repeat{
estNC3 = try(FMMCNCRint.EM(Yc1=Yc1, Yc2=Yc2, X=X, cen=cen, g=3, distr='MVN', tol=1e-3, per=10, SE=TRUE), silent=F)
          if(class(estNC3) != "try-error") break;
}

# Fit FM-MVTCR #
repeat{
estTC1 = try(FMMVTCRint.EM(Yc1=Yc1, Yc2=Yc2, X=X, cen=cen, g=1, distr='MVT', tol=1e-3, per=10, SE=TRUE), silent=F)
          if(class(estTC1) != "try-error") break;
}
repeat{
estTC2 = try(FMMVTCRint.EM(Yc1=Yc1, Yc2=Yc2, X=X, cen=cen, g=2, distr='MVT', tol=1e-3, per=10, SE=TRUE), silent=F)
          if(class(estTC2) != "try-error") break;
}
repeat{
estTC3 = try(FMMVTCRint.EM(Yc1=Yc1, Yc2=Yc2, X=X, cen=cen, g=3, distr='MVT', tol=1e-3, per=10, SE=TRUE), silent=F)
          if(class(estTC3) != "try-error") break;
}

save.image(paste(PATH, 'data/StdWageFit.RData', sep=""))
