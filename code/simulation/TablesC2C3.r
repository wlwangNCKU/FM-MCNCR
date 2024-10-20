##################################################################################
#
#   Filename: TablesC2C3.r
#   Purpose: produce Table C.2 and Table C.3 for SD and SE for parameter estimates 
#            in Experiment 2
#   Input data files: intermediate results (estSD.txt) sorted in the subfolders 
#                    'results/Experiment1/SIM1a/'; 'results/Experiment1/SIM1b/'; 
#                    'results/Experiment1/SIM2a/'; 'results/Experiment1/SIM2b/'; 
#                    'results/Experiment1/SIM3a/'; 'results/Experiment1/SIM3b/'; 
#                    'results/Experiment1/SIM4a/'; 'results/Experiment1/SIM4b/'; 
#   Output data files: results/Experiment1/TableC2.csv
#                      results/Experiment1/TableC3.csv
#
##################################################################################

PATH1a = paste(PATH, 'results/Experiment1/SIM1a/', sep='')
PATH1b = paste(PATH, 'results/Experiment1/SIM1b/', sep='')
PATH2a = paste(PATH, 'results/Experiment1/SIM2a/', sep='')
PATH2b = paste(PATH, 'results/Experiment1/SIM2b/', sep='')
PATH3a = paste(PATH, 'results/Experiment1/SIM3a/', sep='')
PATH3b = paste(PATH, 'results/Experiment1/SIM3b/', sep='')
PATH4a = paste(PATH, 'results/Experiment1/SIM4a/', sep='')
PATH4b = paste(PATH, 'results/Experiment1/SIM4b/', sep='')

# MCN #
Tab1a3 = read.table(paste(PATH1a, 'estSD.txt',sep=""))
Tab2a3 = read.table(paste(PATH2a, 'estSD.txt',sep=""))
Tab3a3 = read.table(paste(PATH3a, 'estSD.txt',sep=""))
Tab4a3 = read.table(paste(PATH4a, 'estSD.txt',sep=""))

# MVN #
Tab1b3 = read.table(paste(PATH1b, 'estSD.txt',sep=""))
Tab2b3 = read.table(paste(PATH2b, 'estSD.txt',sep=""))
Tab3b3 = read.table(paste(PATH3b, 'estSD.txt',sep=""))
Tab4b3 = read.table(paste(PATH4b, 'estSD.txt',sep=""))

name2 = c('Rep', 'rate', 'model',
          paste(c('w','mu1','mu2','mu3','sig11','sig21','sig22',
          'sig31','sig32','sig33','nu','rho'), rep(1:2, each=12), seq=''),
          paste('SD.', c('w','mu1','mu2','mu3','sig11','sig21',
          'sig22','sig31','sig32','sig33','nu','rho'), rep(1:2, each=12), seq=''))
          
colnames(Tab1a3) = colnames(Tab1b3) = colnames(Tab2a3) = colnames(Tab2b3) =
colnames(Tab3a3) = colnames(Tab3b3) = colnames(Tab4a3) = colnames(Tab4b3) = name2

# Table: Parameter estimation
vech.posi=function(dim) cbind(rep(1:dim, 1:dim), unlist(mapply(':', 1, 1:dim)))
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

par.true = NULL
for(i in 1: g) par.true = c(par.true, c(w[i], mu[,i], Sigma[,,i][vech.posi(p)], nu[i], rho[i]))

ESTSE = function(Tab3, rate, par.true)
{
M1 = Tab3[which(Tab3$rate == rate & Tab3$model == 'MCN'), 4:27]
M2 = Tab3[which(Tab3$rate == rate & Tab3$model == 'MVNC'), 4:27]
M3 = Tab3[which(Tab3$rate == rate & Tab3$model == 'MCNC'), 4:27]
rep1 = nrow(M1)

MCNout = rbind(colMeans(M1),
      c(colMeans(M1) - par.true),
      colMeans((M1 - matrix(rep(par.true, rep1), nrow=rep1, byrow=T))^2),
      apply(M1, 2, sd),
      colMeans(Tab3[which(Tab3$rate == rate & Tab3$model == 'MCN'), 28:51]))

MVNCout = rbind(colMeans(M2),
      c(colMeans(M2) - par.true),
      colMeans((M2 - matrix(rep(par.true, rep1), nrow=rep1, byrow=T))^2),
      apply(M2, 2, sd),
      colMeans(Tab3[which(Tab3$rate == rate & Tab3$model == 'MVNC'), 28:51]))

MCNCout = rbind(colMeans(M3),
      c(colMeans(M3) - par.true),
      colMeans((M3 - matrix(rep(par.true, rep1), nrow=rep1, byrow=T))^2),
      apply(M3, 2, sd),
#      colMeans(Tab3[which(Tab3$rate == rate & Tab3$model == 'MCNC'), 28:51]))
      apply(Tab3[which(Tab3$rate == rate & Tab3$model == 'MCNC'), 28:51], 2, mean, trim=0.05))

rownames(MCNout) = rownames(MVNCout) = rownames(MCNCout) = c('EST', 'Bias', 'MSE', 'SD', 'SE')
return(list(MCNout = MCNout, MVNCout = MVNCout, MCNCout=MCNCout))
}

# Table C.2 #
# MCN data with censoring rate=10%
est1a1 = ESTSE(Tab1a3, rate=0.1, par.true=par.true)
est2a1 = ESTSE(Tab2a3, rate=0.1, par.true=par.true)
est3a1 = ESTSE(Tab3a3, rate=0.1, par.true=par.true)
est4a1 = ESTSE(Tab4a3, rate=0.1, par.true=par.true)

# MCN data with censoring rate=20%
est1a2 = ESTSE(Tab1a3, rate=0.2, par.true=par.true)
est2a2 = ESTSE(Tab2a3, rate=0.2, par.true=par.true)
est3a2 = ESTSE(Tab3a3, rate=0.2, par.true=par.true)
est4a2 = ESTSE(Tab4a3, rate=0.2, par.true=par.true)

TabC2 = cbind(
round(cbind(t(est1a1$MCNCout)[,4:5], t(est2a1$MCNCout)[,4:5], t(est3a1$MCNCout)[,4:5], t(est4a1$MCNCout)[,4:5]), 3),
round(cbind(t(est1a2$MCNCout)[,4:5], t(est2a2$MCNCout)[,4:5], t(est3a2$MCNCout)[,4:5], t(est4a2$MCNCout)[,4:5]), 3))
colnames(TabC2) = paste(c(rep(c('C10%-n150', 'C10%-n300', 'C10%-n600', 'C10%-n900',
                    'C20%-n150', 'C20%-n300', 'C20%-n600', 'C20%-n900'), each=2)),
                  c('SD', 'SE'), sep=' ')
TabC2 = TabC2[-13, ]
write.csv(TabC2, paste(PATH, 'results/Experiment1/TableC2.csv',sep=""))

# Table C.3 #
# MVN data with censoring rate=10%
est1b1 = ESTSE(Tab1b3, rate=0.1, par.true=par.true)
est2b1 = ESTSE(Tab2b3, rate=0.1, par.true=par.true)
est3b1 = ESTSE(Tab3b3, rate=0.1, par.true=par.true)
est4b1 = ESTSE(Tab4b3, rate=0.1, par.true=par.true)

# MVN data with censoring rate=20%
est1b2 = ESTSE(Tab1b3, rate=0.2, par.true=par.true)
est2b2 = ESTSE(Tab2b3, rate=0.2, par.true=par.true)
est3b2 = ESTSE(Tab3b3, rate=0.2, par.true=par.true)
est4b2 = ESTSE(Tab4b3, rate=0.2, par.true=par.true)

TabC3 = cbind(
round(cbind(t(est1b1$MVNCout)[,4:5], t(est2b1$MVNCout)[,4:5], t(est3b1$MVNCout)[,4:5], t(est4b1$MVNCout)[,4:5]), 3),
round(cbind(t(est1b2$MVNCout)[,4:5], t(est2b2$MVNCout)[,4:5], t(est3b2$MVNCout)[,4:5], t(est4b2$MVNCout)[,4:5]), 3))
colnames(TabC3) = paste(c(rep(c('C10%-n150', 'C10%-n300', 'C10%-n600', 'C10%-n900',
                    'C20%-n150', 'C20%-n300', 'C20%-n600', 'C20%-n900'), each=2)),
                  c('SD', 'SE'), sep=' ')
TabC3 = TabC3[-11, ]
write.csv(TabC3, paste(PATH, 'results/Experiment1/TableC3.csv',sep=""))

save.image(paste(PATH, 'results/Experiment1/TablesC2C3.RData', sep=''))
