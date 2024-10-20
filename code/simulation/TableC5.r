##################################################################################
#
#   Filename: TableC5.r
#   Purpose: produce Table C.4 for SD and SE of parameter estimates in Experiment 2
#   Input data files: intermediate results (estSD.txt) sorted in the subfolders 
#                     'results/Experiment2/SIM1/'; 'results/Experiment2/SIM2/'; 
#                     'results/Experiment2/SIM3/'; 'results/Experiment2/SIM4/'; 
#   Output data files: results/Experiment2/TableC5.csv
#
##################################################################################

PATH1 = paste(PATH, 'results/Experiment2/SIM1/', sep='')
PATH2 = paste(PATH, 'results/Experiment2/SIM2/', sep='')
PATH3 = paste(PATH, 'results/Experiment2/SIM3/', sep='')
PATH4 = paste(PATH, 'results/Experiment2/SIM4/', sep='')

### Read outputs ###
Tab13 = read.table(paste(PATH1, 'estSD.txt',sep=""))
Tab23 = read.table(paste(PATH2, 'estSD.txt',sep=""))
Tab33 = read.table(paste(PATH3, 'estSD.txt',sep=""))
Tab43 = read.table(paste(PATH4, 'estSD.txt',sep=""))

name3 = c('Rep', 'rate', 'model',
          paste(c('w','beta1','beta2','sig11','sig21','sig22','nu','rho'), rep(1:3, each=8), seq=''),
          paste('SD.', c('w','beta1','beta2','sig11','sig21','sig22','nu','rho'), rep(1:3, each=8), seq=''))
names(Tab13) = names(Tab23) = names(Tab33) = names(Tab43) = name3

g = 3; p = 2
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
vech.posi=function(dim) cbind(rep(1:dim, 1:dim), unlist(mapply(':', 1, 1:dim)))
par.true = NULL
for(i in 1: g) par.true = c(par.true, c(w[i], Beta[,i], Sigma[,,i][vech.posi(p)], nu[i], rho[i]))

ESTSE = function(Tab3, crate, par.true)
{
M1 = Tab3[which(Tab3$rate == crate & Tab3$model == 'MCN'), 4:27]
M2 = Tab3[which(Tab3$rate == crate & Tab3$model == 'MCNC'), 4:27]
rep1 = nrow(M1)

MCNout = rbind(colMeans(M1),
      c(colMeans(M1) - par.true),
      colMeans((M1 - matrix(rep(par.true, rep1), nrow=rep1, byrow=T))^2),
      apply(M1, 2, sd),
      apply(Tab3[which(Tab3$rate == crate & Tab3$model == 'MCN'), 28:51], 2, mean, trim=0.9))

MCNCout = rbind(colMeans(M2),
      c(colMeans(M2) - par.true),
      colMeans((M2 - matrix(rep(par.true, rep1), nrow=rep1, byrow=T))^2),
      apply(M2, 2, sd),
      apply(Tab3[which(Tab3$rate == crate & Tab3$model == 'MCNC'), 28:51], 2, mean, trim=0.9))

rownames(MCNout) = rownames(MCNCout) = c('EST', 'Bias', 'MSE', 'SD', 'SE')
return(list(MCNout = MCNout, MCNCout=MCNCout))
}

out11 = ESTSE(Tab13, crate=0.1, par.true)
out21 = ESTSE(Tab23, crate=0.1, par.true)
out31 = ESTSE(Tab33, crate=0.1, par.true)
out41 = ESTSE(Tab43, crate=0.1, par.true)

out12 = ESTSE(Tab13, crate=0.2, par.true)
out22 = ESTSE(Tab23, crate=0.2, par.true)
out32 = ESTSE(Tab33, crate=0.2, par.true)
out42 = ESTSE(Tab43, crate=0.2, par.true)

TabC5 = t(rbind(out11$MCNCout[4:5, ], out21$MCNCout[4:5, ], 
                out31$MCNCout[4:5, ], out41$MCNCout[4:5, ],
                out12$MCNCout[4:5, ], out22$MCNCout[4:5, ], 
                out32$MCNCout[4:5, ], out42$MCNCout[4:5, ]))
colnames(TabC5) = paste(c(rep(c('C10%-n150', 'C10%-n300', 'C10%-n600', 'C10%-n900',
                    'C20%-n150', 'C20%-n300', 'C20%-n600', 'C20%-n900'), each=2)),
                    c('SD', 'SE'), sep=' ')
TabC5 = TabC5[-17, ]
print(round(TabC5, 3)) 

write.csv(TabC5, paste(PATH, 'results/Experiment2/TableC5.csv',sep=""))
save.image(paste(PATH, 'results/Experiment2/TableC5.RData', sep=''))
