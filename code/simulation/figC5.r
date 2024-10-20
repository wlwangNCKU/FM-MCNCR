##################################################################################
#
#   Filename: figC5.r
#   Purpose: produce Figure C.5 for MSE for MLE across sample size in Experiment 2
#   Input data files: results/Experiment2/TableC5.RData 
#   Output data files: results/Experiment2/figC5.eps
#
##################################################################################

load(paste(PATH, 'results/Experiment2/TableC5.RData', sep=''))

Size = c(150, 300, 600, 900)
# crate = 10%
MSE1 = cbind(out11$MCNCout['MSE', ], out21$MCNCout['MSE', ], out31$MCNCout['MSE', ], out41$MCNCout['MSE', ])

# crate = 20%
MSE2 = cbind(out12$MCNCout['MSE', ], out22$MCNCout['MSE', ], out32$MCNCout['MSE', ], out42$MCNCout['MSE', ])

colnames(MSE1) = colnames(MSE2) = paste('n=',Size,sep='')
m = nrow(MSE1)

PATH = paste(getwd(),"/Data_and_Code_FMMCNCR/",sep="")
postscript(paste(PATH, 'results/Experiment2/figC5.eps', sep=''), width=40, height=6)
layout(matrix(c(0,rep(12,8), 1,13:20, 2,21:28, 3,29:36, 0,4:11), 5, 9, byrow=T), c(1.5,rep(5,8)), c(1,rep(4.5, 3),1))
par(mar=c(0,0,0,0))
for(i in 1: 3){
plot(c(0, 1), c(0, 3), type='n', xaxt='n', yaxt='n', xlab='', ylab='', bty='n')
mtext('MSE', 2, line=-1.5, cex=1)
}
par(mar=c(0.5,0,0,0))
for(i in 1: 8){
plot(c(0, 1), c(0, 3), type='n', xaxt='n', yaxt='n', xlab='', ylab='', bty='n')
mtext('      Sample Size', 1, line=-2.5, cex=0.75)
}
par(mar=c(0,0,0.5,0))
plot(c(0, 1), c(0, 3), type='n', xaxt='n', yaxt='n', xlab='', ylab='', bty='n')
legend("topright", c('10% Censoring', '20% Censoring'), pch=c(1,19), lty=c(2,1), bty='n', ncol=2, cex=1.5, lwd=1.5)

theta.name = c(expression(w[1]), expression(beta[10]), expression(beta[11]), expression(sigma[111]), expression(sigma[121]), expression(sigma[122]), expression(nu[1]), expression(rho[1]),
               expression(w[2]), expression(beta[20]), expression(beta[21]), expression(sigma[211]), expression(sigma[221]), expression(sigma[222]), expression(nu[2]), expression(rho[2]),
               expression(w[3]), expression(beta[30]), expression(beta[31]), expression(sigma[311]), expression(sigma[321]), expression(sigma[322]), expression(nu[3]), expression(rho[3]))
par(mar=c(2.5,2.5,1.5,0))
for(r in 1: m){
 ylim1 = range(c(MSE1[r, 1:4], MSE2[r, 1:4]))
 plot(1:4, MSE1[r, 1:4], type='n', ylab='', xlab='', xlim=c(0.75, 4.25), ylim=ylim1, las=1, xaxt='n', cex.axis=0.75)
 axis(1, 1:4, labels=Size, cex.axis=0.75)
 mtext(theta.name[r], 3, line=0, cex=1, font.main=3)
 lines(1:4, MSE1[r, 1:4], lty=2, lwd=1.2)
 points(1:4, MSE1[r, 1:4], pch=1, cex=1.5)
 lines(1:4, MSE2[r, 1:4], lty=1, lwd=1.2)
 points(1:4, MSE2[r, 1:4], pch=19, cex=1.2)
}
dev.off()
