##################################################################################
#
#   Filename: figC2.r
#   Purpose: produce Figure C.3a and C.3b for MSE of MLE of each parameters 
#                    across sample size in Experiment 1
#   Input data files: results/Experiment1/TablesC2C3.RData
#                     intermediate results (cls.txt) sorted in the subfolders 
#                    'results/Experiment1/SIM1a/'; 'results/Experiment1/SIM1b/'; 
#                    'results/Experiment1/SIM2a/'; 'results/Experiment1/SIM2b/'; 
#                    'results/Experiment1/SIM3a/'; 'results/Experiment1/SIM3b/'; 
#                    'results/Experiment1/SIM4a/'; 'results/Experiment1/SIM4b/'; 
#   Output data files: results/Experiment1/figC2.eps
#
##################################################################################

load(paste(PATH, 'results/Experiment1/TablesC2C3.RData', sep=''))

Size = c(150, 300, 600, 900)
# MCN, crate = 10%
MSE1 = cbind(est1a1$MCNout['MSE', ], est2a1$MCNout['MSE', ], est3a1$MCNout['MSE', ], est4a1$MCNout['MSE', ],
         est1a1$MVNCout['MSE', ], est2a1$MVNCout['MSE', ], est3a1$MVNCout['MSE', ], est4a1$MVNCout['MSE', ],
         est1a1$MCNCout['MSE', ], est2a1$MCNCout['MSE', ], est3a1$MCNCout['MSE', ], est4a1$MCNCout['MSE', ])

# MCN, crate = 20%
MSE2 = cbind(est1a2$MCNout['MSE', ], est2a2$MCNout['MSE', ], est3a2$MCNout['MSE', ], est4a2$MCNout['MSE', ],
         est1a2$MVNCout['MSE', ], est2a2$MVNCout['MSE', ], est3a2$MVNCout['MSE', ], est4a2$MVNCout['MSE', ],
         est1a2$MCNCout['MSE', ], est2a2$MCNCout['MSE', ], est3a2$MCNCout['MSE', ], est4a2$MCNCout['MSE', ])

# MVN, crate = 10%
MSE3 = cbind(est1b1$MCNout['MSE', ], est2b1$MCNout['MSE', ], est3b1$MCNout['MSE', ], est4b1$MCNout['MSE', ],
         est1b1$MVNCout['MSE', ], est2b1$MVNCout['MSE', ], est3b1$MVNCout['MSE', ], est4b1$MVNCout['MSE', ],
         est1b1$MCNCout['MSE', ], est2b1$MCNCout['MSE', ], est3b1$MCNCout['MSE', ], est4b1$MCNCout['MSE', ])

# MVN, crate = 20%
MSE4 = cbind(est1b2$MCNout['MSE', ], est2b2$MCNout['MSE', ], est3b2$MCNout['MSE', ], est4b2$MCNout['MSE', ],
         est1b2$MVNCout['MSE', ], est2b2$MVNCout['MSE', ], est3b2$MVNCout['MSE', ], est4b2$MVNCout['MSE', ],
         est1b2$MCNCout['MSE', ], est2b2$MCNCout['MSE', ], est3b2$MCNCout['MSE', ], est4b2$MCNCout['MSE', ])

MSE3[c(11,12,23,24), 5:8] = MSE4[c(11,12,23,24), 5:8] = NA

colnames(MSE1) = colnames(MSE2) = colnames(MSE3) = colnames(MSE4) = rep(paste('n=',Size,sep=''), 3)
m = nrow(MSE1)

### Figure C.3a ###
PATH = paste(getwd(),"/Data_and_Code_FMMCNCR/",sep="")
postscript(paste(PATH, 'results/Experiment1/figC3a.eps', sep=''), width=25, height=18)
layout(matrix(c(rep(1,7), 2,12:17, 3,18:23, 4,24:29, 5,30:35, 0,6:11), 6, 7, byrow=T), c(1.5, rep(5,6)), c(2, rep(4,4), 1))
par(mar=c(0,0,0.5,0))
plot(c(0, 1), c(0, 3), type='n', xaxt='n', yaxt='n', xlab='', ylab='', bty='n')
mtext('(a) MCN Censored Data', 1, line=-3.5, cex=1.2)
legend("bottomright", c('10% Censoring, FM-MCN', '10% Censoring, FM-MNC', '10% Censoring, FM-MCNC',
                        '20% Censoring, FM-MCN', '20% Censoring, FM-MNC', '20% Censoring, FM-MCNC'),
 col=rep(c("blue","green4","red"), 2), pch=c(2,5,1,17,18,19), lty=rep(c(2,3,1),2), bty='n', ncol=6, cex=1)

par(mar=c(0,0.5,0,0))
for(i in 1: 4){
plot(c(0, 1), c(0, 3), type='n', xaxt='n', yaxt='n', xlab='', ylab='', bty='n')
mtext('MSE', 2, line=-2, cex=1)
}

par(mar=c(0.5,0,0,0))
for(i in 1: 5){
plot(c(0, 1), c(0, 3), type='n', xaxt='n', yaxt='n', xlab='', ylab='', bty='n')
mtext('Sample Size', 1, line=-2.5, cex=0.75)
}
par(mar=c(0.5,0,0,0.5))
plot(c(0, 1), c(0, 3), type='n', xaxt='n', yaxt='n', xlab='', ylab='', bty='n')
mtext('Sample Size', 1, line=-2.5, cex=0.75)

theta.name = c(expression(w[1]), expression(mu[10]), expression(mu[11]), expression(mu[12]),
               expression(sigma[111]), expression(sigma[121]), expression(sigma[122]), expression(sigma[131]),
               expression(sigma[132]), expression(sigma[133]), expression(nu[1]), expression(rho[1]),
               expression(w[2]), expression(mu[20]), expression(mu[21]), expression(mu[22]),
               expression(sigma[211]), expression(sigma[221]), expression(sigma[222]), expression(sigma[231]),
               expression(sigma[232]), expression(sigma[233]), expression(nu[2]), expression(rho[2]))
par(mar=c(2.5,2.25,1.5,0.25))
for(r in 1: m){
 ylim1 = range(c(MSE1[r, ], MSE2[r, ]))
 plot(1:4, MSE1[r, 1:4], type='n', ylab='', xlab='', xlim=c(0.75, 4.25), ylim=ylim1, las=1, xaxt='n', cex.axis=0.75)
 axis(1, 1:4, labels=Size, cex.axis=0.75)
 mtext(theta.name[r], 3, line=0, cex=1, font.main=3)
 lines(1:4, MSE1[r, 1:4], lty=2, col='blue')
 points(1:4, MSE1[r, 1:4], pch=2, col='blue')
 lines(1:4, MSE1[r, 5:8], lty=3, col='green4')
 points(1:4, MSE1[r, 5:8], pch=5, cex=1.2, col='green4')
 lines(1:4, MSE1[r, 9:12], lty=1, col='red')
 points(1:4, MSE1[r, 9:12], pch=1, cex=1.2, col='red')

 lines(1:4, MSE2[r, 1:4], lty=2, col='blue')
 points(1:4, MSE2[r, 1:4], pch=17, col='blue')
 lines(1:4, MSE2[r, 5:8], lty=3, col='green4')
 points(1:4, MSE2[r, 5:8], pch=18, cex=1.2, col='green4')
 lines(1:4, MSE2[r, 9:12], lty=1, col='red')
 points(1:4, MSE2[r, 9:12], pch=19, cex=1.2, col='red')
}
dev.off()

### Figure C.3b ###
PATH = paste(getwd(),"/Data_and_Code_FMMCNCR/",sep="")
postscript(paste(PATH, 'results/Experiment1/figC3b.eps', sep=''), width=25, height=18)
layout(matrix(c(rep(1,7), 2,12:17, 3,18:23, 4,24:29, 5,30:35, 0,6:11), 6, 7, byrow=T), c(1.5, rep(5,6)), c(2, rep(4,4), 1))
par(mar=c(0,0,0.5,0))
plot(c(0, 1), c(0, 3), type='n', xaxt='n', yaxt='n', xlab='', ylab='', bty='n')
mtext('(b) MN Censored Data', 1, line=-3.5, cex=1.2)
legend("bottomright", c('10% Censoring, FM-MCN', '10% Censoring, FM-MNC', '10% Censoring, FM-MCNC',
                        '20% Censoring, FM-MCN', '20% Censoring, FM-MNC', '20% Censoring, FM-MCNC'),
 col=rep(c("blue","green4","red"), 2), pch=c(2,5,1,17,18,19), lty=rep(c(2,3,1),2), bty='n', ncol=6, cex=1)

par(mar=c(0,0.5,0,0))
for(i in 1: 4){
plot(c(0, 1), c(0, 3), type='n', xaxt='n', yaxt='n', xlab='', ylab='', bty='n')
mtext('MSE', 2, line=-2.5, cex=1)
}

par(mar=c(0.5,0,0,0))
for(i in 1: 4){
plot(c(0, 1), c(0, 3), type='n', xaxt='n', yaxt='n', xlab='', ylab='', bty='n')
mtext('Sample Size', 1, line=-2.5, cex=0.75)
}
par(mar=c(0.5,0,0,0))
plot(c(0, 1), c(0, 3), type='n', xaxt='n', yaxt='n', xlab='', ylab='', bty='n')
mtext('Sample Size', 1, line=-15, cex=0.75)

par(mar=c(0.5,0,0,0.5))
plot(c(0, 1), c(0, 3), type='n', xaxt='n', yaxt='n', xlab='', ylab='', bty='n')
mtext('Sample Size', 1, line=-15, cex=0.75)

theta.name = c(expression(w[1]), expression(mu[10]), expression(mu[11]), expression(mu[12]),
               expression(sigma[111]), expression(sigma[121]), expression(sigma[122]), expression(sigma[131]),
               expression(sigma[132]), expression(sigma[133]), expression(nu[1]), expression(rho[1]),
               expression(w[2]), expression(mu[20]), expression(mu[21]), expression(mu[22]),
               expression(sigma[211]), expression(sigma[221]), expression(sigma[222]), expression(sigma[231]),
               expression(sigma[232]), expression(sigma[233]), expression(nu[2]), expression(rho[2]))

par(mar=c(2.5,2.25,1.5,0.25))
for(r in 1: (m-2)){
 ylim2 = range(c(MSE3[r, ], MSE4[r, ]), na.rm=T)
 plot(1:4, MSE3[r, 1:4], type='n', ylab='', xlab='', xlim=c(0.75, 4.25), ylim=ylim2, las=1, xaxt='n', cex.axis=0.75)
 axis(1, 1:4, labels=Size, cex.axis=0.75)
 mtext(theta.name[r], 3, line=0, cex=1, font.main=3)
 lines(1:4, MSE3[r, 1:4], lty=2, col='blue')
 points(1:4, MSE3[r, 1:4], pch=2, col='blue')
 lines(1:4, MSE3[r, 5:8], lty=3, col='green4')
 points(1:4, MSE3[r, 5:8], pch=5, cex=1.2, col='green4')
 lines(1:4, MSE3[r, 9:12], lty=1, col='red')
 points(1:4, MSE3[r, 9:12], pch=1, cex=1.2, col='red')

 lines(1:4, MSE4[r, 1:4], lty=2, col='blue')
 points(1:4, MSE4[r, 1:4], pch=17, col='blue')
 lines(1:4, MSE4[r, 5:8], lty=3, col='green4')
 points(1:4, MSE4[r, 5:8], pch=18, cex=1.2, col='green4')
 lines(1:4, MSE4[r, 9:12], lty=1, col='red')
 points(1:4, MSE4[r, 9:12], pch=19, cex=1.2, col='red')
}
dev.off()
