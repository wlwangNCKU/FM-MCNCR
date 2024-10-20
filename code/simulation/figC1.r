##################################################################################
#
#   Filename: figC1.r
#   Purpose: produce Figure C.1 for MSE of recovered censored data of Experiment 1
#   Input data files: intermediate results (MSEyc.txt) sorted in the subfolders 
#                    'results/Experiment1/SIM1a/'; 'results/Experiment1/SIM1b/'; 
#                    'results/Experiment1/SIM2a/'; 'results/Experiment1/SIM2b/'; 
#                    'results/Experiment1/SIM3a/'; 'results/Experiment1/SIM3b/'; 
#                    'results/Experiment1/SIM4a/'; 'results/Experiment1/SIM4b/'; 
#   Output data files: results/Experiment1/figC1.eps
#   Required R packages: vioplot      
#
##################################################################################

library(vioplot)

PATH1a = paste(PATH, 'results/Experiment1/SIM1a/', sep='')
PATH1b = paste(PATH, 'results/Experiment1/SIM1b/', sep='')
PATH2a = paste(PATH, 'results/Experiment1/SIM2a/', sep='')
PATH2b = paste(PATH, 'results/Experiment1/SIM2b/', sep='')
PATH3a = paste(PATH, 'results/Experiment1/SIM3a/', sep='')
PATH3b = paste(PATH, 'results/Experiment1/SIM3b/', sep='')
PATH4a = paste(PATH, 'results/Experiment1/SIM4a/', sep='')
PATH4b = paste(PATH, 'results/Experiment1/SIM4b/', sep='')

# MCN #
Tab1a5 = read.table(paste(PATH1a, 'MSEyc.txt', sep=""))
Tab2a5 = read.table(paste(PATH2a, 'MSEyc.txt', sep=""))
Tab3a5 = read.table(paste(PATH3a, 'MSEyc.txt', sep=""))
Tab4a5 = read.table(paste(PATH4a, 'MSEyc.txt', sep=""))

# MVN #
Tab1b5 = read.table(paste(PATH1b, 'MSEyc.txt', sep=""))
Tab2b5 = read.table(paste(PATH2b, 'MSEyc.txt', sep=""))
Tab3b5 = read.table(paste(PATH3b, 'MSEyc.txt', sep=""))
Tab4b5 = read.table(paste(PATH4b, 'MSEyc.txt', sep=""))

name3 = c('Rep', 'rate', 'MVNC1', 'MCNC1', 'MVNC2', 'MCNC2')
colnames(Tab1a5) = colnames(Tab1b5) = colnames(Tab2a5) = colnames(Tab2b5) =
colnames(Tab3a5) = colnames(Tab3b5) = colnames(Tab4a5) = colnames(Tab4b5) = name3

postscript(paste(PATH, 'results/Experiment1/figC1.eps', sep=''), width=18, height=5)
layout(matrix(c(12,1,1,2,2, 12,3:11), 3, 5, byrow=T), c(2,4,4,5,4), c(0.5,0.5,3))
par(mar=c(0,0,0.5,0))
plot(c(0, 1), c(0, 3), type='n', xaxt='n', yaxt='n', xlab='', ylab='')
mtext('(a) MCN Censored Data', 1, line=-2.5, cex=1.75, font=2)
par(mar=c(0,3,0.5,0.5))
plot(c(0, 1), c(0, 3), type='n', xaxt='n', yaxt='n', xlab='', ylab='')
mtext('(b) MVN Censored Data', 1, line=-2.5, cex=1.75, font=2)
par(mar=c(0,0,0,0))
plot(c(0, 1), c(0, 3), type='n', xaxt='n', yaxt='n', xlab='', ylab='')
mtext('10% Censoring', 1, line=-2.75, cex=1.25, font=2)
plot(c(0, 1), c(0, 3), type='n', xaxt='n', yaxt='n', xlab='', ylab='')
mtext('20% Censoring', 1, line=-2.75, cex=1.25, font=2)
par(mar=c(0,3,0,0))
plot(c(0, 1), c(0, 3), type='n', xaxt='n', yaxt='n', xlab='', ylab='')
mtext('10% Censoring', 1, line=-2.75, cex=1.25, font=2)
par(mar=c(0,0,0,0.5))
plot(c(0, 1), c(0, 3), type='n', xaxt='n', yaxt='n', xlab='', ylab='')
mtext('20% Censoring', 1, line=-2.75, cex=1.25, font=2)

par(mar=c(3.5,0.5,0,0))
plot(c(0, 1), c(0, 3), type='n', xaxt='n', yaxt='n', xlab='', ylab='', bty='n')
mtext(expression(paste('MSE(',hat(y)^c,')',sep='')), 2, line=-3.5, cex=1.5, font=2)
nsize = c(150, 300, 600, 900)
MCN10 = rbind(cbind(Size=nsize[1], Tab1a5[which(Tab1a5$rate == 0.1), 3:4]),
              cbind(Size=nsize[2], Tab2a5[which(Tab2a5$rate == 0.1), 3:4]),
              cbind(Size=nsize[3], Tab3a5[which(Tab3a5$rate == 0.1), 3:4]),
              cbind(Size=nsize[4], Tab4a5[which(Tab4a5$rate == 0.1), 3:4]))
MCN20 = rbind(cbind(Size=nsize[1], Tab1a5[which(Tab1a5$rate == 0.2), 3:4]),
              cbind(Size=nsize[2], Tab2a5[which(Tab2a5$rate == 0.2), 3:4]),
              cbind(Size=nsize[3], Tab3a5[which(Tab3a5$rate == 0.2), 3:4]),
              cbind(Size=nsize[4], Tab4a5[which(Tab4a5$rate == 0.2), 3:4]))

MVN10 = rbind(cbind(Size=nsize[1], Tab1b5[which(Tab1b5$rate == 0.1), 3:4]),
              cbind(Size=nsize[2], Tab2b5[which(Tab2b5$rate == 0.1), 3:4]),
              cbind(Size=nsize[3], Tab3b5[which(Tab3b5$rate == 0.1), 3:4]),
              cbind(Size=nsize[4], Tab4b5[which(Tab4b5$rate == 0.1), 3:4]))
MVN20 = rbind(cbind(Size=nsize[1], Tab1b5[which(Tab1b5$rate == 0.2), 3:4]),
              cbind(Size=nsize[2], Tab2b5[which(Tab2b5$rate == 0.2), 3:4]),
              cbind(Size=nsize[3], Tab3b5[which(Tab3b5$rate == 0.2), 3:4]),
              cbind(Size=nsize[4], Tab4b5[which(Tab4b5$rate == 0.2), 3:4]))

par(mar=c(3.5,0,0,0))
ylim1 = range(rbind(MCN10, MCN20)[,-1])
ylim2 = range(rbind(MVN10, MVN20)[,-1])

vioplot(MVNC1~Size, data=MCN10, ylim=ylim1, col = "yellow", rectCol="green", lineCol="green4", border="green4", plotCentre = "line", side = "left", xlab='', las=1, cex.axis=1.2)
vioplot(MCNC1~Size, data=MCN10, ylim=ylim1, col = "pink", rectCol="red", lineCol="red4", border="red4", plotCentre = "line", side = "right", las=1, add=T, cex.axis=1.2)
mtext('Sample Size', 1, line=2.5, cex=0.75)

vioplot(MVNC1~Size, data=MCN20, ylim=ylim1, col = "yellow", rectCol="green", lineCol="green4", border="green4", plotCentre = "line", side = "left", xlab='', las=1, yaxt='n', cex.axis=1.2)
vioplot(MCNC1~Size, data=MCN20, ylim=ylim1, col = "pink", rectCol="red", lineCol="red4", border="red4", plotCentre = "line", side = "right", las=1, add=T, yaxt='n', cex.axis=1.2)
mtext('Sample Size', 1, line=2.5, cex=0.75)
axis(2, seq(1, 5, 1), labels=F)

par(mar=c(3.5,3,0,0))
vioplot(MVNC1~Size, data=MVN10, ylim=ylim2, col = "yellow", rectCol="green", lineCol="green4", border="green4", plotCentre = "line", side = "left", xlab='', las=1, cex.axis=1.2)
vioplot(MCNC1~Size, data=MVN10, ylim=ylim2, col = "pink", rectCol="red", lineCol="red4", border="red4", plotCentre = "line", side = "right", las=1, add=T, cex.axis=1.5)
mtext('Sample Size', 1, line=2.5, cex=0.75)

par(mar=c(3.5,0,0,0.5))
vioplot(MVNC1~Size, data=MVN20, ylim=ylim2, col = "yellow", rectCol="green", lineCol="green4", border="green4", plotCentre = "line", side = "left", xlab='', las=1, yaxt='n', cex.axis=1.2)
vioplot(MCNC1~Size, data=MVN20, ylim=ylim2, col = "pink", rectCol="red", lineCol="red4", border="red4", plotCentre = "line", side = "right", las=1, add=T, yaxt='n', cex.axis=1.2)
axis(2, seq(0.1, 0.5, 0.1), labels=F)
mtext('Sample Size', 1, line=2.5, cex=0.75)

par(mar=c(0,0.5,0.5,0))
plot(c(0, 1), c(0, 3), type='n', xaxt='n', yaxt='n', xlab='', ylab='', bty='n')
legend("center", c('FM-MNC','FM-MCNC'), fill=c('yellow','pink'), border=c('green4','red4'), cex=1.2, bty='n')
dev.off()
