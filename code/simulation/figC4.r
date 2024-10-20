##################################################################################
#
#   Filename: figC4.r
#   Purpose: produce Figure C.4 for MSE of theta.hat and y.fit in Experiment 2
#   Input data files: intermediate results (MSEest.txt and MSEyfit.txt) 
#                     sorted in the subfolders 
#                     'results/Experiment2/SIM1/'; 'results/Experiment2/SIM2/'; 
#                     'results/Experiment2/SIM3/'; 'results/Experiment2/SIM4/'; 
#   Output data files: results/Experiment1/figC4.eps
#   Required R packages: vioplot      
#
##################################################################################

library(vioplot)

PATH1 = paste(PATH, 'results/Experiment2/SIM1/', sep='')
PATH2 = paste(PATH, 'results/Experiment2/SIM2/', sep='')
PATH3 = paste(PATH, 'results/Experiment2/SIM3/', sep='')
PATH4 = paste(PATH, 'results/Experiment2/SIM4/', sep='')

Tab14 = read.table(paste(PATH1, 'MSEest.txt',sep=""))
Tab15 = read.table(paste(PATH1, 'MSEyfit.txt',sep=""))
Tab24 = read.table(paste(PATH2, 'MSEest.txt',sep=""))
Tab25 = read.table(paste(PATH2, 'MSEyfit.txt',sep=""))
Tab34 = read.table(paste(PATH3, 'MSEest.txt',sep=""))
Tab35 = read.table(paste(PATH3, 'MSEyfit.txt',sep=""))
Tab44 = read.table(paste(PATH4, 'MSEest.txt',sep=""))
Tab45 = read.table(paste(PATH4, 'MSEyfit.txt',sep=""))

name4 = c('Rep','rate','MCNR','MCNCR')
name5 = c('Rep','rate','yfitC11','yfitCC11','yfitC12','yfitCC12')
names(Tab14) = names(Tab24) = names(Tab34) = names(Tab44) = name4
names(Tab15) = names(Tab25) = names(Tab35) = names(Tab45) = name5

Size = c(150, 300, 600, 900)
# MSE: Parameter estimates #
MSEest1 = rbind(
cbind(Size=Size[1], Tab14[which(Tab14$rate == 0.1), 3:4]),
cbind(Size=Size[2], Tab24[which(Tab24$rate == 0.1), 3:4]),
cbind(Size=Size[3], Tab34[which(Tab34$rate == 0.1), 3:4]),
cbind(Size=Size[4], Tab44[which(Tab44$rate == 0.1), 3:4]))

MSEest2 = rbind(
cbind(Size=Size[1], Tab14[which(Tab14$rate == 0.2), 3:4]),
cbind(Size=Size[2], Tab24[which(Tab24$rate == 0.2), 3:4]),
cbind(Size=Size[3], Tab34[which(Tab34$rate == 0.2), 3:4]),
cbind(Size=Size[4], Tab44[which(Tab44$rate == 0.2), 3:4]))

# MSE: Fitted values #
MSEyfit1 = rbind(
cbind(Size=Size[1], Tab15[which(Tab15$rate == 0.1), 3:4]),
cbind(Size=Size[2], Tab25[which(Tab25$rate == 0.1), 3:4]),
cbind(Size=Size[3], Tab35[which(Tab35$rate == 0.1), 3:4]),
cbind(Size=Size[4], Tab45[which(Tab45$rate == 0.1), 3:4]))

MSEyfit2 = rbind(
cbind(Size=Size[1], Tab15[which(Tab15$rate == 0.2), 3:4]),
cbind(Size=Size[2], Tab25[which(Tab25$rate == 0.2), 3:4]),
cbind(Size=Size[3], Tab35[which(Tab35$rate == 0.2), 3:4]),
cbind(Size=Size[4], Tab45[which(Tab45$rate == 0.2), 3:4]))

postscript(paste(PATH, 'results/Experiment2/figC4.eps', sep=''), width=20, height=15)
layout(matrix(c(0,2,3,1,4,5,0,7,8,6,9,10), 4, 3, byrow=T), c(1,5,5), c(1,5,1,5))
# MSE: MLE #
par(mar=c(3.5,0.5,0,0))
plot(c(0, 1), c(0, 3), type='n', xaxt='n', yaxt='n', xlab='', ylab='', bty='n')
mtext(expression(paste('MSE(',hat(theta),')',sep='')), 2, line=-3, cex=1.5, font=2)
par(mar=c(0,0,0.5,2))
plot(c(0, 1), c(0, 3), type='n', xaxt='n', yaxt='n', xlab='', ylab='')
mtext('10% Censoring', 1, line=-2, cex=1.5, font=2)
par(mar=c(0,1,0.5,0.5))
plot(c(0, 1), c(0, 3), type='n', xaxt='n', yaxt='n', xlab='', ylab='')
mtext('20% Censoring', 1, line=-2, cex=1.5, font=2)
par(mar=c(3.5,0,0,2))
ylim1 = c(0, 2)
vioplot(MCNR~Size, data=MSEest1, ylim=ylim1, col = "cyan", rectCol="blue", lineCol="royalblue", border="royalblue", plotCentre = "line", side = "left", xlab='', las=1, cex.axis=1.5)
vioplot(MCNCR~Size, data=MSEest1, ylim=ylim1, col = "pink", rectCol="red", lineCol="red4", border="red4", plotCentre = "line", side = "right", las=1, add=T, cex.axis=1.5)
mtext('Sample Size', 1, line=2.5, cex=0.9)
legend("topright", c('FM-MCNR','FM-MCNCR'), fill=c('cyan','pink'), border=c('royalblue','red4'), cex=1.2, bty='n')
par(mar=c(3.5,1,0,0.5))
ylim2 = c(0, 8)
vioplot(MCNR~Size, data=MSEest2, ylim=ylim2, col = "cyan", rectCol="blue", lineCol="royalblue", border="royalblue", plotCentre = "line", side = "left", xlab='', las=1, cex.axis=1.5)
vioplot(MCNCR~Size, data=MSEest2, ylim=ylim2, col = "pink", rectCol="red", lineCol="red4", border="red4", plotCentre = "line", side = "right", las=1, add=T)
mtext('Sample Size', 1, line=2.5, cex=0.9)
legend("topright", c('FM-MCNR','FM-MCNCR'), fill=c('cyan','pink'), border=c('royalblue','red4'), cex=1.2, bty='n')

# MSE: yfit #
par(mar=c(3.5,0.5,0,0))
plot(c(0, 1), c(0, 3), type='n', xaxt='n', yaxt='n', xlab='', ylab='', bty='n')
mtext(expression(paste('MSE(',hat(y),')',sep='')), 2, line=-3, cex=1.5, font=2)
par(mar=c(0,0,0.5,2))
plot(c(0, 1), c(0, 3), type='n', xaxt='n', yaxt='n', xlab='', ylab='')
mtext('10% Censoring', 1, line=-2, cex=1.5, font=2)
par(mar=c(0,1,0.5,0.5))
plot(c(0, 1), c(0, 3), type='n', xaxt='n', yaxt='n', xlab='', ylab='')
mtext('20% Censoring', 1, line=-2, cex=1.5, font=2)
par(mar=c(3.5,0,0,2))
ylim1 = range(MSEyfit1[,2:3])
vioplot(yfitC11~Size, data=MSEyfit1, ylim=ylim1, col = "cyan", rectCol="blue", lineCol="royalblue", border="royalblue", plotCentre = "line", side = "left", xlab='', las=1, cex.axis=1.5)
vioplot(yfitCC11~Size, data=MSEyfit1, ylim=ylim1, col = "pink", rectCol="red", lineCol="red4", border="red4", plotCentre = "line", side = "right", las=1, add=T)
mtext('Sample Size', 1, line=2.5, cex=0.9)
legend("topright", c('FM-MCNR','FM-MCNCR'), fill=c('cyan','pink'), border=c('royalblue','red4'), cex=1.2, bty='n')
par(mar=c(3.5,1,0,0.5))
ylim2 = range(MSEyfit2[,2:3])
vioplot(yfitC11~Size, data=MSEyfit2, ylim=ylim2, col = "cyan", rectCol="blue", lineCol="royalblue", border="royalblue", plotCentre = "line", side = "left", xlab='', las=1, cex.axis=1.5)
vioplot(yfitCC11~Size, data=MSEyfit2, ylim=ylim2, col = "pink", rectCol="red", lineCol="red4", border="red4", plotCentre = "line", side = "right", las=1, add=T)
mtext('Sample Size', 1, line=2.5, cex=0.9)
legend("topright", c('FM-MCNR','FM-MCNCR'), fill=c('cyan','pink'), border=c('royalblue','red4'), cex=1.2, bty='n')
dev.off()
