#install.packages("sampleSelection")

################################################################################
#
#   Filename: fig2.r
#   Purpose: produce Figure 2 for the Mroz87 data
#   Input data files: StdWageFit.RData
#   Output data files: results/fig2.eps
#
################################################################################

load(paste(PATH, 'data/StdWageFit.RData', sep=''))

PATH = paste(getwd(),"/Data_and_Code_FMMCNCR/",sep="")
postscript(paste(PATH, 'results/fig2.eps', sep=''), width=8, height=20) 
layout(matrix(c(1,5,6,7,11,2,8,9,12,14,3,10,13,15,16,4), 4, 4, byrow=T), c(5.5, 5, 5, 5.5), c(5.5, 5, 5, 5.5))
yname = c("Wife's Work Hours", "Husband's Work Hours", "Wife's Wage", "Husband's Wage")
Yfit = matrix(estCC2$yfit$yfit1, nrow=n, ncol=p, byrow=T)
Range1 = as.list(p)
for(k in 1: p) Range1[[k]] = range(Yc1[,k])
par(mar=c(0.25, 2, 2, 0.25))
a = hist(Yc1[,1], nclass=20, plot=F)
plot(Range1[[1]], c(0, max(a$density)), type='n', xaxt='n', yaxt='n', xlab='', ylab='')
axis(3, seq(-1, 5, 1), cex.axis=1.15)
hist(Yc1[,1], nclass=20, prob=T, las=1, cex.axis=0.5, main='', xaxt='n', yaxt='n', col=c('gray',rep('white',100)), add=T)
mtext(yname[1], 3, line=-2, cex=0.9, font=2)

par(mar=rep(0.25, 4))
a = hist(Yc1[,2], nclass=20, plot=F)
plot(Range1[[2]], c(0, max(a$density)), type='n', xaxt='n', yaxt='n', xlab='', ylab='')
hist(Yc1[,2], nclass=25, prob=T, las=1, cex.axis=0.5, main='', xaxt='n', yaxt='n', col=rep('white',100), add=T)
mtext(yname[2], 3, line=-2, cex=0.9, font=2)

a = hist(Yc1[,3], nclass=20, plot=F)
plot(Range1[[3]], c(0, max(a$density)), type='n', xaxt='n', yaxt='n', xlab='', ylab='')
hist(Yc1[,3], nclass=25, prob=T, las=1, cex.axis=0.5, main='', xaxt='n', yaxt='n', col=c('gray',rep('white',100)), add=T)
mtext(yname[3], 3, line=-2, cex=0.9, font=2)

par(mar=c(2, 0.25, 0.25, 2))
a = hist(Yc1[,4], nclass=20, plot=F)
plot(Range1[[4]], c(0, max(a$density)), type='n', xaxt='n', yaxt='n', xlab='', ylab='')
axis(1, seq(0, 8, 2), cex.axis=1.15)
hist(Yc1[,4], nclass=25, prob=T, las=1, cex.axis=0.5, main='', xaxt='n', yaxt='n', col=rep('white',100), add=T)
mtext(yname[4], 3, line=-2, cex=0.9, font=2)

out1 = which(estCC2$xi$xi1 > 0.5)
cenID = which(rowSums(cen) != 0)
outID = out1[which(out1 %in% cenID == F)] 

par(mar=c(0.25, 0.25, 2, 0.25))
plot(Yc1[,2], Yc1[,1], xlim=Range1[[2]], ylim=Range1[[1]], cex=0.9, pch=20, las=1, cex.axis=1, lwd=0.25, xaxt='n', yaxt='n')
points(Yc1[outID, 2], Yc1[outID, 1], cex=1.5, pch=1, col='green4')
axis(3, seq(-2,8,2), cex.axis=1.15)
abline(h = c(u1.std[1], u2.std[1]), lty=2, col=gray(.5))

plot(Yc1[,3], Yc1[,1], xlim=Range1[[3]], ylim=Range1[[1]], cex=0.9, pch=20, las=1, cex.axis=1, lwd=0.25, xaxt='n', yaxt='n')
points(Yc1[outID, 3], Yc1[outID, 1], cex=1.5, pch=1, col='green4')
axis(3, seq(-2,8,2), cex.axis=1.15)
abline(v = c(u1.std[3], u2.std[3]), lty=2, col=gray(.5))
abline(h = c(u1.std[1], u2.std[1]), lty=2, col=gray(.5))

par(mar=c(0.25, 0.25, 2, 2))
plot(Yc1[,4], Yc1[,1], xlim=Range1[[4]], ylim=Range1[[1]], cex=0.9, pch=20, las=1, cex.axis=1, lwd=0.25, xaxt='n', yaxt='n')
points(Yc1[outID, 4], Yc1[outID, 1], cex=1.5, pch=1, col='green4')
axis(3, seq(-2,8,2), cex.axis=1.15)
axis(4, seq(-2,5,1), cex.axis=1.15, las=1)
abline(h = c(u1.std[1], u2.std[1]), lty=2, col=gray(.5))

par(mar=c(0.25, 0.25, 0.25, 0.25))
plot(Yc1[,3], Yc1[,2], xlim=Range1[[3]], ylim=Range1[[2]], cex=0.9, pch=20, las=1, cex.axis=1, lwd=0.25, xaxt='n', yaxt='n')
points(Yc1[outID, 3], Yc1[outID, 2], cex=1.5, pch=1, col='green4')
abline(v = c(u1.std[3], u2.std[3]), lty=2, col=gray(.5))

par(mar=c(0.25, 0.25, 0.25, 2))
plot(Yc1[,4], Yc1[,2], xlim=Range1[[4]], ylim=Range1[[2]], cex=0.9, pch=20, las=1, cex.axis=1, lwd=0.25, xaxt='n', yaxt='n')
points(Yc1[outID, 4], Yc1[outID, 2], cex=1.5, pch=1, col='green4')
axis(4, seq(-2,8,2), cex.axis=1.15, las=1)

plot(Yc1[,4], Yc1[,3], xlim=Range1[[4]], ylim=Range1[[3]], cex=0.9, pch=20, las=1, cex.axis=1, lwd=0.25, xaxt='n', yaxt='n')
points(Yc1[outID, 4], Yc1[outID, 3], cex=1.5, pch=1, col='green4')
axis(4, seq(-2,4,2), cex.axis=1.15, las=1)
abline(h = c(u1.std[3], u2.std[3]), lty=2, col=gray(.5))

Range2 = as.list(p)
for(k in 1: p) Range2[[k]] = range(Yfit[,k])
par(mar=c(0.25, 2, 0.25, 0.25))
plot(Yfit[,1], Yfit[,2], xlim=Range2[[1]], ylim=Range2[[2]], col=c('red','blue')[unclass(estCC2$pre.cls$post.cls)], pch=c(3,2)[unclass(estCC2$pre.cls$post.cls)], cex=0.9, lwd=0.4, xaxt='n', yaxt='n')
axis(2, seq(-0.4,0.4,0.2), cex.axis=0.9, las=1)

plot(Yfit[,1], Yfit[,3], xlim=Range2[[1]], ylim=Range2[[3]], col=c('red','blue')[unclass(estCC2$pre.cls$post.cls)], pch=c(3,2)[unclass(estCC2$pre.cls$post.cls)], cex=0.9, lwd=0.4, xaxt='n', yaxt='n')
axis(2, seq(-1,1,0.5), cex.axis=0.9, las=1)

par(mar=c(2, 2, 0.25, 0.25))
plot(Yfit[,1], Yfit[,4], xlim=Range2[[1]], ylim=Range2[[4]], col=c('red','blue')[unclass(estCC2$pre.cls$post.cls)], pch=c(3,2)[unclass(estCC2$pre.cls$post.cls)], cex=0.9, lwd=0.4, xaxt='n', yaxt='n')
axis(1, seq(-1.5,1,0.5), cex.axis=1.15, las=1)
axis(2, seq(-1,0.5,0.5), cex.axis=0.9, las=1)

par(mar=c(0.25, 0.25, 0.25, 0.25))
plot(Yfit[,2], Yfit[,3], xlim=Range2[[2]], ylim=Range2[[3]], col=c('red','blue')[unclass(estCC2$pre.cls$post.cls)], pch=c(3,2)[unclass(estCC2$pre.cls$post.cls)], cex=0.9, lwd=0.4, xaxt='n', yaxt='n')

par(mar=c(2, 0.25, 0.25, 0.25))
plot(Yfit[,2], Yfit[,4], xlim=Range2[[2]], ylim=Range2[[4]], col=c('red','blue')[unclass(estCC2$pre.cls$post.cls)], pch=c(3,2)[unclass(estCC2$pre.cls$post.cls)], cex=0.9, lwd=0.4, xaxt='n', yaxt='n')
axis(1, seq(-0.4,0.4,0.2), cex.axis=1.15)

plot(Yfit[,3], Yfit[,4], xlim=Range2[[3]], ylim=Range2[[4]], col=c('red','blue')[unclass(estCC2$pre.cls$post.cls)], pch=c(3,2)[unclass(estCC2$pre.cls$post.cls)], cex=0.9, lwd=0.4, xaxt='n', yaxt='n')
axis(1, seq(-1,1,0.5), cex.axis=1.15)
dev.off()
