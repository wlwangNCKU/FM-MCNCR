################################################################################
#
#   Filename: fig1.r
#   Purpose: produce Figure 1 for the Wholesale data
#   Input data files: data/SaleFit.RData
#   Output data files: results/fig1.eps
#
################################################################################

load(paste(PATH, 'data/SaleFit.RData', sep=''))

# Visualization #
est = estCC2$para
m = 200
lim1 = c(-2, 12)
xx = expand.grid(x1<-seq(lim1[1], lim1[2], length=m), x2<-seq(lim1[1], lim1[2], length=m))
d12 = d13 = d14 = d23 = d24 = d34 = matrix(NA, m, m)
for(j in 1: m)for(k in 1: m){
 d12[j, k] = dfmmcn(c(x1[j], x2[k]), g=2, w=est$w, mu=est$mu[c(1,2),], Sigma=est$Sigma[c(1,2), c(1,2), ], nu=est$nu, rho=est$rho, distr='MCN')
 d13[j, k] = dfmmcn(c(x1[j], x2[k]), g=2, w=est$w, mu=est$mu[c(1,3),], Sigma=est$Sigma[c(1,3), c(1,3), ], nu=est$nu, rho=est$rho, distr='MCN')
 d14[j, k] = dfmmcn(c(x1[j], x2[k]), g=2, w=est$w, mu=est$mu[c(1,4),], Sigma=est$Sigma[c(1,4), c(1,4), ], nu=est$nu, rho=est$rho, distr='MCN')
 d23[j, k] = dfmmcn(c(x1[j], x2[k]), g=2, w=est$w, mu=est$mu[c(2,3),], Sigma=est$Sigma[c(2,3), c(2,3), ], nu=est$nu, rho=est$rho, distr='MCN')
 d24[j, k] = dfmmcn(c(x1[j], x2[k]), g=2, w=est$w, mu=est$mu[c(2,4),], Sigma=est$Sigma[c(2,4), c(2,4), ], nu=est$nu, rho=est$rho, distr='MCN')
 d34[j, k] = dfmmcn(c(x1[j], x2[k]), g=2, w=est$w, mu=est$mu[c(3,4),], Sigma=est$Sigma[c(3,4), c(3,4), ], nu=est$nu, rho=est$rho, distr='MCN')
}

PATH = paste(getwd(),"/Data_and_Code_FMMCNCR/",sep="")
postscript(paste(PATH, 'results/fig1.eps', sep=''), width=8, height=20)
layout(matrix(c(1,5,6,7,11,2,8,9,12,14,3,10,13,15,16,4), 4, 4, byrow=T), c(5.5, 5, 5, 5.5), c(5.5, 5, 5, 5.5))

# Diagonal
PAR = as.list(p)
PAR[[1]] = c(0.25, 2, 2, 0.25)
PAR[[2]] = PAR[[3]] = rep(0.25, 4)
PAR[[4]] = c(2, 0.25, 0.25, 2)
for(k in 1: p){
par(mar = PAR[[k]])
yk = list(yk1 = Yc2[cls==1, k], yk2 = Yc2[cls==2, k])
h1 = lapply(yk, hist, breaks = seq(min(Yc2[,k]), max(Yc2[,k]), length=30), plot=F)
t1 = rbind(h1[[1]]$counts, h1[[2]]$counts)
rownames(t1)=names(h1)
colnames(t1)=h1[[1]]$mids
a1 = 1: length(h1$yk1$breaks) - 1
plot(c(0,30),c(0,max(t1)), type='n',xlab='',ylab='',xaxt='n',yaxt='n')
barplot(t1[1,a1], ylim=c(0,max(t1)), col=0, border="red", lwd=2, space=0, xaxt='n', yaxt='n', add=T)
barplot(t1[2,a1], ylim=c(0,max(t1)), col='cyan', border='blue', las=1, space=0, add=T, xaxt='n', yaxt='n')
mtext(colnames(dataset)[k], 1, line=-10, font=2, cex=1.15)
legend('center', c('Retail; Cluster 1', 'Horeca; Cluster 2', 'Outlier'), pch=c(19, 4, 1), col=c('red','blue', 'green4'), bty='n', cex=1.15)
legend('right', c('', '', ''), fill=c('white', 'cyan', 'white'), border=c('red', 'blue', 'white'), bty='n', cex=1.15)
}

# Upper triangular
Range = as.list(p)
for(k in 1: p) Range[[k]] = range(c(Yc2[,k], estCC2$yhat$yhat1[,k]))

out1 = which(estCC2$xi$xi1 > 0.5)
cenID = which(rowSums(cen) != 0)
outID = out1[which(out1 %in% cenID == F)]
stdubd = as.numeric(p)
for(k in 1:p) stdubd[k] = unique(Yc2[which(cen[,k]==1), k])

par(mar=c(0.25, 0.25, 2, 0.25))
for(j in 2:3){
 plot(Yc2[,j], Yc2[,1], xlim=Range[[j]], ylim=Range[[1]], cex=0.6, pch=c(19, 4)[unclass(cls)], col=c('red','blue')[unclass(cls)], las=1, cex.axis=1.15, lwd=0.5, xaxt='n', yaxt='n')
 points(Yc2[outID, j], Yc2[outID, 1], cex=1.5, pch=1, col='green4')
 axis(3, seq(0,12,2), cex.axis=1.15)
 abline(v = stdubd[j], lty=2, col=gray(.5))
 abline(h = stdubd[1], lty=2, col=gray(.5))
}
par(mar=c(0.25, 0.25, 2, 2))
plot(Yc2[,4], Yc2[,1], xlim=Range[[4]], ylim=Range[[1]], cex=0.6, pch=c(19, 4)[unclass(cls)], col=c('red','blue')[unclass(cls)], las=1, cex.axis=1.15, lwd=0.5, xaxt='n', yaxt='n')
points(Yc2[outID, 4], Yc2[outID, 1], cex=1.5, pch=1, col='green4')
axis(3, seq(0,8,2), cex.axis=1.15)
axis(4, seq(0,8,2), cex.axis=1.15, las=1)
abline(v = stdubd[4], lty=2, col=gray(.5))
abline(h = stdubd[1], lty=2, col=gray(.5))

par(mar=c(0.25, 0.25, 0.25, 0.25))
plot(Yc2[,3], Yc2[,2], xlim=Range[[3]], ylim=Range[[2]], cex=0.6, pch=c(19, 4)[unclass(cls)], col=c('red','blue')[unclass(cls)], las=1, cex.axis=1.15, lwd=0.5, xaxt='n', yaxt='n')
points(Yc2[outID, 3], Yc2[outID, 2], cex=1.5, pch=1, col='green4')
abline(v = stdubd[3], lty=2, col=gray(.5))
abline(h = stdubd[2], lty=2, col=gray(.5))

par(mar=c(0.25, 0.25, 0.25, 2))
plot(Yc2[,4], Yc2[,2], xlim=Range[[4]], ylim=Range[[2]], cex=0.6, pch=c(19, 4)[unclass(cls)], col=c('red','blue')[unclass(cls)], las=1, cex.axis=1.15, lwd=0.5, xaxt='n', yaxt='n')
points(Yc2[outID, 4], Yc2[outID, 2], cex=1.5, pch=1, col='green4')
axis(4, seq(0,12,2), cex.axis=1.15, las=1)
abline(v = stdubd[4], lty=2, col=gray(.5))
abline(h = stdubd[2], lty=2, col=gray(.5))

plot(Yc2[,4], Yc2[,3], xlim=Range[[4]], ylim=Range[[3]], cex=0.6, pch=c(19, 4)[unclass(cls)], col=c('red','blue')[unclass(cls)], las=1, cex.axis=1.15, lwd=0.5, xaxt='n', yaxt='n')
points(Yc2[outID, 4], Yc2[outID, 3], cex=1.5, pch=1, col='green4')
axis(4, seq(0,12,2), cex.axis=1.15, las=1)
abline(v = stdubd[4], lty=2, col=gray(.5))
abline(h = stdubd[3], lty=2, col=gray(.5))

# Lower triangular
pre.cls = unclass(estCC2$pre.cls$post.cls)
Yhat = estCC2$yhat$yhat1
setlevel = c(0.2,0.1,0.05,0.005,0.002,0.001,0.0003,5e-5,5e-6)

par(mar=c(0.25, 2, 0.25, 0.25))
plot(Yhat[,1], Yhat[,2], xlim=Range[[1]], ylim=Range[[2]], type='n', xaxt='n', yaxt='n')
axis(2, seq(0,12,2), cex.axis=1.15, las=1)
abline(v = stdubd[1], lty=2, col=gray(.5))
abline(h = stdubd[2], lty=2, col=gray(.5))
contour(x1, x2, d12, add=T, col=gray(.6), lwd=0.01, drawlabels=F, levels=setlevel)
points(Yhat[,1], Yhat[,2], col=c('red','blue')[pre.cls], pch=c(19,4)[pre.cls], cex=0.5, lwd=0.5)

plot(Yhat[,1], Yhat[,3], xlim=Range[[1]], ylim=Range[[3]], type='n', xaxt='n', yaxt='n')
axis(2, seq(0,12,2), cex.axis=1.15, las=1)
abline(v = stdubd[1], lty=2, col=gray(.5))
abline(h = stdubd[3], lty=2, col=gray(.5))
contour(x1, x2, d13, add=T, col=gray(.6), lwd=0.01, drawlabels=F, levels=setlevel)
points(Yhat[,1], Yhat[,3], col=c('red','blue')[pre.cls], pch=c(19,4)[pre.cls], cex=0.5, lwd=0.5)

par(mar=c(2, 2, 0.25, 0.25))
plot(Yhat[,1], Yhat[,4], xlim=Range[[1]], ylim=Range[[4]], type='n', xaxt='n', yaxt='n')
abline(v =stdubd[1], lty=2, col=gray(.5))
abline(h = stdubd[4], lty=2, col=gray(.5))
axis(1, seq(0,8,2), cex.axis=1.15, las=1)
axis(2, seq(0,8,2), cex.axis=1.15, las=1)
contour(x1, x2, d14, add=T, col=gray(.6), lwd=0.01, drawlabels=F, levels=setlevel)
points(Yhat[,1], Yhat[,4], col=c('red','blue')[pre.cls], pch=c(19,4)[pre.cls], cex=0.5, lwd=0.5)

par(mar=c(0.25, 0.25, 0.25, 0.25))
plot(Yhat[,2], Yhat[,3], xlim=Range[[2]], ylim=Range[[3]], type='n', xaxt='n', yaxt='n')
abline(v = stdubd[2], lty=2, col=gray(.5))
abline(h = stdubd[3], lty=2, col=gray(.5))
contour(x1, x2, d23, add=T, col=gray(.6), lwd=0.01, drawlabels=F, levels=setlevel)
points(Yhat[,2], Yhat[,3], col=c('red','blue')[pre.cls], pch=c(19,4)[pre.cls], cex=0.5, lwd=0.5)

par(mar=c(2, 0.25, 0.25, 0.25))
plot(Yhat[,2], Yhat[,4], xlim=Range[[2]], ylim=Range[[4]], type='n', xaxt='n', yaxt='n')
axis(1, seq(0,8,2), cex.axis=1.15)
abline(v = stdubd[2], lty=2, col=gray(.5))
abline(h = stdubd[4], lty=2, col=gray(.5))
contour(x1, x2, d24, add=T, col=gray(.6), lwd=0.01, drawlabels=F, levels=setlevel)
points(Yhat[,2], Yhat[,4], col=c('red','blue')[pre.cls], pch=c(19,4)[pre.cls], cex=0.5, lwd=0.5)

plot(Yhat[,3], Yhat[,4], xlim=Range[[3]], ylim=Range[[4]], type='n', xaxt='n', yaxt='n')
axis(1, seq(0,12,2), cex.axis=1.15)
abline(v = stdubd[3], lty=2, col=gray(.5))
abline(h = stdubd[4], lty=2, col=gray(.5))
contour(x1, x2, d34, add=T, col=gray(.6), lwd=0.01, drawlabels=F, levels=setlevel)
points(Yhat[,3], Yhat[,4], col=c('red','blue')[pre.cls], pch=c(19,4)[pre.cls], cex=0.5, lwd=0.5)
dev.off()
