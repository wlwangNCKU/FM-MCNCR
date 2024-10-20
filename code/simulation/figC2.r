##################################################################################
#
#   Filename: figC2.r
#   Purpose: produce Figure C.2 for CCR and ARI of Experiment 1 
#   Input data files: intermediate results (cls.txt) sorted in the subfolders 
#                    'results/Experiment1/SIM1a/'; 'results/Experiment1/SIM1b/'; 
#                    'results/Experiment1/SIM2a/'; 'results/Experiment1/SIM2b/'; 
#                    'results/Experiment1/SIM3a/'; 'results/Experiment1/SIM3b/'; 
#                    'results/Experiment1/SIM4a/'; 'results/Experiment1/SIM4b/'; 
#   Output data files: results/Experiment1/figC2.eps
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
Tab1a2 = read.table(paste(PATH1a, 'cls.txt',sep=""))
Tab2a2 = read.table(paste(PATH2a, 'cls.txt',sep=""))
Tab3a2 = read.table(paste(PATH3a, 'cls.txt',sep=""))
Tab4a2 = read.table(paste(PATH4a, 'cls.txt',sep=""))

# MVN #
Tab1b2 = read.table(paste(PATH1b, 'cls.txt',sep=""))
Tab2b2 = read.table(paste(PATH2b, 'cls.txt',sep=""))
Tab3b2 = read.table(paste(PATH3b, 'cls.txt',sep=""))
Tab4b2 = read.table(paste(PATH4b, 'cls.txt',sep=""))

name1 = c('Rep', 'rate', 'criterion', 'MCN', 'MVNC', 'MCNC')

colnames(Tab1a2) = colnames(Tab1b2) = colnames(Tab2a2) = colnames(Tab2b2) =
colnames(Tab3a2) = colnames(Tab3b2) = colnames(Tab4a2) = colnames(Tab4b2) = name1

postscript(paste(PATH, 'results/Experiment1/figC2.eps', sep=''), width=12, height=18)
layout(matrix(c(1,rep(3,8),rep(4,8), 1,rep(5,4),rep(6,4),rep(7,4),rep(8,4), 1,9:24, 2,25:40, 41,rep(43,8),rep(44,8), 41,rep(45,4),rep(46,4),rep(47,4),rep(48,4), 41,49:64, 42,65:80), 8, 17, byrow=T),
       c(3,rep(2,7),2.05, 4, rep(2,6),2.05), rep(c(1,1,1,5),2))
# CCR
par(mar=c(0,0.5,0.5,0))
plot(c(0, 1), c(0, 3), type='n', xaxt='n', yaxt='n', xlab='', ylab='')
legend('center', c('1:FM-MCN', '2:FM-MNC', '3:FM-MCNC'), border=c('blue','green4','red'), fill=c('cyan','yellow','pink'), cex=0.9, bty='n')
par(mar=c(3.5,0.5,0,0))
plot(c(0, 1), c(0, 3), type='n', xaxt='n', yaxt='n', xlab='', ylab='')
mtext('CCR', 2, line=-3, cex=1.5, font=2)
par(mar=c(0,0,0.5,0))
plot(c(0, 1), c(0, 3), type='n', xaxt='n', yaxt='n', xlab='', ylab='')
mtext('(a) MCN Censored Data', 1, line=-2, cex=1.5, font=2)
par(mar=c(0,3.5,0.5,0.5))
plot(c(0, 1), c(0, 3), type='n', xaxt='n', yaxt='n', xlab='', ylab='')
mtext('(c) MN Censored Data', 1, line=-2, cex=1.5, font=2)
par(mar=c(0,0,0,0))
plot(c(0, 1), c(0, 3), type='n', xaxt='n', yaxt='n', xlab='', ylab='')
mtext('10% Censoring', 1, line=-2.25, cex=1.2, font=2)
plot(c(0, 1), c(0, 3), type='n', xaxt='n', yaxt='n', xlab='', ylab='')
mtext('20% Censoring', 1, line=-2.25, cex=1.2, font=2)
par(mar=c(0,3.5,0,0))
plot(c(0, 1), c(0, 3), type='n', xaxt='n', yaxt='n', xlab='', ylab='')
mtext('10% Censoring', 1, line=-2.25, cex=1.2, font=2)
par(mar=c(0,0,0,0.5))
plot(c(0, 1), c(0, 3), type='n', xaxt='n', yaxt='n', xlab='', ylab='')
mtext('20% Censoring', 1, line=-2.25, cex=1.2, font=2)
par(mar=c(0,0,0,0))
n = c(150, 300, 600, 900)
for(j in 1: 2){
for(i in 1: 4){
plot(c(0, 1), c(0, 3), type='n', xaxt='n', yaxt='n', xlab='', ylab='')
mtext(paste('n = ', n[i], sep=''), 1, line=-2.5, cex=0.9, font=3)
}}
par(mar=c(0,3.5,0,0))
plot(c(0, 1), c(0, 3), type='n', xaxt='n', yaxt='n', xlab='', ylab='')
mtext(paste('n = ', n[1], sep=''), 1, line=-2.5, cex=0.9, font=3)
par(mar=c(0,0,0,0))
for(i in 2: 4){
plot(c(0, 1), c(0, 3), type='n', xaxt='n', yaxt='n', xlab='', ylab='')
mtext(paste('n = ', n[i], sep=''), 1, line=-2.5, cex=0.9, font=3)
}
for(i in 1: 3){
plot(c(0, 1), c(0, 3), type='n', xaxt='n', yaxt='n', xlab='', ylab='')
mtext(paste('n = ', n[i], sep=''), 1, line=-2.5, cex=0.9, font=3)
}
par(mar=c(0,0,0,0.5))
plot(c(0, 1), c(0, 3), type='n', xaxt='n', yaxt='n', xlab='', ylab='')
mtext(paste('n = ', n[4], sep=''), 1, line=-2.5, cex=0.9, font=3)

par(mar=c(3.5,0,0,0))
idx = c('MCN', 'MVNC', 'MCNC')
ylim1 = c(0.8, 1)
seq1 = seq(0.8, 1, 0.05)

boxplot(Tab1a2[which(Tab1a2$rate == 0.1 & Tab1a2$criterion == 'CCR'), idx], names=c(1,2,3), ylim=ylim1, las=1, border=c('blue','green4','red'), col=c('cyan','yellow','pink'), cex.axis=1)
boxplot(Tab2a2[which(Tab2a2$rate == 0.1 & Tab2a2$criterion == 'CCR'), idx], names=c(1,2,3), ylim=ylim1, yaxt='n', border=c('blue','green4','red'), col=c('cyan','yellow','pink'), cex.axis=1)
axis(2, seq1, labels=F)
boxplot(Tab3a2[which(Tab3a2$rate == 0.1 & Tab3a2$criterion == 'CCR'), idx], names=c(1,2,3), ylim=ylim1, yaxt='n', border=c('blue','green4','red'), col=c('cyan','yellow','pink'), cex.axis=1)
axis(2, seq1, labels=F)
boxplot(Tab4a2[which(Tab4a2$rate == 0.1 & Tab4a2$criterion == 'CCR'), idx], names=c(1,2,3), ylim=ylim1, yaxt='n', border=c('blue','green4','red'), col=c('cyan','yellow','pink'), cex.axis=1)
axis(2, seq1, labels=F)

boxplot(Tab1a2[which(Tab1a2$rate == 0.2 & Tab1a2$criterion == 'CCR'), idx], names=c(1,2,3), ylim=ylim1, yaxt='n', border=c('blue','green4','red'), col=c('cyan','yellow','pink'), cex.axis=1)
axis(2, seq1, labels=F)
boxplot(Tab2a2[which(Tab2a2$rate == 0.2 & Tab2a2$criterion == 'CCR'), idx], names=c(1,2,3), ylim=ylim1, yaxt='n', border=c('blue','green4','red'), col=c('cyan','yellow','pink'), cex.axis=1)
axis(2, seq1, labels=F)
boxplot(Tab3a2[which(Tab3a2$rate == 0.2 & Tab3a2$criterion == 'CCR'), idx], names=c(1,2,3), ylim=ylim1, yaxt='n', border=c('blue','green4','red'), col=c('cyan','yellow','pink'), cex.axis=1)
boxplot(Tab4a2[which(Tab4a2$rate == 0.2 & Tab4a2$criterion == 'CCR'), idx], names=c(1,2,3), ylim=ylim1, yaxt='n', border=c('blue','green4','red'), col=c('cyan','yellow','pink'), cex.axis=1)
axis(2, seq1, labels=F)

par(mar=c(3.5,3.5,0,0))
ylim2 = ylim1
seq2 = seq1

boxplot(Tab1b2[which(Tab1b2$rate == 0.1 & Tab1b2$criterion == 'CCR'), idx], names=c(1,2,3), ylim=ylim2, las=1, border=c('blue','green4','red'), col=c('cyan','yellow','pink'), cex.axis=1)
par(mar=c(3.5,0,0,0))
boxplot(Tab2b2[which(Tab2b2$rate == 0.1 & Tab2b2$criterion == 'CCR'), idx], names=c(1,2,3), ylim=ylim2, yaxt='n', border=c('blue','green4','red'), col=c('cyan','yellow','pink'), cex.axis=1)
axis(2, seq2, labels=F)
boxplot(Tab3b2[which(Tab3b2$rate == 0.1 & Tab3b2$criterion == 'CCR'), idx], names=c(1,2,3), ylim=ylim2, yaxt='n', border=c('blue','green4','red'), col=c('cyan','yellow','pink'), cex.axis=1)
axis(2, seq2, labels=F)
boxplot(Tab4b2[which(Tab4b2$rate == 0.1 & Tab4b2$criterion == 'CCR'), idx], names=c(1,2,3), ylim=ylim2, yaxt='n', border=c('blue','green4','red'), col=c('cyan','yellow','pink'), cex.axis=1)
axis(2, seq2, labels=F)

boxplot(Tab1b2[which(Tab1b2$rate == 0.2 & Tab1b2$criterion == 'CCR'), idx], names=c(1,2,3), ylim=ylim2, yaxt='n', border=c('blue','green4','red'), col=c('cyan','yellow','pink'), cex.axis=1)
axis(2, seq2, labels=F)
boxplot(Tab2b2[which(Tab2b2$rate == 0.2 & Tab2b2$criterion == 'CCR'), idx], names=c(1,2,3), ylim=ylim2, yaxt='n', border=c('blue','green4','red'), col=c('cyan','yellow','pink'), cex.axis=1)
axis(2, seq2, labels=F)
boxplot(Tab3b2[which(Tab3b2$rate == 0.2 & Tab3b2$criterion == 'CCR'), idx], names=c(1,2,3), ylim=ylim2, yaxt='n', border=c('blue','green4','red'), col=c('cyan','yellow','pink'), cex.axis=1)
par(mar=c(3.5,0,0,0.5))
boxplot(Tab4b2[which(Tab4b2$rate == 0.2 & Tab4b2$criterion == 'CCR'), idx], names=c(1,2,3), ylim=ylim2, yaxt='n', border=c('blue','green4','red'), col=c('cyan','yellow','pink'), cex.axis=1)
axis(2, seq2, labels=F)

# AIR
par(mar=c(0,0.5,0.5,0))
plot(c(0, 1), c(0, 3), type='n', xaxt='n', yaxt='n', xlab='', ylab='')
legend('center', c('1:FM-MCN', '2:FM-MNC', '3:FM-MCNC'), border=c('blue','green4','red'), fill=c('cyan','yellow','pink'), cex=0.9, bty='n')
par(mar=c(3.5,0.5,0,0))
plot(c(0, 1), c(0, 3), type='n', xaxt='n', yaxt='n', xlab='', ylab='')
mtext('ARI', 2, line=-3, cex=1.5, font=2)
par(mar=c(0,0,0.5,0))
plot(c(0, 1), c(0, 3), type='n', xaxt='n', yaxt='n', xlab='', ylab='')
mtext('(b) MCN Censored Data', 1, line=-2, cex=1.5, font=2)
par(mar=c(0,3.5,0.5,0.5))
plot(c(0, 1), c(0, 3), type='n', xaxt='n', yaxt='n', xlab='', ylab='')
mtext('(d) MN Censored Data', 1, line=-2, cex=1.5, font=2)
par(mar=c(0,0,0,0))
plot(c(0, 1), c(0, 3), type='n', xaxt='n', yaxt='n', xlab='', ylab='')
mtext('10% Censoring', 1, line=-2.25, cex=1.2, font=2)
plot(c(0, 1), c(0, 3), type='n', xaxt='n', yaxt='n', xlab='', ylab='')
mtext('20% Censoring', 1, line=-2.25, cex=1.2, font=2)
par(mar=c(0,3.5,0,0))
plot(c(0, 1), c(0, 3), type='n', xaxt='n', yaxt='n', xlab='', ylab='')
mtext('10% Censoring', 1, line=-2.25, cex=1.2, font=2)
par(mar=c(0,0,0,0.5))
plot(c(0, 1), c(0, 3), type='n', xaxt='n', yaxt='n', xlab='', ylab='')
mtext('20% Censoring', 1, line=-2.25, cex=1.2, font=2)
par(mar=c(0,0,0,0))
n = c(150, 300, 600, 900)
for(j in 1: 2){
for(i in 1: 4){
plot(c(0, 1), c(0, 3), type='n', xaxt='n', yaxt='n', xlab='', ylab='')
mtext(paste('n = ', n[i], sep=''), 1, line=-2.5, cex=0.9, font=3)
}}
par(mar=c(0,3.5,0,0))
plot(c(0, 1), c(0, 3), type='n', xaxt='n', yaxt='n', xlab='', ylab='')
mtext(paste('n = ', n[1], sep=''), 1, line=-2.5, cex=0.9, font=3)
par(mar=c(0,0,0,0))
for(i in 2: 4){
plot(c(0, 1), c(0, 3), type='n', xaxt='n', yaxt='n', xlab='', ylab='')
mtext(paste('n = ', n[i], sep=''), 1, line=-2.5, cex=0.9, font=3)
}
for(i in 1: 3){
plot(c(0, 1), c(0, 3), type='n', xaxt='n', yaxt='n', xlab='', ylab='')
mtext(paste('n = ', n[i], sep=''), 1, line=-2.5, cex=0.9, font=3)
}
par(mar=c(0,0,0,0.5))
plot(c(0, 1), c(0, 3), type='n', xaxt='n', yaxt='n', xlab='', ylab='')
mtext(paste('n = ', n[4], sep=''), 1, line=-2.5, cex=0.9, font=3)

par(mar=c(3.5,0,0,0))
ylim3 = c(0.3, 1)
seq3 = seq(0.3, 1, 0.1)

boxplot(Tab1a2[which(Tab1a2$rate == 0.1 & Tab1a2$criterion == 'ARI'), idx], names=c(1,2,3), ylim=ylim3, las=1, border=c('blue','green4','red'), col=c('cyan','yellow','pink'), cex.axis=1)
par(mar=c(3.5,0,0,0))
boxplot(Tab2a2[which(Tab2a2$rate == 0.1 & Tab2a2$criterion == 'ARI'), idx], names=c(1,2,3), ylim=ylim3, yaxt='n', border=c('blue','green4','red'), col=c('cyan','yellow','pink'), cex.axis=1)
axis(2, seq3, labels=F)
boxplot(Tab3a2[which(Tab3a2$rate == 0.1 & Tab3a2$criterion == 'ARI'), idx], names=c(1,2,3), ylim=ylim3, yaxt='n', border=c('blue','green4','red'), col=c('cyan','yellow','pink'), cex.axis=1)
axis(2, seq3, labels=F)
boxplot(Tab4a2[which(Tab4a2$rate == 0.1 & Tab4a2$criterion == 'ARI'), idx], names=c(1,2,3), ylim=ylim3, yaxt='n', border=c('blue','green4','red'), col=c('cyan','yellow','pink'), cex.axis=1)
axis(2, seq3, labels=F)

boxplot(Tab1a2[which(Tab1a2$rate == 0.2 & Tab1a2$criterion == 'ARI'), idx], names=c(1,2,3), ylim=ylim3, yaxt='n', border=c('blue','green4','red'), col=c('cyan','yellow','pink'), cex.axis=1)
axis(2, seq3, labels=F)
boxplot(Tab2a2[which(Tab2a2$rate == 0.2 & Tab2a2$criterion == 'ARI'), idx], names=c(1,2,3), ylim=ylim3, yaxt='n', border=c('blue','green4','red'), col=c('cyan','yellow','pink'), cex.axis=1)
axis(2, seq3, labels=F)
boxplot(Tab3a2[which(Tab3a2$rate == 0.2 & Tab3a2$criterion == 'ARI'), idx], names=c(1,2,3), ylim=ylim3, yaxt='n', border=c('blue','green4','red'), col=c('cyan','yellow','pink'), cex.axis=1)
axis(2, seq3, labels=F)
boxplot(Tab4a2[which(Tab4a2$rate == 0.2 & Tab4a2$criterion == 'ARI'), idx], names=c(1,2,3), ylim=ylim3, yaxt='n', border=c('blue','green4','red'), col=c('cyan','yellow','pink'), cex.axis=1)
axis(2, seq3, labels=F)

par(mar=c(3.5,3.5,0,0))
ylim4 = ylim3
seq4 = seq3

boxplot(Tab1b2[which(Tab1b2$rate == 0.1 & Tab1b2$criterion == 'ARI'), idx], names=c(1,2,3), ylim=ylim4, las=1, border=c('blue','green4','red'), col=c('cyan','yellow','pink'), cex.axis=1)
par(mar=c(3.5,0,0,0))
boxplot(Tab2b2[which(Tab2b2$rate == 0.1 & Tab2b2$criterion == 'ARI'), idx], names=c(1,2,3), ylim=ylim4, yaxt='n', border=c('blue','green4','red'), col=c('cyan','yellow','pink'), cex.axis=1)
axis(2, seq4, labels=F)
boxplot(Tab3b2[which(Tab3b2$rate == 0.1 & Tab3b2$criterion == 'ARI'), idx], names=c(1,2,3), ylim=ylim4, yaxt='n', border=c('blue','green4','red'), col=c('cyan','yellow','pink'), cex.axis=1)
axis(2, seq4, labels=F)
boxplot(Tab4b2[which(Tab4b2$rate == 0.1 & Tab4b2$criterion == 'ARI'), idx], names=c(1,2,3), ylim=ylim4, yaxt='n', border=c('blue','green4','red'), col=c('cyan','yellow','pink'), cex.axis=1)
axis(2, seq4, labels=F)

boxplot(Tab1b2[which(Tab1b2$rate == 0.2 & Tab1b2$criterion == 'ARI'), idx], names=c(1,2,3), ylim=ylim4, yaxt='n', border=c('blue','green4','red'), col=c('cyan','yellow','pink'), cex.axis=1)
axis(2, seq4, labels=F)
boxplot(Tab2b2[which(Tab2b2$rate == 0.2 & Tab2b2$criterion == 'ARI'), idx], names=c(1,2,3), ylim=ylim4, yaxt='n', border=c('blue','green4','red'), col=c('cyan','yellow','pink'), cex.axis=1)
axis(2, seq4, labels=F)
boxplot(Tab3b2[which(Tab3b2$rate == 0.2 & Tab3b2$criterion == 'ARI'), idx], names=c(1,2,3), ylim=ylim4, yaxt='n', border=c('blue','green4','red'), col=c('cyan','yellow','pink'), cex.axis=1)
axis(2, seq4, labels=F)
par(mar=c(3.5,0,0,0.5))
boxplot(Tab4b2[which(Tab4b2$rate == 0.2 & Tab4b2$criterion == 'ARI'), idx], names=c(1,2,3), ylim=ylim4, yaxt='n', border=c('blue','green4','red'), col=c('cyan','yellow','pink'), cex.axis=1)
axis(2, seq4, labels=F)
dev.off()
