##################################################################################
#
#   Filename: TableC1.r
#   Purpose: produce Table C.1 for averages of AIC and BIC scores and frequencies 
#            in Experiment 1
#   Input data files: intermediate results (modelfit.txt) sorted in the subfolders 
#                    'results/Experiment1/SIM1a/'; 'results/Experiment1/SIM1b/'; 
#                    'results/Experiment1/SIM2a/'; 'results/Experiment1/SIM2b/'; 
#                    'results/Experiment1/SIM3a/'; 'results/Experiment1/SIM3b/'; 
#                    'results/Experiment1/SIM4a/'; 'results/Experiment1/SIM4b/'; 
#   Output data files: results/Experiment1/TableC1.csv
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
Tab1a1 = read.table(paste(PATH1a, 'modelfit.txt',sep=""))
Tab2a1 = read.table(paste(PATH2a, 'modelfit.txt',sep=""))
Tab3a1 = read.table(paste(PATH3a, 'modelfit.txt',sep=""))
Tab4a1 = read.table(paste(PATH4a, 'modelfit.txt',sep=""))

# MVN #
Tab1b1 = read.table(paste(PATH1b, 'modelfit.txt',sep=""))
Tab2b1 = read.table(paste(PATH2b, 'modelfit.txt',sep=""))
Tab3b1 = read.table(paste(PATH3b, 'modelfit.txt',sep=""))
Tab4b1 = read.table(paste(PATH4b, 'modelfit.txt',sep=""))

name1 = c('Rep', 'rate', 'criterion', 'MCN', 'MVNC', 'MCNC')

colnames(Tab1a1) = colnames(Tab1b1) = colnames(Tab2a1) = colnames(Tab2b1) = 
colnames(Tab3a1) = colnames(Tab3b1) = colnames(Tab4a1) = colnames(Tab4b1) = name1

selection = function(cri=c('aic','bic'), Tab1, Tab2, Tab3, Tab4)
{
cri = cri[1]
A1 = Tab1[which(Tab1$rate == 0.1 & Tab1$criterion == cri), 5:6]
A2 = Tab2[which(Tab2$rate == 0.1 & Tab2$criterion == cri), 5:6]
A3 = Tab3[which(Tab3$rate == 0.1 & Tab3$criterion == cri), 5:6]
A4 = Tab4[which(Tab4$rate == 0.1 & Tab4$criterion == cri), 5:6]
B1 = Tab1[which(Tab1$rate == 0.2 & Tab1$criterion == cri), 5:6]
B2 = Tab2[which(Tab2$rate == 0.2 & Tab2$criterion == cri), 5:6]
B3 = Tab3[which(Tab3$rate == 0.2 & Tab3$criterion == cri), 5:6]
B4 = Tab4[which(Tab4$rate == 0.2 & Tab4$criterion == cri), 5:6]

out1 = cbind(colMeans(A1), colMeans(A2), colMeans(A3), colMeans(A4),
      colMeans(B1), colMeans(B2), colMeans(B3), colMeans(B4))

a1 = table(matrix(apply(A1, 1, order), nrow=2)[1,])
a2 = table(matrix(apply(A2, 1, order), nrow=2)[1,])
a3 = table(matrix(apply(A3, 1, order), nrow=2)[1,])
a4 = table(matrix(apply(A4, 1, order), nrow=2)[1,])

b1 = table(matrix(apply(B1, 1, order), nrow=2)[1,])
b2 = table(matrix(apply(B2, 1, order), nrow=2)[1,])
b3 = table(matrix(apply(B3, 1, order), nrow=2)[1,])
b4 = table(matrix(apply(B4, 1, order), nrow=2)[1,])

n11 = n12 = n13 = n14 = n21 = n22 = n23 = n24 = numeric(2)
for(i in 1: 2){
 if(length(which(names(a1)==i)) != 0) n11[i] = a1[which(names(a1)==i)]
 if(length(which(names(a2)==i)) != 0) n12[i] = a2[which(names(a2)==i)]
 if(length(which(names(a3)==i)) != 0) n13[i] = a3[which(names(a3)==i)]
 if(length(which(names(a4)==i)) != 0) n14[i] = a4[which(names(a4)==i)]

 if(length(which(names(b1)==i)) != 0) n21[i] = b1[which(names(b1)==i)]
 if(length(which(names(b2)==i)) != 0) n22[i] = b2[which(names(b2)==i)]
 if(length(which(names(b3)==i)) != 0) n23[i] = b3[which(names(b3)==i)]
 if(length(which(names(b4)==i)) != 0) n24[i] = b4[which(names(b4)==i)]
}

out2 = cbind(n11, n12, n13, n14, n21, n22, n23, n24)
return(list(Mean = round(out1, 2), Freq = out2))
}

MCNaic = selection(cri='aic', Tab1 = Tab1a1, Tab2 = Tab2a1, Tab3 = Tab3a1, Tab4 = Tab4a1)
MCNbic = selection(cri='bic', Tab1 = Tab1a1, Tab2 = Tab2a1, Tab3 = Tab3a1, Tab4 = Tab4a1)
MNaic  = selection(cri='aic', Tab1 = Tab1b1, Tab2 = Tab2b1, Tab3 = Tab3b1, Tab4 = Tab4b1)
MNbic  = selection(cri='bic', Tab1 = Tab1b1, Tab2 = Tab2b1, Tab3 = Tab3b1, Tab4 = Tab4b1)

TabC1 = data.frame(rep(c('MCN', 'MN'), each=8), rep(rep(c('AIC', 'BIC'), each=4), 2),
        rep(rep(c('FM-MNC', 'FM-MCNC'), each=2), 4),
        rbind(MCNaic$Mean[1,], MCNaic$Freq[1,], MCNaic$Mean[2,], MCNaic$Freq[2,],
              MCNbic$Mean[1,], MCNbic$Freq[1,], MCNbic$Mean[2,], MCNbic$Freq[2,],
              MNaic$Mean[1,],  MNaic$Freq[1,],  MNaic$Mean[2,],  MNaic$Freq[2,],
              MNbic$Mean[1,],  MNbic$Freq[1,],  MNbic$Mean[2,],  MNbic$Freq[2,]))

colnames(TabC1) = c('Data', 'Criteria', 'Model', 
                    'C10%-n150', 'C10%-n300', 'C10%-n600', 'C10%-n900',
                    'C20%-n150', 'C20%-n300', 'C20%-n600', 'C20%-n900')

write.csv(TabC1, paste(PATH, 'results/Experiment1/TableC1.csv',sep=""), row.names=FALSE)
