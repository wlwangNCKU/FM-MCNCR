##################################################################################
#
#   Filename: TableC4.r
#   Purpose: produce Table C.4 for averages of AIC and BIC scores and frequencies 
#            in Experiment 1
#   Input data files: intermediate results (modelfit.txt) sorted in the subfolders 
#                     'results/Experiment2/SIM1/'; 'results/Experiment2/SIM2/'; 
#                     'results/Experiment2/SIM3/'; 'results/Experiment2/SIM4/'; 
#   Output data files: results/Experiment2/TableC4.csv
#
##################################################################################

PATH1 = paste(PATH, 'results/Experiment2/SIM1/', sep='')
PATH2 = paste(PATH, 'results/Experiment2/SIM2/', sep='')
PATH3 = paste(PATH, 'results/Experiment2/SIM3/', sep='')
PATH4 = paste(PATH, 'results/Experiment2/SIM4/', sep='')

### Read outputs ###
Tab11 = read.table(paste(PATH1, 'modelfit.txt',sep=""))
Tab21 = read.table(paste(PATH2, 'modelfit.txt',sep=""))
Tab31 = read.table(paste(PATH3, 'modelfit.txt',sep=""))
Tab41 = read.table(paste(PATH4, 'modelfit.txt',sep=""))

name1 = c('Rep','rate','criterion','g1','g2','g3','g4','g1','g2','g3','g4')
names(Tab11) = names(Tab21) = names(Tab31) = names(Tab41) = name1

### model selection ###
selection = function(cri=c('aic','bic'), Tab1, Tab2, Tab3, Tab4)
{
cri = cri[1]
A1 = Tab1[which(Tab1$rate == 0.1 & Tab1$criterion == cri), 4:7]
AC1 = Tab1[which(Tab1$rate == 0.1 & Tab1$criterion == cri), 8:11]
A2 = Tab2[which(Tab2$rate == 0.1 & Tab2$criterion == cri), 4:7]
AC2 = Tab2[which(Tab2$rate == 0.1 & Tab2$criterion == cri), 8:11]
A3 = Tab3[which(Tab3$rate == 0.1 & Tab3$criterion == cri), 4:7]
AC3 = Tab3[which(Tab3$rate == 0.1 & Tab3$criterion == cri), 8:11]
A4 = Tab4[which(Tab4$rate == 0.1 & Tab4$criterion == cri), 4:7]
AC4 = Tab4[which(Tab4$rate == 0.1 & Tab4$criterion == cri), 8:11]

B1 = Tab1[which(Tab1$rate == 0.2 & Tab1$criterion == cri), 4:7]
BC1 = Tab1[which(Tab1$rate == 0.2 & Tab1$criterion == cri), 8:11]
B2 = Tab2[which(Tab2$rate == 0.2 & Tab2$criterion == cri), 4:7]
BC2 = Tab2[which(Tab2$rate == 0.2 & Tab2$criterion == cri), 8:11]
B3 = Tab3[which(Tab3$rate == 0.2 & Tab3$criterion == cri), 4:7]
BC3 = Tab3[which(Tab3$rate == 0.2 & Tab3$criterion == cri), 8:11]
B4 = Tab4[which(Tab4$rate == 0.2 & Tab4$criterion == cri), 4:7]
BC4 = Tab4[which(Tab4$rate == 0.2 & Tab4$criterion == cri), 8:11]

out1 = rbind(t(cbind(colMeans(A1), colMeans(AC1), colMeans(A2), colMeans(AC2), colMeans(A3), colMeans(AC3), colMeans(A4), colMeans(AC4))),
             t(cbind(colMeans(B1), colMeans(BC1), colMeans(B2), colMeans(BC2), colMeans(B3), colMeans(BC3), colMeans(B4), colMeans(BC4))))

a1 = table(matrix(apply(A1, 1, order), nrow=4)[1,])
a2 = table(matrix(apply(A2, 1, order), nrow=4)[1,])
a3 = table(matrix(apply(A3, 1, order), nrow=4)[1,])
a4 = table(matrix(apply(A4, 1, order), nrow=4)[1,])

ac1 = table(matrix(apply(AC1, 1, order), nrow=4)[1,])
ac2 = table(matrix(apply(AC2, 1, order), nrow=4)[1,])
ac3 = table(matrix(apply(AC3, 1, order), nrow=4)[1,])
ac4 = table(matrix(apply(AC4, 1, order), nrow=4)[1,])

b1 = table(matrix(apply(B1, 1, order), nrow=4)[1,])
b2 = table(matrix(apply(B2, 1, order), nrow=4)[1,])
b3 = table(matrix(apply(B3, 1, order), nrow=4)[1,])
b4 = table(matrix(apply(B4, 1, order), nrow=4)[1,])

bc1 = table(matrix(apply(BC1, 1, order), nrow=4)[1,])
bc2 = table(matrix(apply(BC2, 1, order), nrow=4)[1,])
bc3 = table(matrix(apply(BC3, 1, order), nrow=4)[1,])
bc4 = table(matrix(apply(BC4, 1, order), nrow=4)[1,])

n11 = nc11 = n12 = nc12 = n13 = nc13 = n14 = nc14 =
n21 = nc21 = n22 = nc22 = n23 = nc23 = n24 = nc24 = numeric(4)
for(i in 1: 4){
 if(length(which(names(a1)==i)) != 0) n11[i] = a1[which(names(a1)==i)]
 if(length(which(names(a2)==i)) != 0) n12[i] = a2[which(names(a2)==i)]
 if(length(which(names(a3)==i)) != 0) n13[i] = a3[which(names(a3)==i)]
 if(length(which(names(a4)==i)) != 0) n14[i] = a4[which(names(a4)==i)]

 if(length(which(names(ac1)==i)) != 0) nc11[i] = ac1[which(names(ac1)==i)]
 if(length(which(names(ac2)==i)) != 0) nc12[i] = ac2[which(names(ac2)==i)]
 if(length(which(names(ac3)==i)) != 0) nc13[i] = ac3[which(names(ac3)==i)]
 if(length(which(names(ac4)==i)) != 0) nc14[i] = ac4[which(names(ac4)==i)]
}

for(i in 1: 4){
 if(length(which(names(b1)==i)) != 0) n21[i] = b1[which(names(b1)==i)]
 if(length(which(names(b2)==i)) != 0) n22[i] = b2[which(names(b2)==i)]
 if(length(which(names(b3)==i)) != 0) n23[i] = b3[which(names(b3)==i)]
 if(length(which(names(b4)==i)) != 0) n24[i] = b4[which(names(b4)==i)]

 if(length(which(names(bc1)==i)) != 0) nc21[i] = bc1[which(names(bc1)==i)]
 if(length(which(names(bc2)==i)) != 0) nc22[i] = bc2[which(names(bc2)==i)]
 if(length(which(names(bc3)==i)) != 0) nc23[i] = bc3[which(names(bc3)==i)]
 if(length(which(names(bc4)==i)) != 0) nc24[i] = bc4[which(names(bc4)==i)]
}

out2 = rbind(t(cbind(n11, nc11, n12, nc12, n13, nc13, n14, nc14)),
             t(cbind(n21, nc21, n22, nc22, n23, nc23, n24, nc24)))
return(list(Mean = round(out1, 2), Freq = out2))
}

# AIC
TabAIC = selection(cri='aic', Tab1=Tab11, Tab2=Tab21, Tab3=Tab31, Tab4=Tab41)
# BIC
TabBIC = selection(cri='bic', Tab1=Tab11, Tab2=Tab21, Tab3=Tab31, Tab4=Tab41)

TabC4 = NULL
for(i in 1: 16){
TabC4 = rbind(TabC4, cbind(rbind(TabAIC$Mean[i, ], TabAIC$Freq[i, ]), 
                           rbind(TabBIC$Mean[i, ], TabBIC$Freq[i, ])))
}

TabC4 = data.frame(rep(c('10%', '20%'), each=16), 
                   rep(rep(c(150, 300, 600, 900), each=4), 2), 
                   rep(rep(c('FM-MCNR', 'FM-MCNCR'), each=2), 8), TabC4)

colnames(TabC4) = c('Crate', 'Size', 'Model', 'AIC:g1', 'AIC:g2', 'AIC:g3', 'AIC:g4',
                    'BIC:g1', 'BIC:g2', 'BIC:g3', 'BIC:g4')

print(round(TabC4, 3))

write.csv(TabC4, paste(PATH, 'results/Experiment2/TableC4.csv',sep=""), row.names=FALSE)
