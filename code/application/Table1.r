################################################################################
#
#   Filename: Table1.r
#   Purpose: produce Table 1 for the Wholesale data
#   Input data files: data/SaleFit.RData
#   Output data files: results/Table1a.csv, results/Table1b.csv
#
################################################################################

load(paste(PATH, 'data/SaleFit.RData', sep=''))

Tab1a = rbind(c(estNC2$model.inf$m,  estTC2$model.inf$m,  estCC2$model.inf$m),
      c(estNC2$model.inf$loglik, estTC2$model.inf$loglik, estCC2$model.inf$loglik),
      c(estNC2$model.inf$aic,    estTC2$model.inf$aic,    estCC2$model.inf$aic),
      c(estNC2$model.inf$bic,    estTC2$model.inf$bic,    estCC2$model.inf$bic),
      c(estNC2$model.inf$run,    estTC2$model.inf$run,    estCC2$model.inf$run),
      c(CCRnc2, CCRtc2, CCRcc2),
      c(ARInc2, ARItc2, ARIcc2))
colnames(Tab1a) = c('FM-MNC', 'FM-MTC', 'FM-MCNC')
rownames(Tab1a) = c('m', 'loglik', 'AIC', 'BIC', 'CPU(sec)', 'CCR', 'ARI')
print(Tab1a)

PATH = paste(getwd(),"/Data_and_Code_FMMCNCR/",sep="")
write.csv(Tab1a, paste(PATH, 'results/Table1a.csv',sep=""))

N1 = table(factor(estNC2$pre.cls$post.cls[which(cls==1)], level=c(1,2)))
N2 = table(factor(estNC2$pre.cls$post.cls[which(cls==2)], level=c(1,2)))
OUTNC2 = rbind(N1, N2)
colnames(OUTNC2)=1:2
       
T1 = table(factor(estTC2$pre.cls$post.cls[which(cls==1)], level=c(1,2)))
T2 = table(factor(estTC2$pre.cls$post.cls[which(cls==2)], level=c(1,2)))
OUTTC2 = rbind(T1, T2)
colnames(OUTTC2)=1:2
       
C1 = table(factor(estCC2$pre.cls$post.cls[which(cls==1)], level=c(1,2)))
C2 = table(factor(estCC2$pre.cls$post.cls[which(cls==2)], level=c(1,2)))
OUTCC2 = rbind(C1, C2)
colnames(OUTCC2)=1:2

Tab1b = cbind(OUTNC2, OUTTC2, OUTCC2)
rownames(Tab1b) = c('Retail', 'Horeca')
print(Tab1b)
write.csv(Tab1b, paste(PATH, 'results/Table1b.csv',sep=""))
