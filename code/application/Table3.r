################################################################################
#
#   Filename: Table3.r
#   Purpose: produce Table 3 for the Mroz87 data
#   Input data files: StdWageFit.RData
#   Output data files: results/Table3.csv
#
################################################################################

load(paste(PATH, 'data/StdWageFit.RData', sep=''))

Tab3 = rbind(as.vector(rbind(c(estNC1$model.inf$m, estNC2$model.inf$m, estNC3$model.inf$m),
                c(estTC1$model.inf$m, estTC2$model.inf$m, estTC3$model.inf$m),
                c(estCC1$model.inf$m, estCC2$model.inf$m, estCC3$model.inf$m))),
as.vector(rbind(c(estNC1$model.inf$loglik, estNC2$model.inf$loglik, estNC3$model.inf$loglik),
                c(estTC1$model.inf$loglik, estTC2$model.inf$loglik, estTC3$model.inf$loglik),
                c(estCC1$model.inf$loglik, estCC2$model.inf$loglik, estCC3$model.inf$loglik))),
as.vector(rbind(c(estNC1$model.inf$aic, estNC2$model.inf$aic, estNC3$model.inf$aic),
                c(estTC1$model.inf$aic, estTC2$model.inf$aic, estTC3$model.inf$aic),
                c(estCC1$model.inf$aic, estCC2$model.inf$aic, estCC3$model.inf$aic))),
as.vector(rbind(c(estNC1$model.inf$bic, estNC2$model.inf$bic, estNC3$model.inf$bic),
                c(estTC1$model.inf$bic, estTC2$model.inf$bic, estTC3$model.inf$bic),
                c(estCC1$model.inf$bic, estCC2$model.inf$bic, estCC3$model.inf$bic))),
as.vector(rbind(c(estNC1$model.inf$run, estNC2$model.inf$run, estNC3$model.inf$run),
                c(estTC1$model.inf$run, estTC2$model.inf$run, estTC3$model.inf$run),
                c(estCC1$model.inf$run, estCC2$model.inf$run, estCC3$model.inf$run)))) 

colnames(Tab3) = c('MNCR-g1', 'MTCR-g1', 'MCNCR-g1', 
                   'MNCR-g2', 'MTCR-g2', 'MCNCR-g2', 
                   'MNCR-g3', 'MTCR-g3', 'MCNCR-g3')
rownames(Tab3) = c('m', 'loglik', 'AIC', 'BIC', 'CPU(sec)')
print(round(Tab3, 2))

PATH = paste(getwd(),"/Data_and_Code_FMMCNCR/",sep="")
write.csv(Tab3, paste(PATH, 'results/Table3.csv',sep=""))
