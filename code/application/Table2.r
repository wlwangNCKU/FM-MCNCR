################################################################################
#
#   Filename: Table2.r
#   Purpose: produce Table 2 for the Wholesale data
#   Input data files: data/SaleFit.RData
#   Output data files: results/Table2.csv
#
################################################################################

load(paste(PATH, 'data/SaleFit.RData', sep=''))

Tab2 = cbind(t(estNC2$IM$theta), 
             t(cbind(estTC2$IM$theta[,1:16], 0, estTC2$IM$theta[,17:31], 0)), 
             t(estCC2$IM$theta))
colnames(Tab2) = c('FM-MNC', 'SD', 'FM-MTC', 'SD', 'FM-MCNC', 'SD')
rownames(Tab2) = c('w1', 'mu11', 'mu12', 'mu13', 'mu14', 
                   'sig111', 'sig121', 'sig122', 'sig131', 'sig132', 'sig133', 
                   'sig141', 'sig142', 'sig143', 'sig144', 'nu1', 'rho1',
                   'mu21', 'mu22', 'mu23', 'mu24', 
                   'sig211', 'sig221', 'sig222', 'sig231', 'sig232', 'sig233', 
                   'sig241', 'sig242', 'sig243', 'sig244', 'nu2', 'rho2')
print(round(Tab2, 3))

PATH = paste(getwd(),"/Data_and_Code_FMMCNCR/",sep="")
write.csv(Tab2, paste(PATH, 'results/Table2.csv',sep=""))
