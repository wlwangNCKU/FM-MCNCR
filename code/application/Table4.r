################################################################################
#
#   Filename: Table4.r
#   Purpose: produce Table 4 for the Mroz87 data
#   Input data files: StdWageFit.RData
#   Output data files: results/Table4.csv
#
################################################################################

load(paste(PATH, 'data/StdWageFit.RData', sep=''))

Tab4 = cbind(t(estNC2$IM$theta.hat)[1:25, ],
       rbind(t(estTC2$IM$theta.hat)[1:24, ], 0),
             t(estCC2$IM$theta.hat)[1:25, ],
       rbind(c(estNC2$para$w[2],0), t(estNC2$IM$theta.hat)[26:49, ]),
       rbind(c(estTC2$para$w[2],0), t(estTC2$IM$theta.hat)[25:47, ], 0),
       rbind(c(estCC2$para$w[2],0), t(estCC2$IM$theta.hat)[26:49, ]))

colnames(Tab4) = c('FM-MNCR(cls1)', 'SD', 'FM-MTCR(cls1)', 'SD', 'FM-MCNCR(cls1)', 'SD',
                   'FM-MNCR(cls2)', 'SD', 'FM-MTCR(cls2)', 'SD', 'FM-MCNCR(cls2)', 'SD')
rownames(Tab4) = c('w1', 'beta.i11', 'beta.i12', 'beta.i13', 'beta.i21', 'beta.i22', 'beta.i23',
                   'beta.i31', 'beta.i32', 'beta.i33', 'beta.i41', 'beta.i42', 'beta.i43',
                   'sig.i11', 'sig.i21', 'sig.i22', 'sig.i31', 'sig.i32', 'sig.i33', 
                   'sig.i41', 'sig.i42', 'sig.i43', 'sig.i44', 'nui', 'rhoi')

print(round(Tab4, 3))

PATH = paste(getwd(),"/Data_and_Code_FMMCNCR/",sep="")
write.csv(Tab4, paste(PATH, 'results/Table4.csv',sep=""))
