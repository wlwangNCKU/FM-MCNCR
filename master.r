rm(list = ls())
PATH = paste(getwd(),"/Data_and_Code_FMMCNCR/",sep="")

### Application ### 
## Example 1 ##
# Re-perform the intermediate results for Figure 1, Table 1 and Table 2 
source(paste(PATH, 'code/application/runsale.r', sep=''))

# Re-produce Figure 1
source(paste(PATH, 'code/application/fig1.r', sep=''))

# Re-produce Table 1
source(paste(PATH, 'code/application/Table1.r', sep=''))

# Re-produce Table 2
source(paste(PATH, 'code/application/Table2.r', sep=''))


## Example 2 ##
# Re-perform the intermediate results for Figure 2, Table 3 and Table 4
source(paste(PATH, 'code/application/runwage.r', sep=''))

# Re-produce Figure 2
source(paste(PATH, 'code/application/fig2.r', sep=''))

# Re-produce Table 3
source(paste(PATH, 'code/application/Table3.r', sep=''))

# Re-produce Table 4
source(paste(PATH, 'code/application/Table4.r', sep=''))


### Simulation ###
## Experiment 1 ##
# Re-generate simulation results for Experiment 1 with MCN data of size n=150
n = 150; DIS = 'MCN'; WD.PATH = paste(PATH, 'results/Experiment1/SIM1a/', sep='')
source(paste(PATH, 'code/simulation/simMCN.r', sep=''))

# Re-generate simulation results for Experiment 1 with MVN data of size n=150
n = 150; DIS = 'MVN'; WD.PATH = paste(PATH, 'results/Experiment1/SIM1b/', sep='')
source(paste(PATH, 'code/simulation/simMCN.r', sep=''))

# Re-generate simulation results for Experiment 1 with MCN data of size n=300
n = 300; DIS = 'MCN'; WD.PATH = paste(PATH, 'results/Experiment1/SIM2a/', sep='')
source(paste(PATH, 'code/simulation/simMCN.r', sep=''))

# Re-generate simulation results for Experiment 1 with MVN data of size n=300
n = 300; DIS = 'MVN'; WD.PATH = paste(PATH, 'results/Experiment1/SIM2b/', sep='')
source(paste(PATH, 'code/simulation/simMCN.r', sep=''))

# Re-generate simulation results for Experiment 1 with MCN data of size n=600
n = 600; DIS = 'MCN'; WD.PATH = paste(PATH, 'results/Experiment1/SIM3a/', sep='')
source(paste(PATH, 'code/simulation/simMCN.r', sep=''))

# Re-generate simulation results for Experiment 1 with MVN data of size n=600
n = 600; DIS = 'MVN'; WD.PATH = paste(PATH, 'results/Experiment1/SIM3b/', sep='')
source(paste(PATH, 'code/simulation/simMCN.r', sep=''))

# Re-generate simulation results for Experiment 1 with MCN data of size n=900
n = 900; DIS = 'MCN'; WD.PATH = paste(PATH, 'results/Experiment1/SIM4a/', sep='')
source(paste(PATH, 'code/simulation/simMCN.r', sep=''))

# Re-generate simulation results for Experiment 1 with MVN data of size n=900
n = 900; DIS = 'MVN'; WD.PATH = paste(PATH, 'results/Experiment1/SIM4b/', sep='')
source(paste(PATH, 'code/simulation/simMCN.r', sep=''))

# Re-produce Figure C.1
source(paste(PATH, 'code/simulation/figC1.r', sep=''))

# Re-produce Figure C.2
source(paste(PATH, 'code/simulation/figC2.r', sep=''))

# Re-produce Figure C.3a and Figure C.3b
source(paste(PATH, 'code/simulation/figC3ab.r', sep=''))

# Re-produce Table C.1
source(paste(PATH, 'code/simulation/TableC1.r', sep=''))

# Re-produce Table C.2 and Table C.3
source(paste(PATH, 'code/simulation/TablesC2C3.r', sep=''))


## Experiment 2 ##
# Re-generate simulation results for Experiment 2 with sample size n=150
n = 150; WD.PATH = paste(PATH, 'results/Experiment2/SIM1/', sep='')
source(paste(PATH, 'code/simulation/simMCNR.r', sep=''))

# Re-generate simulation results for Experiment 2 with sample size n=300
n = 300; WD.PATH = paste(PATH, 'results/Experiment2/SIM2/', sep='')
source(paste(PATH, 'code/simulation/simMCNR.r', sep=''))

# Re-generate simulation results for Experiment 2 with sample size n=600
n = 600; WD.PATH = paste(WD.PATH, 'results/Experiment2/SIM3/', sep='')
source(paste(PATH, 'code/simulation/simMCNR.r', sep=''))

# Re-generate simulation results for Experiment 2 with sample size n=900
n = 900; WD.PATH = paste(PATH, 'results/Experiment2/SIM4/', sep='')
source(paste(PATH, 'code/simulation/simMCNR.r', sep=''))

# Re-produce Figure C.4
source(paste(PATH, 'code/simulation/figC4.r', sep=''))

# Re-produce Figure C.5
source(paste(PATH, 'code/simulation/figC5.r', sep=''))

# Re-produce Table C.4
source(paste(PATH, 'code/simulation/TableC4.r', sep=''))

# Re-produce Table C.5
source(paste(PATH, 'code/simulation/TableC5.r', sep=''))
