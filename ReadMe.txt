README

#################################################################################

Source code and data for the manuscript 
"Finite Mixtures of Multivariate Contaminated Normal Censored Regression Model"
by Tsung-I Lin and Wan-Lun Wang*

#################################################################################

# Author responsible for the code #
For questions, comments, or remarks about the code please contact the responsible author, Wan-Lun Wang (wangwl@gs.ncku.edu.tw).

# Configurations #
The code was written/evaluated in R with the following software versions:
R version 4.3.1 (2023-06-16 ucrt)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 19045)

Matrix products: default
locale:
[1] LC_COLLATE=Chinese (Traditional)_Taiwan.utf8  LC_CTYPE=Chinese (Traditional)_Taiwan.utf8    LC_MONETARY=Chinese (Traditional)_Taiwan.utf8
[4] LC_NUMERIC=C                                  LC_TIME=Chinese (Traditional)_Taiwan.utf8    

time zone: Asia/Taipei
tzcode source: internal

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

loaded via a namespace (and not attached):
[1] compiler_4.3.1

# Descriptions of the codes # 
Please extract the files to the "current working directory" of the R package.
The getwd() function shall determine an absolute pathname of the "current working directory".

Before running the codes 'runsale.r', 'runwage.r' and 'TableC1.r', one needs to install the following R packages 

    install.packages("combinat") Version：0.0-8
    install.packages("cubature") Version：2.0.4.2
    install.packages("mvtnorm") Version：1.1-2
    install.packages("sampleSelection")  Version: 1.2-12
    install.packages("tmvtnorm") Version：1.4-10
    install.packages("vioplot") Version: 0.4.0

R codes for the implementation of our methodology are provided.

## Subfolder: ./function ##
./function
	contains the program (function) of 
 
    (1) 'TMmoment.r' for evaluating the first two moments of truncated multivariate normal (MVN) and multivariate contaminated normal (MCN) distributions;
    (2) 'FMMCN.EM.r' for running the ECM algorithm for fitting the FM-MCN model;
    (3) 'FMMCNC.EM.r' for running the ECM algorithm for fitting the FM-MCNC model with left censoring;
    (4) 'FMMCNCint.EM.r' for running the ECM algorithm for fitting the FM-MCNC model with interval censoring;
    (5) 'FMMVTC.EM.r' for running the ECM algorithm for fitting the FM-MVTC model with left censoring;
    (6) 'FMMVTCint.EM.r' for running the ECM algorithm for fitting the FM-MVTC model with interval censoring;
    (7) 'FMMCNR.EM.r' for running the ECM algorithm for fitting the FM-MCNR model;
    (8) 'FMMCNCR.EM.r' for running the ECM algorithm for fitting the FM-MCNCR model with left censoring;
    (9) 'FMMCNCRint.EM.r' for running the ECM algorithm for fitting the FM-MCNCR model with interval censoring;
    (10) 'FMMVTCR.EM.r' for running the ECM algorithm for fitting the FM-MVTCR model with left censoring;
    (11) 'FMMCNCRint.EM.r' for running the ECM algorithm for fitting the FM-MCNCR model with interval censoring;
    (12) 'mclust.fn.r' for calculating the correct classication rate (CCR) and the adjusted Rand index (ARI) for classification.

## Subfolder: ./code ##
./code 
       contains two subfolders
       
  "./application" which collects

    (1) 'runsale.r' main script for re-fitting the FM-MNC, FM-MTC, and FM-MCNC models to the Wholesale dataset;
    (2) 'fig1.r' main script for reproducing Figure 1 in the 'Wholesale data' example;
    (3) 'Table1.r' main script for reproducing Table 1 in the 'Wholesale data' example;
    (4) 'Table2.r' main script for reproducing Table 2 in the 'Wholesale data' example;
    (5) 'runwage.r' main script for re-fitting the FM-MNCR, FM-MTCR, and FM-MCNCR models with g=1~3 components to the Mroz87 data;
    (6) 'fig2.r' main script for reproducing Figure 2 in the 'Wage' example;
    (7) 'Table3.r' main script for reproducing Table 3 in the 'Wage' example;
    (8) 'Table4.r' main script for reproducing Table 4 in the 'Wage' example.

 "./simulation" which collects
 
    (1) 'simMCN.r' main script for re-generate simulation results for Experiment 1 with MCN and MN datasets of various sample sizes;
    (2) 'figC1.r' main script for reproducing Figure C.1 that shows the censored-data recovery in Experiment 1;
    (3) 'figC2.r' main script for reproducing Figure C.2 that shows the classification performance in Experiment 1;
    (4) 'figC3ab.r' main script for reproducing Figure C.3 that shows the mean squared errors (MSE) for parameter estimates across sample sizes in Experiment 1;
    (5) 'TableC1.r' main script for reproducing Table C.1 that reports the results of model selection in Experiment 1;
    (6) 'TablesC2C3.r' main script for reproducing Tables C.2 and C.3 that report Monte Carlo standard deivations (SD) and averages of information-based standard errors (SE) of estimated parameters for the 'MCN' and 'MN' sample data, respectively, in Experiment 1;

    (7) 'simMCNR.r' main script for re-generate simulation results for Experiment 2 under various sample sizes;
    (8) 'figC4.r' main script for reproducing Figure C.4 that shows split violin plots of the MSE scores for parameter estimates and fitted values in Experiment 2;
    (9) 'figC5.r' main script for reproducing Figure C.5 that shows the MSE scores for estimates of each parameter obtained by the 3-component FM-MCNCR model in Experiment 2;
    (10) 'TableC4.r' main script for reproducing Table C.4 that reports the results of model selection in Experiment 2;
    (11) 'TableC5.r' main script for reproducing Table C.5 that report Monte Carlo SD and averages of information-based SE of estimated parameters obtained by the 3-component FM-MCNCR model over 100 replications in Experiment 2.

## Subfolder: ./data ##
./data
      contains
      
      (1) 'Wholesale.RData': the wholesale dataset;
      (2) 'SaleFit.RData': the fitting results of the considered three models to the wholesale data;
      (3) 'mroz87.csv': the U.S. women's labor force participation dataset;
      (4) 'StdWageFit.RData': the fitting results of the FM-MNCR, FM-MTCR, and FM-MCNCR models with g=1~3 components to the 'mroz87' dataset.

## Subfolder: ./results ##
./results
      contains 
      
      (1) 'fig1.eps': visualization plot for the wholesale data and the fitted results based on the 2-component FM-MCNC model;
      (2) 'fig2.eps': visualization plot for the Mroz87 data and the fitted results based on the 2-component FM-MCNCR model;
      (3) 'Table1a.csv': summary of model selection criteria under the fitted 2-component FM-MNC, FM-MTC, and FM-MCNC models for the wholesale data;
      (4) 'Table1b.csv': summary of true groups against classification results under the fitted 2-component FM-MNC, FM-MTC, and FM-MCNC models for the wholesale data;
      (5) 'Table2.csv': estimation results of parameters under the fitted 2-component FM-MNC, FM-MTC, and FM-MCNC models for the wholesale data;
      (6) 'Table3.csv': summary of model selection criteria for the nine candidate models to the Mroz87 data;
      (7) 'Table4.csv': estimation results of parameters under the fitted 2-component FM-MNCR, FM-MTCR, and FM-MCNCR models for the Mroz87 data;

and two subfolders

"./Experiment1" which collects the intermediately numerical results for Experiment 1:

      (1) subsubfolders "./SIM1a", "./SIM1b", "./SIM2a", "./SIM2b", "./SIM3a", "./SIM3b", "./SIM4a", and "./SIM4b",  storing the 'cls.txt', 'estSD.txt', 'modelfit.txt', 'MSEest.txt', and 'MSEyc.txt' text files;
      (2) 'figC1.eps': split violin plots of the MSE for the recovered censored measurements obtained by fitting the 2-component FM-MNC and FM-MCNC models to the generated MCN and MN datasets;
      (3) 'figC2.eps': boxplots for the CCR and ARI scores obtained by the fitted 2-component FM-MCN, FM-MNC and FM-MCNC models to the generated MCN and MN datasets;
      (4) 'figC3a.eps' & 'fig3b.eps': plots of MSE for the estimates of model parameters obtained by the fitted 2-component FM-MCN, FM-MNC and FM-MCNC models to the generated MCN and MN datasets, respectively;
      (5) 'TableC1.csv': averages of AIC and BIC scores and frequencies supported by the two criteria for the FM-MNC and FM-MCNC models fitted to the simulated datasets under all scenarios considered;
      (6) 'TableC2.csv': Monte Carlo SD and averages of information-based SE of estimated parameters obtained by fitting the FM-MCNC model to the MCN datasets;
      (7) 'TableC3.csv': Monte Carlo SD and averages of information-based SE of estimated parameters obtained by fitting the FM-MCNC model to the MN datasets;
      (8) 'TablesC2C3.RData': calculation results of Tables C.2 and C.3;

"./Experiment2" which collects the intermediately numerical results for Experiment 2:

      (1) subsubfolders "./SIM1", "./SIM2", "./SIM3", and "./SIM4", storing the 'cls.txt', 'estSD.txt', 'modelfit.txt', 'MSEest.txt', 'MSEyc.txt', and 'MSEyfit.txt' text files;
      (2) 'figC4.eps': split violin plots of the MSE for the estimates of entire model parameters and fitted responses by fitting the 3-component FM-MCNR and FM-MCNCR models to the simulated datasets;
      (3) 'figC5.eps': MSE scores for the estimates of model parameters obtained by fitting the 3-component FM-MCNCR model to the simulated data across sample sizes;
      (4) 'TableC4.csv': averages of AIC and BIC scores and frequencies supported by the two criteria for the FM-MCNR and FM-MCNCR models with g=1~4 components to the simulated datasets;
      (5) 'TableC5.csv': Monte Carlo SD and averages of information-based SE of estimates obtained by fitting the 3-component FM-MCNCR model;
      (6) 'TableC5.RData': calculation results of Table C.5.

### Additional Remark ###
      (1) One can directly run each "source(.)" described in the 'master.r' file in the separate R session to obtain the results.
      (2) Numerical results have been stored in "./result/s", and the fitting results for real-data examples have been stored in './data/SaleFit.Rdata and './data/StdWageFit.Rdata'.
      (3) To reproduce the results presented in Figure 1 and Tables 1-2, firstly run the script 'runsale.r' in the subfolder "./code/application/", which produces the results stored in 'SaleFit.RData' file.
      (4) To reproduce Figure 1 and Tables 1-2, just load 'SaleFit.RData' file in the "./data/", and then run 'fig1.r', 'Table1.r', and 'Table2.r' R scripts.
      (5) To reproduce the results presented in Figure 2 and Tables 3-4, firstly run the script 'runwage.r' in the subfolder "./code/application/", which produces the results stored in 'StdWageFit.RData' file.
      (6) To reproduce Figure 2 and Tables 3-4, just load 'StdWageFit.RData' file in the "./data/", and then run 'fig2.r', 'Table3.r', and 'Table4.r' R scripts.
      (7) Because the 'simMCN.r' and 'simMCNR.r' R scripts take a huge amount of time to run, we record these intermediately numerical results so that one can use the R codes 'figC1.r', 'figC2.r', 'figC3ab.r', 'TableC1.r', 'TablesC2C3.r', 'figC4.r', 'figC5.r', 'TableC4.r', and 'TableC5.r' to obtain the final results based on files stored in "./results/Experiment1/SIM1a", "./results/Experiment1/SIM1b", "./results/Experiment1/SIM2a", "./results/Experiment1/SIM2b", "./results/Experiment1/SIM3a", "./results/Experiment1/SIM3b", "./results/Experiment1/SIM4a", "./results/Experiment1/SIM4b", "./results/Experiment2/SIM1", "./results/Experiment2/SIM2", "./results/Experiment2/SIM3", and "./results/Experiment2/SIM4" subfolders.
