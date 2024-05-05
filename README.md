# Testing-for-Selection
This repository contains R code to perform nonparametric (distribution free) tests for sample selection. The code replicates the empirical analysis and the simulation from the paper "Testing for Quantile Sample Selection" by V. Corradi and Daniel Gutknecht (permanent link to working paper: https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3805043; published version: https://doi.org/10.1093/ectj/utac027). Note that all page references below are w.r.t. the current version of the paper from this link. Please refer to the paper for further details.

There are two subfolders contained in this repository:

1. Empirical Illustration.
2. Simulations.

# 1. Empirical Illustration

This folder contains all R codes for the empirical illustration (the files use data-driven bandwiths and do not reproduce the tables of the published paper, see below for the latter):

 - LSupply95-97_20-04-22.R - R file of that reads in the data, estimates all required quantities, and carries out tests for sample selection for the nonparametric conditional   mean (not in the paper) and quantile on the subsample 1995-1997. 
 - LSupply98-00_20-04-22.R - R file of that reads in the data, estimates all required quantities, and carries out tests for sample selection for the nonparametric conditional   mean (not in the paper) and quantile on the subsample 1998-2000. 
 - CM-A_Tests_20-04-22.R - R file containing relevant functions for the sample selection tests for the conditional mean function. 
 - CQ-A_Tests_20-04-22.R - R file containing relevant functions for the sample selection tests for the conditional quantile function(s). 
 - data_1.csv - CSV file containing wage data in  .csv format from the UK Expenditure Survey. Source of data: Arellano, M. and S. Bonhomme (2017). Quantile selection models with an application to understanding changes in wage inequality. Econometrica 85 (1), 1-28. (http://dx.doi.org/10.3982/ECTA14030)
 - data_1.csv - CSV file containing wage data in  .csv format from the UK Expenditure Survey. Source of data: Arellano, M. and S. Bonhomme (2017). Quantile selection models with an application to understanding changes in wage inequality. Econometrica 85 (1), 1-28. (http://dx.doi.org/10.3982/ECTA14030)

# published version (replicates all Tables from the empirical illustration of the published paper):
 - LSupply95-97_08-09-22.R - R file reproduces the results from the empirical illustration in Tables 1 and 2. 
 - LSupply98-00_08-09-22.R - R file reproduces the results from the empirical illustration in Tables 1 and 2.
 - BWs-1995-1997.csv: CSV file containing the cross-validated bandwidths of different nonparametric objects used to construct the test statistics in Tables 1 and 2 (1995-1997).
 - BWs-1998-2000.csv - CSV file containing the cross-validated bandwidths of different nonparametric objects used to construct the test statistics in Tables 1 and 2 (1998-2000).


# 2. Simulations

This folder contains all R codes for the Monte Carlo simulations (the files do not reproduce the tables of the published paper, see below for the latter). The codes can be found in the following subfolders:

# 2.1 Mean:
 - CM-SimsFINAL.R - Master R file of the simulations for the conditional mean. 
 - CMTest1FINAL.R - R file containing routines for the first test. 
 - CMTest2FINAL.R - R file containing routines for the second test.
# 2.2 Quantile:
 - CQ-SimsFINAL-TEST-1_20-04-22 - R file containing routines for the first (conditional quantile) test. 
 - CQ-SimsFINAL-TEST-2-BW-Select_20-04-22 - R file containing routines for the second (conditional quantile) test (with data-driven bandwidth selection). 
 - CQ-SimsFINAL-TEST-2-NO-BW-Select_20-04-22 - R file containing routines for the second (conditional quantile) test (without data-driven bandwidth selection).

# published version (replicates all Tables from the Monte Carlo simulations of the published paper):
 - CG-Sims-Table-1-FINAL-08-09-22- R file containing routines for the first (conditional) quantile test that reproduces TABLE 1 of supplementary material.
 - CG-Sims-Table-2&3-FINAL-08-09-22 - the R file reproduces the results of TABLES 2 and 3 of the supplementary material.
