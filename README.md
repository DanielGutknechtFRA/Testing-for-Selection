# Testing-for-Selection
This repository contains R code to perform nonparametric (distribution free) tests for sample selection. The code replicates (parts of) the empirical illustration and the simulation results from "Testing for Quantile Sample Selection" by V. Corradi and Daniel Gutknecht (permanent link: https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3805043). Note that all page references below are w.r.t. the current version of the paper from this link. Please refer to the paper for further details.

There are two subfolders contained in this repository:

1. Empirical Illustration.
2. Simulations.

# 1. Empirical Illustration

This folder contains all R codes required to reproduce (parts of) the empirical illustration in Tables 3 (p.30) and 4 (p.31):

 - LSupply95-97_20-04-22.R - Master R file of which reads in the data, estimates all required quantities, and carries out tests for sample selection for the nonparametric conditional   mean (not in the paper) and quantile on the subsample 1995-1997. 
 - LSupply98-00_20-04-22.R - Master R file of which reads in the data, estimates all required quantities, and carries out tests for sample selection for the nonparametric conditional   mean (not in the paper) and quantile on the subsample 1998-2000. 
 - CM-A_Tests_20-04-22.R - R file containing routines for the sample selection tests for the conditional mean function. 
 - CQ-A_Tests_20-04-22.R - R file containing routines for the sample selection tests for the conditional quantile function(s). 
 - data_1.csv - CSV file containing wage data in  .csv format from the UK Expenditure Survey. Source of data: Arellano, M. and S. Bonhomme (2017). Quantile selection models with an application to understanding changes in wage inequality. Econometrica 85 (1), 1-28. (http://dx.doi.org/10.3982/ECTA14030)
 - data_1.csv - CSV file containing wage data in  .csv format from the UK Expenditure Survey. Source of data: Arellano, M. and S. Bonhomme (2017). Quantile selection models with an application to understanding changes in wage inequality. Econometrica 85 (1), 1-28. (http://dx.doi.org/10.3982/ECTA14030)

# 2. Simulations

This folder contains all R codes required to reproduce the simulation results in Tables 1 (p.23), 2 (p.24), and 5 (p.40). The codes can be found in the following subfolders:

# 2.1 Mean:
 - CM-SimsFINAL.R - Master R file of the simulations for the conditional mean. 
 - CMTest1FINAL.R - R file containing routines for the first test. 
 - CMTest2FINAL.R - R file containing routines for the second test.
# 2.2 Quantile:
 - CQ-SimsFINAL-TEST-1_20-04-22 - R file containing routines for the first (conditional quantile) test. 
 - CQ-SimsFINAL-TEST-2-BW-Select_20-04-22 - R file containing routines for the second (conditional quantile) test (with data-driven bandwidth selection). 
 - CQ-SimsFINAL-TEST-2-NO-BW-Select_20-04-22 - R file containing routines for the second (conditional quantile) test (without data-driven bandwidth selection). 
