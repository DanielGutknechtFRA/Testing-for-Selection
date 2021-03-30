# Testing-for-Selection
This repository contains R code to perform nonparametric (distribution free) tests for sample selection. The code replicates (parts of) the empirical illustration and the simulation results from "Testing for Quantile Sample Selection" by V. Corradi and Daniel Gutknecht (permanent link: https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3805043). Note that all page references below are w.r.t. the current version of the paper from this link. Please refer to the paper for further details.

There are two subfolders contained in this repository:

1. Empirical Illustration.
2. Simulations.

# 1. Empirical Illustration

This folder contains all R codes required to reproduce (parts of) the empirical illustration in Tables 3 (p.30) and 4 (p.31). The codes can be found in the following subfolders:

# 1.1 Mean:
 - CM-SimsFINAL.R - Master R file of the simulations for the conditional mean. 
 - CMTest1FINAL.R - R file containing routines for the first test. 
 - CMTest2FINAL.R - R file containing routines for the second test.
# 1.2 Quantile:
 - CQ-SimsFINAL-Test-1 - R file containing routines for the first (conditional quantile) test. 
 - CQ-SimsFINAL-Test-2-BW-Select - R file containing routines for the second (conditional quantile) test (with data-driven bandwidth selection). 
 - CQ-SimsFINAL-Test-2-NO-BW-Select - R file containing routines for the second (conditional quantile) test (without data-driven bandwidth selection). 


# 2. Simulations

This folder contains all R codes required to reproduce the simulation results in Tables 1 (p.23), 2 (p.24), and 5 (p.40). The codes can be found in the following subfolders:

# 2.1 Mean:
 - CM-SimsFINAL.R - Master R file of the simulations for the conditional mean. 
 - CMTest1FINAL.R - R file containing routines for the first test. 
 - CMTest2FINAL.R - R file containing routines for the second test.
# 2.2 Quantile:
 - CQ-SimsFINAL-Test-1 - R file containing routines for the first (conditional quantile) test. 
 - CQ-SimsFINAL-Test-2-BW-Select - R file containing routines for the second (conditional quantile) test (with data-driven bandwidth selection). 
 - CQ-SimsFINAL-Test-2-NO-BW-Select - R file containing routines for the second (conditional quantile) test (without data-driven bandwidth selection). 
