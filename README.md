# Grenander confidence intervals

This repository gives an R script for constructing confidence intervals for a decreasing density, based on the Grenander estimator of a decreasing density.

The method is discussed in "Nonparametric confidence intervals for monotone functions", Annals of Statistics Volume 43, Number 5 (2015), 2019-2054, of Piet Groeneboom and Geurt Jongbloed. The paper is given here as Isotonic.pdf. Confidence intervals of the present type were first introduced by Banerjee and Wellner for current status data: "Likelihood ratio tests for monotone functions", Ann. Statist. 29 (2001), 1699–1731, and "Confidence intervals for current status data", Scand. J. Statist. 32 (2005), 405–424.

The difficulty in constructing the intervals is the condition that the integral over the density is equal to 1, also if one fixes the value of the density at a paqrticular point. This difficulty is solved by Lemma 3.2 in "Nonparametric confidence intervals for monotone functions", where a equation is given for determining the needed Lagrange parameter mu. This equation is changed to an equation for the second parameter lambda. The parameters have the relation mu=(lambda-1)/a, where $a$ is the value of the density at the point of evaluation, see (3.8) of the manuscript. The value of lambda is found in the program by using golden section search.

The R script GJ2015_CI.R gives an example of a sample of size 100 from the standard truncated exponential distribution on the interval [0,2]. One only needs to have running the R package Rcpp to be able to run the R script.
