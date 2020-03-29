# perm_test_ci
This function performs a permutation test for the difference of two group means, 
see for example:

+ Ernst, M. D. (2004). Permutation methods: a basis for exact inference.                            
*Statistical Science*, 19(4), 676-685.
                                                              
It calculates an exact two-sided p-value for the Null hypothesis of no difference between groups.
It creates a histogram of the simulated permutation distribution of mean differences.
Confidence intervals for the mean difference are approximated via the stochastic search 
method described in:

+ Garthwaite, P. H. (1996). Confidence intervals from randomization tests.                          
*Biometrics*, 1387-1393.

Diagnostic plots for the search process are provided.
For feedback, suggestions and errors, please contact stefan-gehrig[at]t-online.de

Support for one-tailed tests and other test statistics is intended for an updated version.