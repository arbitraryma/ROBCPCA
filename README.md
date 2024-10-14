# ROBCPCA
This page provides the main code for ROBCPCA paper.
Files: algorithm 1, 2, 3 are the three algorithms for ROBCPCA method

CPCA.R is coded based on the description of 'Multivariate time series clustering based on common principal component analysis' by Hailin Li as a benchmark of our method. The algorithm is called Mc2PCA

We provided the code for the EEG simulation and also the Japanese vowels application. Before running these, you have to run 'algorithm1.R', 'algorithm2.R', 'algorithm3.R', 'cpca.R', 'standarize_cross_covaraince.R', 'true accuracy.R', and 'noise function AR(2) mixing.R'

'algorithm1.R' helps to construct the common projection axes
'algorithm2.R' is the initial grouping 
'algorithm3.R' does the iteration 

'noise function AR(2) mixing.R' added noises from t-distribution with df=3 to some time points as a way to generate potential outliers since t-distribution has heavy tails and is more likely to contain extreme values.

'true accuracy.R' is another measuring tool we use to compare the performance of these methods. It computes the percentage of correctly grouped MTS objects.

