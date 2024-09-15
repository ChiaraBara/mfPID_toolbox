mfPID - Matlab toolbox for the model-free estimation of Partial Information Decomposition (PID)

References:
- Williams PL, Beer RD, "Nonnegative decomposition of multivariate information". arXiv preprint arXiv:1004.2515, 2010.
- Barà C et al., "artial information decomposition for mixed discrete and continuous random variables". 

Demonstration scripts
_____________________
Simulation - PID measures and significance on Boolean logic AND-gate 
Application - PID measures evaluated considering ventilatory phases as target and cardiac and vascular dynamics as sources for a representative subject (load data_application.mat)

Main computational functions
_____________________________
mfPID_H - Calculate entropy for discrete multidimensional variable by using the frequentistic approach
mfPID_2sources_th - Compute theoretical PID measures for two discrete sources and one discrete target from joint probabilities
mfPID_redundancy_lattice_2sources - Compute partial information quantities of the four atoms of the redundancy lattice
mfPID_2sources_dicrete - Compute PID terms for two discrete sources and one discrete target from observation matrix
mfPID_2sources_mixed - Compute PID terms for two continuous sources and one discrete target from observation matrix
mfPID_2sources_mixed_mex - Compute PID terms for two continuous sources and one discrete target from observation matrix (this makes use of closed mex functions)
mfPID_3sources_th - Compute theoretical PID measures for three discrete sources and one discrete target from joint probabilities
mfPID_redundancy_lattice_3sources - Compute partial information quantities of the eighteen atoms of the redundancy lattice
mfPID_3sources_dicrete - Compute PID terms for three discrete sources and one discrete target from observation matrix
mfPID_3sources_mixed - Compute PID terms for three continuous sources and one discrete target from observation matrix
mfPID_3sources_mixed_mex - Compute PID terms for three continuous sources and one discrete target from observation matrix (this makes use of closed mex functions)
mfPID_quantization - Discretize uniformly a series by using the binning approach (fixing number of bins)
mfPID_nuquantization - Discretize non uniformly a series (fixing the size of the quantization levels)
mfPID_B - Generate the observation matrix indicating number of samples in the past for each series
mfPID_B_lags - Generate the observation matrix indicating number lags of the accounted samples (both past and future)
mfPID_ObsMat_V, mfPID_ObsMat_lags, mfPID_SetLag - Form tmp observation matrix 

Other functions
________________
surrshuf - Generate surrogate data shuffling time series samples randomly
nn_search, nn_prepare, range_search - mex functions for nearest neighbors searching and counting
