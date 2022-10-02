
# BHMSMAfMRI v2.1 

* Changed package URL.

* Added welcome message.

* Added NEWS file properly.


# BHMSMAfMRI v2.0

* Major computational upgrade: The core computational parts of most of the functions were written in C++ using Rcpp and RcppArmadillo libraries. This led to a major boost in computational speed.

* Function name changes: read.fmridata to readfmridata, glmcoeff to glmcoef, waveletcoeff to waveletcoef, pikljbar to postmixprob, postwaveletcoeff to postwaveletcoef, postglmcoeff to postglmcoef, postgroupcoeff to postgroupglmcoef.

* Function added: One new R function, substituteWaveletCoef, was added for convenience.

* Argument changes: The argument 'nsubject' (number of subjects) previously used in multiple functions was renamed as 'n'. The main function BHMSMA received a new argument(k) denoting which regressor to consider for analysis.

* Function scope changes: The function glmcoef was broadened in scope to include any number of regressors.

* Bug fixes: A bug related to the computation of posterior median of the wavelet coefficients in the function postwaveletcoeff was fixed. A bug related to return of the hyperparameter variance estimates in the function BHMSMA was fixed.

* Others: The package compiler was removed from 'depends' as automatic byte-compilation of packages on installation became default in latest R versions.
		  The vignette is updated.


# BHMSMAfMRI v1.3

* Fixed some issues in the function manual files.


# BHMSMAfMRI v1.2

* Changes in main author (corresponding author) email and affiliation.

* Changes in package dependency: Dependency on the packages "fmri" and "ANALYZEfmri" were removed, and in that place, dependency on the package "oro.nifti" was added. The only function updated according with the above changes is read.fmridata.


# BHMSMAfMRI v1.1

* Changes in main author (corresponding author) email and affiliation. 

* Changes in the function BHMSMA(): The function cmpfun from the package Compiler is incorporated within the function BHMSMA to increase computaional speed. Package Compiler is included in the list of dependencies in the DESCRIPTION file.
