\name{BHMSMAfMRI-package}

\alias{BHMSMAfMRI}

\docType{package}

\title{Bayesian Hierarchical Multi-Subject Multiscale Analysis (BHMSMA) of Functional MRI Data}

\description{The \pkg{BHMSMAfMRI} package performs BHMSMA (Sanyal & Ferreira, 2012) of fMRI data, or other multiscale data, using wavelet based prior that borrows strength across subjects and provides posterior smoothed images of the effect sizes and samples from the posterior distribution. The package currently considers analysis of 2D slices/grids only.
}

\details{
\tabular{ll}{
Package: \tab BHMSMAfMRI\cr
Type: \tab Package\cr
Version: \tab 2.0\cr
Date: \tab 2022-09-06\cr
License: \tab GPL (>= 2)\cr
}

Import fMRI data using: \cr
\code{\link{readfmridata}}

The main analysis function, which provides subject-specific posterior estimates, is: \cr
\code{\link{BHMSMA}}

The main function sucessively calls the following functions: \cr
\code{\link{glmcoef}} (get regression coefficients) \cr
\code{\link{waveletcoef}} (get wavelet coefficients) \cr
\code{\link{hyperparamest}} (estimate model hyperparameters) \cr
\code{\link{postmixprob}} (estimate posterior mixture probabilities of wavelet coefficients) \cr
\code{\link{postwaveletcoef}} (compute posterior estimates of wavelet coefficients) \cr
\code{\link{postglmcoef}} (compute posterior estimates of regression coefficients)

For posterior group estimates of regression coefficients use: \cr
\code{\link{postgroupglmcoef}}

For posterior uncertainty estimates use: \cr
\code{\link{postsamples}}

Internal sample data: \cr
\code{\link{fmridata}}

Miscellaneous: \cr
\code{\link{substituteWaveletCoef}}
}

\author{Nilotpal Sanyal <nilotpal.sanyal@gmail.com>, Marco Ferreira <marf@vt.edu>

Maintainer: Nilotpal Sanyal <nilotpal.sanyal@gmail.com>
}

\references{
Sanyal, Nilotpal, and Ferreira, Marco A.R. (2012). Bayesian hierarchical multi-subject multiscale analysis of functional MRI data. Neuroimage, 63, 3, 1519-1531.
}

%~~ Optionally other standard keywords, one per line, from file KEYWORDS in ~~
%~~ the R documentation directory ~~
%\keyword{ package }
%\seealso{
%~~ Optional links to other man pages, e.g. ~~
%~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
%}
