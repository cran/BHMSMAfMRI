\name{glmcoeff}

\alias{glmcoeff}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{Fit GLM to the data time-series and obtain GLM coefficients along with standard error estimates
}

\description{Fits General Linear Model to the time-series corresponding to each voxel in the data and returns the standardized GLM coefficients and their standard error estimates.
}

\usage{
glmcoeff(nsubject, grid, Data, DesignMatrix)
}

\arguments{
  \item{nsubject}{Number of subjects included in the analysis.}
  \item{grid}{The number of voxels in one row (or, one column) of the brain slice of interest. Must be a power of 2. The total number of voxels is grid^2. The maximum grid value for this package is 512.}
  \item{Data}{The data in form of an array with dimension (nsubject,grid,grid,ntime), where ntime is the size of the time series for each voxel.}
  \item{DesignMatrix}{The design matrix used to generate the data.}
}

\value{A list containing the following.
  \item{GLMCoeffStandardized }{An array of dimension (nsubject, grid, grid), containing for each subject the standardized GLM coefficients obtained by fitting GLM to the time-series corresponding to the voxels.}
  \item{GLMEstimatedSE }{An array of dimension (nsubject, grid, grid), containing for each subject the estimated standard errors of the GLM coefficients.}
}

\references{
Friston, K.J., Holmes, A.P., Worsley, K.J., Poline, J., Frith, C.D., Frackowiak, R.S.J., 1994. Statistical parametric maps in functional imaging: a general linear approach. Hum. Brain Mapp. 2 (4), 189-210.
}

\author{Nilotpal Sanyal <nisanyal@ucsd.edu>, Marco Ferreira <marf@vt.edu>
}

%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}

\examples{
nsubject <- 3
grid <- 8
Data <- array(dim=c(3,8,8,10),rnorm(3*8*8*10))
DesignMatrix <- cbind( c(rep(c(1,0),5)), rep(1,10) )
glm.fit <- glmcoeff(nsubject, grid, Data, DesignMatrix)
dim(glm.fit$GLMCoeffStandardized)
#[1] 3 8 8
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%%\keyword{ ~kwd1 }
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
%}