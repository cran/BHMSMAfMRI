\name{fmridata}
\docType{data}
\alias{fmridata}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{A simulated fMRI data for 3 subjects
}

\description{A simulated fMRI data containing true regression coefficients images for three subjects and design matrix
}

\usage{
data(fmridata)
}

\format{A list containing the following. 
 \describe{
  \item{grid }{ =32. The image dimension is 32 by 32.}
  \item{nsubject }{ =3.}
  \item{TrueCoeff }{ An array of dimension (3,32,32), containing the true regression coefficients for the 3 subjects.}
  \item{DesignMatrix }{ A matrix with 9 columns and 2 rows. The first column is a column of ones.}
 }
}

\details{This dataset contains only the true coefficients. The noisy fMRI data, which are generated by adding Gaussian random noise to these true coefficients, are included in the extdata directory within the package directory. The function \emph{readfmridata} can be used to read those data files. The true coefficients and the noisy data both are generated using the R package neuRosim. The following specifications were used to generate the data: totaltime=18, onsets=seq(1,18,by=8), durations=1, TR=2, effectsize=1, hrf="double-gamma", regions=3, radius=c(1,1,1), form="sphere", fading=1, SNR=1.5, noise="white". The centers of the activation regions were chosen manually. For information regarding the specifications, see neuRosim help.}

