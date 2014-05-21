
waveletcoeff <- function( nsubject, grid, GLMCoeffStandardized, wave.family="DaubLeAsymm", filter.number=6, bc="periodic" )
{

 # ....................... Wavelet Transform ................

 WaveletCoefficientMatrix <- matrix( nrow=nsubject, ncol=grid^2-1 )

 for(subject in 1:nsubject)
 {
  dwt <- imwd(GLMCoeffStandardized[subject,,], type="wavelet", family=wave.family, filter.number=filter.number, bc=bc, RetFather=TRUE, verbose=FALSE)
 
  WaveletCoefficientMatrix[subject,] <- c(dwt$w0L1,dwt$w0L2,dwt$w0L3,dwt$w1L1,dwt$w1L2,dwt$w1L3,dwt$w2L1,dwt$w2L2,dwt$w2L3,dwt$w3L1,dwt$w3L2,dwt$w3L3,dwt$w4L1,dwt$w4L2,dwt$w4L3,dwt$w5L1,dwt$w5L2,dwt$w5L3,dwt$w6L1,dwt$w6L2,dwt$w6L3,dwt$w7L1,dwt$w7L2,dwt$w7L3,dwt$w8L1,dwt$w8L2,dwt$w8L3)     #.....Note: Using up to w8. So, applicable only up to 2^9 by 2^9 data 
 }
 return(list(WaveletCoefficientMatrix=WaveletCoefficientMatrix))
 
}