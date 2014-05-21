
postwaveletcoeff <- function( nsubject, grid, WaveletCoefficientMatrix, hyperparam, pklj.bar, analysis )
{

 # ......................... Posterior Mean of Wavelet Coefficients ..............................
 
 if(analysis == "multi")
 {
   C4 <- hyperparam[5]
   C5 <- hyperparam[6]
   
   PostMeanWaveletCoeff <- matrix( nrow=nsubject, ncol=grid^2-1 )
   PostMedianWaveletCoeff <- matrix( nrow=nsubject, ncol=grid^2-1 )
   for(subject in 1:nsubject)
   {
    
    PostMean <- PostMedian <- c()
    for(l in c(0:(log2(grid)-1)))
    for(j in c(1:(2^(2*l)*3)))
    {   
      if(l==0) d <- WaveletCoefficientMatrix[subject,j] else d <-  WaveletCoefficientMatrix[subject, sum(2^(2*(0:(l-1)))*3) + j] 
      if(l==0) p <- pklj.bar[subject,j] else p <- pklj.bar[subject,sum(2^(2*(0:(l-1)))*3) + j]   
      c.l <- C4 * 2^(-C5 * l)
      mean <- p * c.l/(1+c.l) * d
      PostMean <- c(PostMean,mean)
      if(p<=0.5) median <- 0 else median <- c.l/(1+c.l) * d - sign(d) * (c.l/(1+c.l))^(.5) * qnorm(0.5/p)
      PostMedian <- c(PostMedian,mean)
    } 
    PostMeanWaveletCoeff[subject,] <- PostMean
    PostMedianWaveletCoeff[subject,] <- PostMedian
   }

   return(list(PostMeanWaveletCoeff=PostMeanWaveletCoeff, PostMedianWaveletCoeff=PostMedianWaveletCoeff))
   
   }
 
 if(analysis == "single")
 {
   
   PostMeanWaveletCoeff <- matrix( nrow=nsubject, ncol=grid^2-1 )
   PostMedianWaveletCoeff <- matrix( nrow=nsubject, ncol=grid^2-1 )
   for(subject in 1:nsubject)
   {
     C4 <- hyperparam[subject,5]
     C5 <- hyperparam[subject,6]
     
     PostMean <- PostMedian <- c()
     for(l in c(0:(log2(grid)-1)))
     for(j in c(1:(2^(2*l)*3)))
     {   
       if(l==0) d <- WaveletCoefficientMatrix[subject,j] else d <-  WaveletCoefficientMatrix[subject, sum(2^(2*(0:(l-1)))*3) + j] 
       if(l==0) p <- pklj.bar[subject,j] else p <- pklj.bar[subject,sum(2^(2*(0:(l-1)))*3) + j]   
       c.l <- C4 * 2^(-C5 * l)
       mean <- p * c.l/(1+c.l) * d
       PostMean <- c(PostMean,mean)
       if(p<=0.5) median <- 0 else median <- c.l/(1+c.l) * d - sign(d) * (c.l/(1+c.l))^(.5) * qnorm(0.5/p)
       PostMedian <- c(PostMedian,mean)
     }
    
     PostMeanWaveletCoeff[subject,] <- PostMean
     PostMedianWaveletCoeff[subject,] <- PostMedian
   }
   
   return(list(PostMeanWaveletCoeff=PostMeanWaveletCoeff, PostMedianWaveletCoeff=PostMedianWaveletCoeff))
 }
 
}
