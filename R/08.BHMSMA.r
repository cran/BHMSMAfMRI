
BHMSMA = function( nsubject, grid, Data, DesignMatrix, TrueCoeff=NULL, analysis, wave.family="DaubLeAsymm", filter.number=6, bc="periodic")
{

 if(analysis == "multi")
 {
   glmDat = cmpfun(glmcoeff)( nsubject, grid, Data, DesignMatrix )

   wavecoeff.glmDat = waveletcoeff(nsubject, grid, glmDat$GLMCoeffStandardized, wave.family, filter.number, bc)

   hyparEst = cmpfun(hyperparamest)( nsubject, grid, wavecoeff.glmDat$WaveletCoefficientMatrix, analysis="multi" )

   pklj.bar = cmpfun(pikljbar)(nsubject, grid, wavecoeff.glmDat$WaveletCoefficientMatrix, hyparEst$hyperparam, analysis="multi")

   postwavelet = cmpfun(postwaveletcoeff)( nsubject, grid, wavecoeff.glmDat$WaveletCoefficientMatrix, hyparEst$hyperparam, pklj.bar=pklj.bar$pklj.bar, analysis="multi" )

   postglm = cmpfun(postglmcoeff)(nsubject, grid, glmDat$GLMCoeffStandardized, postwavelet$PostMeanWaveletCoeff, wave.family, filter.number, bc)
 }
 
 if(analysis == "single")
 {
   glmDat = cmpfun(glmcoeff)( nsubject, grid, Data, DesignMatrix )

   wavecoeff.glmDat = waveletcoeff(nsubject, grid, glmDat$GLMCoeffStandardized, wave.family, filter.number, bc)

   hyparEst = cmpfun(hyperparamest)( nsubject, grid, wavecoeff.glmDat$WaveletCoefficientMatrix, analysis="single" )

   pklj.bar = cmpfun(pikljbar)(nsubject, grid, wavecoeff.glmDat$WaveletCoefficientMatrix, hyparEst$hyperparam, analysis="single")

   postwavelet = cmpfun(postwaveletcoeff)( nsubject, grid, wavecoeff.glmDat$WaveletCoefficientMatrix, hyparEst$hyperparam, pklj.bar=pklj.bar$pklj.bar, analysis="single" )

   postglm = cmpfun(postglmcoeff)(nsubject, grid, glmDat$GLMCoeffStandardized, postwavelet$PostMeanWaveletCoeff, wave.family, filter.number, bc)
 }
 
 if(! is.null(TrueCoeff))
 {
  MSE = c()
  for(i in 1:nsubject)
  {
   MSE[i] = sum ( ( as.vector(TrueCoeff[i,,]/glmDat$GLMEstimatedSE[i,,]) - as.vector(postglm$GLMcoeffposterior[i,,]) )^2 )
  }
 }
 
 if(! is.null(TrueCoeff))
  return(list( GLMCoeffStandardized = glmDat$GLMCoeffStandardized, GLMEstimatedSE = glmDat$GLMEstimatedSE, WaveletCoefficientMatrix = wavecoeff.glmDat$WaveletCoefficientMatrix, hyperparam = hyparEst$hyperparam, hyperparamVar = hyparEst$VarMLE, pklj.bar=pklj.bar$pklj.bar,  Waveletcoeffposterior = postwavelet$PostMeanWaveletCoeff, GLMcoeffposterior = postglm$GLMcoeffposterior, MSE=MSE))

 if(is.null(TrueCoeff))
  return(list( GLMCoeffStandardized = glmDat$GLMCoeffStandardized, GLMEstimatedSE = glmDat$GLMEstimatedSE, WaveletCoefficientMatrix = wavecoeff.glmDat$WaveletCoefficientMatrix, hyperparam = hyparEst$hyperparam, hyperparamVar = hyparEst$VarMLE, pklj.bar=pklj.bar$pklj.bar,  Waveletcoeffposterior = postwavelet$PostMeanWaveletCoeff, GLMcoeffposterior = postglm$GLMcoeffposterior))
 
}



