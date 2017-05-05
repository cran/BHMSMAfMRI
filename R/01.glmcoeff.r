
glmcoeff = function( nsubject, grid, Data, DesignMatrix )
{

 GLMCoefficientsNonstandardized = GLMCoefficientsStandardized = GLMEstimatedVariance = array( 0, dim = c( nsubject, grid, grid ) ) 

 for(subject in 1:nsubject)
 {
 
  for(i in 1:grid)   
  for(j in 1:grid)
  { 
   y = Data[subject,i,j,]    
   if (! sum(y)==0)
    {
      glm.fit = lm(y ~ 0 + DesignMatrix)

      GLMCoefficientsNonstandardized[subject,i,j] = glm.fit $coefficients[-1]

      mse = anova(glm.fit)[,3][2]

      GLMEstimatedVariance[subject,i,j] = mse * solve( t(DesignMatrix) %*% DesignMatrix )[2,2]

      GLMCoefficientsStandardized[subject,i,j] = GLMCoefficientsNonstandardized[subject,i,j] / sqrt(GLMEstimatedVariance[subject,i,j])
    }
  }
 
 }
 return( list( GLMCoeffStandardized = GLMCoefficientsStandardized, GLMEstimatedSE = sqrt(GLMEstimatedVariance) ) )
 
}