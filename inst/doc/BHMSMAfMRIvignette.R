## ---- out.lines = 10, eval=F--------------------------------------------------
#  library(BHMSMAfMRI)
#  BHMSMAmulti <- BHMSMA(n, grid, data, designmat, k, "multi", truecoef)
#  names(BHMSMAmulti)
#  [1] "GLMCoefStandardized"      "GLMCoefSE"
#  [3] "WaveletCoefficientMatrix" "hyperparam"
#  [5] "hyperparamVar"            "posteriorMixProb"
#  [7] "Waveletcoefposterior"     "GLMcoefposterior"

## ---- eval=T------------------------------------------------------------------
library(BHMSMAfMRI)
fpath <- system.file("extdata", package="BHMSMAfMRI")
untar(paste0(fpath,"/fmridata.tar"), exdir=tempdir())
n <- 3
grid <- 32
ntime <- 9
data <- array(dim=c(n,grid,grid,ntime))
for(subject in 1:n)
{
  directory <- paste0(tempdir(),"/fmridata","/s0",subject,"/")
  a <- readfmridata(directory, format="Analyze", 
  prefix=paste0("s0",subject,"_t"), nimages=ntime, dim.image=c(grid,grid,1))
  data[subject,,,] <- a[,,1,]
}
dim(a)

## ---- eval=T------------------------------------------------------------------
data(fmridata)
names(fmridata)
truecoef <- fmridata$TrueCoeff
designmat <- fmridata$DesignMatrix
dim(truecoef)
dim(designmat)

## ----TrueCoef, fig.cap = "True regression coefficient images for the 3 subjects", fig.width=12, fig.height=4.2----
par(mfrow=c(1,n), cex=1)
for(subject in 1:n)
  image(truecoef[subject,,], main=paste0("Subject ",subject), 
         col=heat.colors(8))

## ---- eval=T------------------------------------------------------------------
glmmap <- glmcoef(n, grid, data, designmat)
names(glmmap)
dim(glmmap$GLMCoefStandardized)
dim(glmmap$GLMCoefSE)

## ----GLMCoef, fig.cap = "Standardized regression coefficient estimates images for the second regressor for all subjects", fig.width=12, fig.height=4.2----
k <- 2
par(mfrow=c(1,n), cex=1)
for(subject in 1:n)
  image(abs(glmmap$GLMCoefStandardized[subject,,,k]), col=heat.colors(8),
         zlim=c(0,6), main=paste0("Subject ",subject))

## ---- eval=T------------------------------------------------------------------
wavecoefglmmap <- waveletcoef(n, grid, glmmap$GLMCoefStandardized[,,,k], 
            wave.family="DaubLeAsymm", filter.number=6, bc="periodic")
names(wavecoefglmmap)
dim(wavecoefglmmap$WaveletCoefficientMatrix)

## ---- eval=T--------------------------------------------------------------------------------------
options(width = 100)
hyperest <- hyperparamest(n, grid, wavecoefglmmap$WaveletCoefficientMatrix, 
                           analysis = "multi")
names(hyperest)
round(hyperest$hyperparam,3)
signif(hyperest$hyperparamVar,4)

## ---- eval=T--------------------------------------------------------------------------------------
a.kl <- hyperest$hyperparam[1] * 2^(-hyperest$hyperparam[2] * (0:4))
b.kl <- hyperest$hyperparam[3] * 2^(-hyperest$hyperparam[4] * (0:4))
c.kl <- hyperest$hyperparam[5] * 2^(-hyperest$hyperparam[6] * (0:4))
round(a.kl,3)

## ---- eval=T--------------------------------------------------------------------------------------
pkljbar <- postmixprob(n, grid, wavecoefglmmap$WaveletCoefficientMatrix, 
                        hyperest$hyperparam, analysis = "multi")
names(pkljbar)
dim(pkljbar$pkljbar)
round(pkljbar$pkljbar[1,1:10],4)

## ---- eval=T--------------------------------------------------------------------------------------
postwavecoefglmmap <- postwaveletcoef(n, grid, 
                         wavecoefglmmap$WaveletCoefficientMatrix, 
                 hyperest$hyperparam, pkljbar$pkljbar, analysis = "multi")
names(postwavecoefglmmap)
dim(postwavecoefglmmap$PostMeanWaveletCoef)
dim(postwavecoefglmmap$PostMedianWaveletCoeff)

## ---- eval=T--------------------------------------------------------------------------------------
postglmmap <- postglmcoef(n, grid, glmmap$GLMCoefStandardized[,,,k], 
         postwavecoefglmmap$PostMeanWaveletCoef, wave.family="DaubLeAsymm", 
           filter.number=6, bc="periodic")
str(postglmmap,vec.len = 3, digits.d = 2)

## ----PostCoef, fig.cap = "Posterior standardized regression coefficient images for the 3 subjects obtained by BHMSMA", fig.width=12, fig.height=4.2----
par(mfrow=c(1,n), cex=1)
for(subject in 1:n)
  image(abs(postglmmap$GLMcoefposterior[subject,,]), col=heat.colors(8),
         zlim=c(0,6), main=paste0("Subject ",subject))

## ---- eval=T--------------------------------------------------------------------------------------
MSE <- c()
for (i in 1:n) 
 MSE[i] <- sum((as.vector(truecoef[i,,]/glmmap$GLMCoefSE[i,,,2])  
                - as.vector(postglmmap$GLMcoefposterior[i,,]))^2)
round(MSE,3)

## ---- eval=T, fig.width=12, fig.height=4.2--------------------------------------------------------
Postsamp <- postsamples( nsample=50, n, grid, glmmap$GLMCoefStandardized[,,,k],
               wavecoefglmmap$WaveletCoefficientMatrix,  hyperest$hyperparam, 
               pkljbar$pkljbar, "multi", seed=123)
names(Postsamp)
dim(Postsamp$samples)
dim(Postsamp$postdiscovery)

## ----PostDiscovery, eval=T, fig.cap = "Posterior discovery images for the 3 subjects", fig.width=12, fig.height=4.2----
par(mfrow=c(1,n))
for(subject in 1:n)
  image(Postsamp$postdiscovery[subject,,], col=heat.colors(8),
   main=paste0("Subject ",subject))

## ---- eval=T--------------------------------------------------------------------------------------
postsd <- array(dim=c(n,grid,grid))
for(subject in 1:n)
 postsd[subject,,] <- apply(Postsamp$samples[subject,,,], 1:2, sd)
                             round(postsd[1,1:5,1:5],3)

## ---- eval=T--------------------------------------------------------------------------------------
postgroup <- postgroupglmcoef( n, grid, glmmap$GLMCoefStandardized[,,,k], 
                               postwavecoefglmmap$PostMeanWaveletCoef)
names(postgroup)
dim(postgroup$groupcoef)

## ----PostGroupCoef, fig.cap = "Posterior group regression coefficient image",  fig.width=3, fig.height=3, fig.align="center"----
par(mfrow=c(1,1),cex=0.5)
image(abs(postgroup$groupcoef),col=heat.colors(8),zlim=c(0,6))

