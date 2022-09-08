############# BHMSMAfMRI package R functions #############

# Function to compute GLM (general linear model) coefficients of the regressor variables for all voxels of a single brain slice for each subject 
#
glmcoef = function( n, grid, data, designmat )
{
    glm_coef_standardized = glm_coef_se = array( 0, dim = c( n, grid, grid, ncol(designmat) ) ) 

    for(i in 1:n)
    {
        out = glmcoef_sub(grid, data[i,,,], designmat)
        glm_coef_standardized[i,,,] = out$GLM_coef_st
        glm_coef_se[i,,,] = out$GLM_coef_se
    }
    return( list( GLMCoefStandardized = glm_coef_standardized, GLMCoefSE = glm_coef_se ) )
}

# Function to compute wavelet transform (coefficients) of a 2D GLM coefficient map (e.g., corresponding to a single brain slice) of the regressor of interest for each subject
#
waveletcoef = function( n, grid, glmcoefstd, wave.family = "DaubLeAsymm", filter.number = 6, bc = "periodic" )
{
   # ....................... Wavelet Transform ................

    WaveletCoefficientMatrix = matrix( nrow = n, ncol = grid^2-1 )

    for(i in 1:n)
    {
        dwt = imwd(glmcoefstd[i,,], type = "wavelet", family = wave.family, filter.number = filter.number, bc = bc, RetFather = TRUE, verbose = FALSE)

        WaveletCoefficientMatrix[i,] = c(dwt$w0L1, dwt$w0L2, dwt$w0L3, dwt$w1L1, dwt$w1L2, dwt$w1L3, dwt$w2L1, dwt$w2L2, dwt$w2L3, dwt$w3L1, dwt$w3L2, dwt$w3L3, dwt$w4L1, dwt$w4L2, dwt$w4L3, dwt$w5L1, dwt$w5L2, dwt$w5L3, dwt$w6L1, dwt$w6L2, dwt$w6L3, dwt$w7L1, dwt$w7L2, dwt$w7L3, dwt$w8L1, dwt$w8L2, dwt$w8L3)     #.....Note: Using up to w8. So, applicable only up to 2^9 by 2^9 data 
    }
    return( list(WaveletCoefficientMatrix = WaveletCoefficientMatrix) )
 
}

# Function to compute estimates of the hyperparameters that appear in the BHMSMA model based on multi-subject or single subject analyses (see References)
#
hyperparamest = function( n, grid, waveletcoefmat, analysis="multi" )
{
    #....................... Estimating C0, C1, C2, C3, C4 and C5..................

    if(analysis == "multi")
    {
        C0 = C1 = C2 = C3 = C4 = C5 = 1
        
        for(j in 1:50)
        {
            max_likelihood = nlminb(start = C5, objective=minus_ll, lower=0, upper=Inf, C0=C0, C1=C1, C2=C2, C3=C3, C4=C4, subs=1:n, grid=grid, waveletcoefmat=waveletcoefmat)
            C5 = max_likelihood$par
            max_likelihood = nlminb(start = C4, objective=minus_ll, lower=0, upper=Inf, C0=C0, C1=C1, C2=C2, C3=C3, C5=C5, subs=1:n, grid=grid, waveletcoefmat=waveletcoefmat)
            C4 = max_likelihood$par
            max_likelihood = nlminb(start = C0, objective=minus_ll, lower=0, upper=Inf, C1=C1, C2=C2, C3=C3, C4=C4, C5=C5, subs=1:n, grid=grid, waveletcoefmat=waveletcoefmat)
            C0 = max_likelihood$par
            max_likelihood = nlminb(start = C1, objective=minus_ll, lower=0, upper=Inf, C0=C0, C2=C2, C3=C3, C4=C4, C5=C5, subs=1:n, grid=grid, waveletcoefmat=waveletcoefmat)
            C1 = max_likelihood$par
            max_likelihood = nlminb(start = C2, objective=minus_ll, lower=0, upper=Inf, C0=C0, C1=C1, C3=C3, C4=C4, C5=C5, subs=1:n, grid=grid, waveletcoefmat=waveletcoefmat)
            C2 = max_likelihood$par
            max_likelihood = nlminb(start = C3, objective=minus_ll, lower=-Inf, upper=Inf, C0=C0, C1=C1, C2=C2, C4=C4, C5=C5, subs=1:n, grid=grid, waveletcoefmat=waveletcoefmat)
            C3 = max_likelihood$par
        }

        # ........... Computing MLE variance estimates................


        h = sqrt(.Machine$double.eps)*c(C0,C1,C2,C3,C4,C5)

        if(C0==0) h[1]= 1.490116e-12
        if(C1==0) h[2]= 1.490116e-12
        if(C2==0) h[3]= 1.490116e-12
        if(C3==0) h[4]= 1.490116e-12
        if(C4==0) h[5]= 1.490116e-12
        if(C5==0) h[6]= 1.490116e-12


        VarMLE = var_mle(C0,C1,C2,C3,C4,C5,n,grid,waveletcoefmat,h)    # Note: taking tol=1e-030. Not setting appropriate tolerance level may show the "system is computationally singular" error.

        return(list(hyperparam=c(C0,C1,C2,C3,C4,C5),hyperparamVar=VarMLE))
    } else
    #
    #
    #
    #
    if(analysis == "single")
    {
        hyperparam = matrix(NA, nrow=n, ncol=6)
        hyperparamVar = array(dim = c(n,6,6))

        for(i in 1:n)
        {
            C0 = C1 = C2 = C3 = C4 = C5 = 1

            for(j in 1:50)
            {
                max_likelihood = nlminb(start = C5, objective=minus_ll, lower=0, upper=Inf, C0=C0, C1=C1, C2=C2, C3=C3, C4=C4, subs=i, grid=grid, waveletcoefmat=waveletcoefmat)
                C5 = max_likelihood$par
                max_likelihood = nlminb(start = C4, objective=minus_ll, lower=0, upper=Inf, C0=C0, C1=C1, C2=C2, C3=C3, C5=C5, subs=i, grid=grid, waveletcoefmat=waveletcoefmat)
                C4 = max_likelihood$par
                max_likelihood = nlminb(start = C0, objective=minus_ll, lower=0, upper=Inf, C1=C1, C2=C2, C3=C3, C4=C4, C5=C5, subs=i, grid=grid, waveletcoefmat=waveletcoefmat)
                C0 = max_likelihood$par
                max_likelihood = nlminb(start = C1, objective=minus_ll, lower=0, upper=Inf, C0=C0, C2=C2, C3=C3, C4=C4, C5=C5, subs=i, grid=grid, waveletcoefmat=waveletcoefmat)
                C1 = max_likelihood$par
                max_likelihood = nlminb(start = C2, objective=minus_ll, lower=0, upper=Inf, C0=C0, C1=C1, C3=C3, C4=C4, C5=C5, subs=i, grid=grid, waveletcoefmat=waveletcoefmat)
                C2 = max_likelihood$par
                max_likelihood = nlminb(start = C3, objective=minus_ll, lower=-Inf, upper=Inf, C0=C0, C1=C1, C2=C2, C4=C4, C5=C5, subs=i, grid=grid, waveletcoefmat=waveletcoefmat)
                C3 = max_likelihood$par
            }


            # ........... Computing MLE variance estimates ..........

            h = sqrt(.Machine$double.eps)*c(C0,C1,C2,C3,C4,C5)
            if(C0==0) h[1]= 1.490116e-12
            if(C1==0) h[2]= 1.490116e-12
            if(C2==0) h[3]= 1.490116e-12
            if(C3==0) h[4]= 1.490116e-12
            if(C4==0) h[5]= 1.490116e-12
            if(C5==0) h[6]= 1.490116e-12

            VarMLE = var_mle(C0,C1,C2,C3,C4,C5,i,grid,waveletcoefmat,h)    # Note: taking tol=1e-030. Not setting appropriate tolerance level may show the "system is computationally singular" error.

            hyperparam[i,] = c(C0,C1,C2,C3,C4,C5)
            hyperparamVar[i,,] = VarMLE
        }

        return(list(hyperparam=hyperparam,hyperparamVar=hyperparamVar))

    }
}

# Function to compute the mixture probabilities, which define the marginal posterior distribution of the wavelet coefficients of the BHMSMA model, for each subject based on multi-subject or single subject analyses (see References)
postmixprob = function(n, grid, waveletcoefmat, hyperparam, analysis="multi")
{
    # ....................... piklj bar by Trapezoidal Rule ...................................

    if(analysis == "multi")
    {
        C0 = hyperparam[1]
        C1 = hyperparam[2]
        C2 = hyperparam[3]
        C3 = hyperparam[4]
        C4 = hyperparam[5]
        C5 = hyperparam[6]

        pkljbar = pklj_bar(grid, n, waveletcoefmat, C0, C1, C2, C3, C4, C5)

        return(list(pkljbar=pkljbar))
    } else
    #
    #
    #
    #
    #
    #
    if(analysis == "single")
    {
        pkljbar = matrix(NA,nrow=n, ncol=grid^2-1)

        for(i in 1:n)
        {
            C0 = hyperparam[i,1]
            C1 = hyperparam[i,2]
            C2 = hyperparam[i,3]
            C3 = hyperparam[i,4]
            C4 = hyperparam[i,5]
            C5 = hyperparam[i,6]
                     
            pkljbar[i,] = pklj_bar(grid, 1, waveletcoefmat[i,,drop=F], C0, C1, C2, C3, C4, C5) 
        }
        return(list(pkljbar=pkljbar))
    } 
}
 
# Function to compute the posterior estimates (mean and median) of the wavelet coefficients of the BHMSMA model for each subject based on multi-subject or single subject analyses (see References)
postwaveletcoef = function( n, grid, waveletcoefmat, hyperparam, pkljbar, analysis="multfif" )
{

    # ......................... Posterior Mean of Wavelet Coeficients ..............................
 
    if(analysis == "multi")
    {
        C4 = hyperparam[5]
        C5 = hyperparam[6]

        out = post_wavelet_coef(grid, n, waveletcoefmat, pkljbar, C4, C5)

        return(out)

    } else
    #
    #
    #
    #
    #
    #
    {
        PostMeanWaveletCoef = matrix( nrow=n, ncol=grid^2-1 )
        PostMedianWaveletCoef = matrix( nrow=n, ncol=grid^2-1 )
        
        for(i in 1:n)
        {
            C4 = hyperparam[i,5]
            C5 = hyperparam[i,6]

            out = post_wavelet_coef(grid, 1, waveletcoefmat[i,,drop=F], pkljbar, C4, C5)

            PostMeanWaveletCoef[i,] = out$PostMeanWaveletCoef[1,]
            PostMedianWaveletCoef[i,] = out$PostMedianWaveletCoef[1,]
        }

        return(list(PostMeanWaveletCoef=PostMeanWaveletCoef, PostMedianWaveletCoef=PostMedianWaveletCoef))
    }
 
}
# Function that substitutes the wavelet coefficients stored in an wavelet object with user given values and returns the modified wavelet object
substituteWaveletCoef = function(grid, waveletobj, values)
{
    out = waveletobj
    out$w0L1 = values[1]
    out$w0L2 = values[2]
    out$w0L3 = values[3]
    if(grid>2)
    {
        out$w1L1 = values[4:7]
        out$w1L2 = values[8:11]
        out$w1L3 = values[12:15]
    }
    if(grid>4)
    {
        out$w2L1 = values[16:31]
        out$w2L2 = values[32:47]
        out$w2L3 = values[48:63]
    }
    if(grid>8)
    {
        out$w3L1 = values[64:127]
        out$w3L2 = values[128:191]
        out$w3L3 = values[192:255]
    }
    if(grid>16)
    {
        out$w4L1 = values[256:511]
        out$w4L2 = values[512:767]
        out$w4L3 = values[768:1023]
    }
    if(grid>32) 
    {
        out$w5L1 = values[1024:2047]
        out$w5L2 = values[2048:3071]
        out$w5L3 = values[3072:4095] 
    }
    if(grid>64) 
    {
        out$w6L1 = values[4096:8191] 
        out$w6L2 = values[8192:12287] 
        out$w6L3 = values[12288:16383]
    }
    if(grid>128)
    {
        out$w7L1 = values[16384:32767]
        out$w7L2 = values[32768:49151]
        out$w7L3 = values[49152:65535]
    }
    if(grid>256)
    {
        out$w8L1 = values[65536:131071]
        out$w8L2 = values[131072:196607]
        out$w8L3 = values[196608:262143]
    }
    out$w0Lconstant = waveletobj$w0Lconstant

    return(out)
}

# Function to compute the posterior estimates (mean and median) of the 2D GLM coefficients map (corresponding to a single brain slice) of the regressor of interest in the BHMSMA model for each subject based on multi-subject or single subject analyses (see References)
postglmcoef = function(n, grid, glmcoefstd, postmeanwaveletcoef, wave.family="DaubLeAsymm", filter.number=6, bc="periodic" )
{
    # ............................. Posterior Mean of GLM coeficients .............................

    PostMeanRecons = array(dim=c(n,grid,grid))
   
    for(i in 1:n)
    {
        dwt = imwd(glmcoefstd[i,,], type="wavelet", family=wave.family, filter.number=filter.number, bc=bc, RetFather=TRUE, verbose=FALSE)

        dwt_new = substituteWaveletCoef(grid,dwt,postmeanwaveletcoef[i,])

        PostMeanRecons[i,,] = imwr(dwt_new)
    }
   return(list(GLMcoefposterior=PostMeanRecons))
}

#1: Function to perform the complete BHMSMA analysis (see References) of a 2D GLM coefficient map (corresponding to a single brain slice) of the regressor of interest based on multi-subject or single subject analyses
BHMSMA = function( n, grid, data, designmat, k, analysis="multi", truecoef=NULL, wave.family="DaubLeAsymm", filter.number=6, bc="periodic")
{

    if(!is.matrix(designmat)) stop("designmat must be a matrix")
    if(n!=nrow(data)) stop(cat("data doesn't have n=",n," rows.\n"))
    if(ncol(designmat)==1) cat("Warning: designmat has only one column. The function doesn't add any intercept column by itself.")
    if(k > ncol(designmat)) stop("Fix the input k.")

    # if( !length( unique(designmat[,k])) == 1)  #if kth column is not an intercept column, add an intercept column
        # designmat = cbind( rep(1,nrow(designmat)), designmat[,k,drop=F] )
        
    glmmap = glmcoef(n, grid, data, designmat)

    wavecoefglmmap = waveletcoef(n, grid, glmmap$GLMCoefStandardized[,,,k], wave.family, filter.number, bc)

    hyperest = hyperparamest( n, grid, wavecoefglmmap$WaveletCoefficientMatrix, analysis )

    pkljbar = postmixprob(n, grid, wavecoefglmmap$WaveletCoefficientMatrix, hyperest$hyperparam, analysis)

    postwavecoefglmmap = postwaveletcoef( n, grid, wavecoefglmmap$WaveletCoefficientMatrix, hyperest$hyperparam, pkljbar$pkljbar, analysis )

    postglmmap = postglmcoef(n, grid, glmmap$GLMCoefStandardized[,,,k], postwavecoefglmmap$PostMeanWaveletCoef, wave.family, filter.number, bc)

    if(! is.null(truecoef))
    {
        MSE = c()
        for(i in 1:n)
        {
            MSE[i] = sum ( ( as.vector(truecoef[i,,]/glmmap$GLMCoefSE[i,,,k]) - as.vector(postglmmap$GLMcoefposterior[i,,]) )^2 )
        }
    }

    if(!is.null(truecoef))
    return(list( GLMCoefStandardized = glmmap$GLMCoefStandardized, GLMCoefSE = glmmap$GLMCoefSE, WaveletCoefficientMatrix = wavecoefglmmap$WaveletCoefficientMatrix, hyperparam = hyperest$hyperparam, hyperparamVar = hyperest$hyperparamVar, posteriorMixProb=pkljbar$pkljbar,  Waveletcoefposterior = postwavecoefglmmap$PostMeanWaveletCoef, GLMcoefposterior = postglmmap$GLMcoefposterior, MSE=MSE))

    if(is.null(truecoef))
    return(list( GLMCoefStandardized = glmmap$GLMCoefStandardized, GLMCoefSE = glmmap$GLMCoefSE, WaveletCoefficientMatrix = wavecoefglmmap$WaveletCoefficientMatrix, hyperparam = hyperest$hyperparam, hyperparamVar = hyperest$hyperparamVar, posteriorMixProb=pkljbar$pkljbar,  Waveletcoefposterior = postwavecoefglmmap$PostMeanWaveletCoef, GLMcoefposterior = postglmmap$GLMcoefposterior))
 
}

# Function to obtain samples from the posterior distribution of a 2D GLM coefficient map (corresponding to a single brain slice) of the regressor of interest in the BHMSMA model for each subject based on multi-subject or single subject analyses (see References)
postsamples = function(nsample, n, grid, glmcoefstd, waveletcoefmat, hyperparam, pkljbar, analysis, wave.family="DaubLeAsymm", filter.number=6, bc="periodic", seed)
{
 
    if(analysis == "multi")
    {
        C4 = hyperparam[5]
        C5 = hyperparam[6] 

        # ................. Simulating  d_iklj..............

        idwt = array( NA, dim=c(n, grid, grid, nsample))
        postdiscovery = array( dim=c(n, grid, grid) )

        d.sim = post_samp(nsample, grid, n, waveletcoefmat, pkljbar, C4, C5, seed)

        for(i in 1:n) 
        {
            dwt = imwd(glmcoefstd[i,,], type="wavelet", family=wave.family, filter.number=filter.number, bc=bc, RetFather=TRUE, verbose=FALSE)            

            for(g in 1:nsample)
            {
                dwt_new = substituteWaveletCoef(grid,dwt,d.sim[i,,g])
                idwt[i,,,g] = imwr(dwt_new) 
            }

            delta = 1
            phi = 0.999

            for(j in 1:grid)
            for(l in 1:grid)
            {
                v = abs(idwt[i,j,l,])
                p =  length(v[v>delta])/length(v)
                if(p>phi)  postdiscovery[i,j,l] = p else postdiscovery[i,j,l] = 0
            }
        }
        return(list(samples = idwt, postdiscovery = postdiscovery))
 
    } else
    #
    #
    #
    #
    #
    #
    if(analysis == "single")
    {
        idwt = array( NA, dim=c(n, grid, grid, nsample))
        postdiscovery = array( dim=c(n, grid, grid) )

        for(i in 1:n)
        {
            C4 = hyperparam[i,5]
            C5 = hyperparam[i,6] 

            dwt = imwd(glmcoefstd[i,,], type="wavelet", family=wave.family, filter.number=filter.number, bc=bc, RetFather=TRUE, verbose=FALSE)

            d.sim = post_samp(nsample, grid, i, waveletcoefmat, pkljbar, C4, C5, seed)

            for(g in 1:nsample)
            {
                dwt_new = substituteWaveletCoef(grid,dwt,d.sim[1,,g])
                idwt[i,,,g] = imwr(dwt_new) 
            }

            delta = 1
            phi = 0.999

            for(j in 1:grid)
            for(l in 1:grid)
            {
                v = abs(idwt[i,j,l,])
                p =  length(v[v>delta])/length(v)
                if(p>phi)  postdiscovery[i,j,l] = p else postdiscovery[i,j,l] = 0
            }

        }
        return(list(samples = idwt, postdiscovery = postdiscovery))
    }

}

# Function to compute the posterior estimates (mean and median) of the 2D GLM coefficient group map (corresponding to a single brain slice) of the regressor of interest in the BHMSMA model based on multi-subject or single subject analyses (see References)
postgroupglmcoef = function( n, grid, glmcoefstd, postmeanwaveletcoef, wave.family="DaubLeAsymm", filter.number=6, bc="periodic" )
{
 
    postmeanwaveletcoef = apply(postmeanwaveletcoef, 2, mean, na.rm=TRUE)

    scaling = c()
    for(i in 1:n)
    {
        dwt = imwd(glmcoefstd[i,,], type="wavelet", family=wave.family, filter.number=filter.number, bc=bc, RetFather=TRUE, verbose=FALSE)
        scaling = c(scaling,dwt$w0Lconstant)
    }

    dwt_new = substituteWaveletCoef(grid,dwt,postmeanwaveletcoef)
    
    dwt_new$w0Lconstant = mean(scaling)

    PostMeanReconsgroup = imwr(dwt_new)

    return(list(groupcoef=PostMeanReconsgroup))

}

# Function to read and import fMRI data from several file types 
readfmridata = function( directory, format, prefix, nimages, dim.image, nii=TRUE )
{
    fmridata = array(NA, dim=c(dim.image,nimages))

    if(format=="Analyze")
    {
        for(i in 1:nimages)
        {
            if (i <= 9) aux = paste0("000",i,".img") 
            if ((i > 9) && (i <= 99)) aux = paste0("00",i,".img") 
            if ((i > 99) && (i <= 999)) aux = paste0("0",i,".img") 
            if (i > 999) aux = paste0(i,".img") 
            a = readANALYZE(paste0(directory, "/", prefix, aux))   
            fmridata[,,,i] = a[,,,1]
            }
        return(fmridata) 
    }

    if(format=="Nifti")
    { 
        for(i in 1:nimages)
        {
            if(nii==TRUE) ext = ".nii" else ext = ".img"
            if (i <= 9) aux = paste0("000",i,ext) 
            if ((i > 9) && (i <= 99)) aux = paste0("00",i,ext) 
            if ((i > 99) && (i <= 999)) aux = paste0("0",i,ext) 
            if (i > 999) aux = paste0(i,ext) 
            fmridata[,,,i] = readNIfTI(paste0(directory, "/", prefix, aux)) 
        }
        return(fmridata) 
    }

    if(format=="Afni")
    {
        fmridata = readAFNI(paste0(directory, "/", prefix))
        return(fmridata) 
    }
 
}




