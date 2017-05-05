
postsamples = function(nsample, nsubject, grid, GLMCoeffStandardized, WaveletCoefficientMatrix, hyperparam, pklj.bar, analysis, wave.family="DaubLeAsymm", filter.number=6, bc="periodic")
{
 
 if(analysis == "multi")
 {
 
   C4 = hyperparam[5]
   C5 = hyperparam[6] 
   
   # ................. Simulating  d_iklj..............

   idwt = array( NA, dim=c(nsubject, grid, grid, nsample))
   postdiscovery = array( dim=c(nsubject, grid, grid) )
   
   for(subject in 1:nsubject) 
   {
     dwt = imwd(GLMCoeffStandardized[subject,,], type="wavelet", family=wave.family, filter.number=filter.number, bc=bc, RetFather=TRUE, verbose=FALSE)
      
     d.sim = array(NA,dim=c(grid^2-1, nsample))
     
     for(g in 1:nsample)
     {
       for(l in c(0:(log2(grid)-1)))
       for(j in c(1:(2^(2*l)*3)))
       { 
        if(l==0) d = WaveletCoefficientMatrix[subject, j]  else d = WaveletCoefficientMatrix[subject, sum(2^(2*(0:(l-1)))*3) + j]
        if(l==0) p = pklj.bar[subject, j] else p =  pklj.bar[subject, sum(2^(2*(0:(l-1)))*3) + j]   
        c_kl = C4 * 2^(-C5*l)
        U = rbinom(1, size=1, prob=min(p,1))
        if (l==0) {if (U==1) d.sim[j,g] = rnorm( 1, mean = c_kl/(1+c_kl) * d , sd = sqrt(c_kl/(1+c_kl)) )  else d.sim[j,g] = 0 } else {if (U==1) d.sim[sum(2^(2*(0:(l-1)))*3) + j,g] = rnorm( 1, mean = c_kl/(1+c_kl) * d , sd = sqrt(c_kl/(1+c_kl)) )  else d.sim[sum(2^(2*(0:(l-1)))*3) + j,g] = 0 }
       }

       w.test1 = dwt
       w.test1$w0L1 = d.sim[1,g]
       w.test1$w0L2 = d.sim[2,g]
       w.test1$w0L3 = d.sim[3,g]
       if(grid>2)
       {
       w.test1$w1L1 = d.sim[4:7,g]
       w.test1$w1L2 = d.sim[8:11,g]
       w.test1$w1L3 = d.sim[12:15,g]
       }
       if(grid>4)
       {
       w.test1$w2L1 = d.sim[16:31,g]
       w.test1$w2L2 = d.sim[32:47,g]
       w.test1$w2L3 = d.sim[48:63,g]
       }
       if(grid>8)
       {
       w.test1$w3L1 = d.sim[64:127,g]
       w.test1$w3L2 = d.sim[128:191,g]
       w.test1$w3L3 = d.sim[192:255,g]
       }
       if(grid>16)
       {
       w.test1$w4L1 = d.sim[256:511,g]
       w.test1$w4L2 = d.sim[512:767,g]
       w.test1$w4L3 = d.sim[768:1023,g]
       }
       if(grid>32)
       {
       w.test1$w5L1 = d.sim[1024:2047,g]
       w.test1$w5L2 = d.sim[2048:3071,g]
       w.test1$w5L3 = d.sim[3072:4095,g]
       }
       if(grid>64) 
       {
       w.test1$w6L1 = d.sim[4096:8191,g] 
       w.test1$w6L2 = d.sim[8192:12287,g] 
       w.test1$w6L3 = d.sim[12288:16383,g]
       }
       if(grid>128)
       {
       w.test1$w7L1 = d.sim[16384:32767,g]
       w.test1$w7L2 = d.sim[32768:49151,g]
       w.test1$w7L3 = d.sim[49152:65535,g]
       }
       if(grid>256)
       {
       w.test1$w8L1 = d.sim[65536:131071,g]
       w.test1$w8L2 = d.sim[131072:196607,g]
       w.test1$w8L3 = d.sim[196608:262143,g]
       }
       w.test1$w0Lconstant = dwt$w0Lconstant  
       
       idwt[subject,,,g] = imwr(w.test1) 
     }
     
     delta = 1
     phi = 0.999
     
     for(i in 1:grid)
     for(j in 1:grid)
     {
       v = abs(idwt[subject,i,j,])
       p =  length(v[v>delta])/length(v)
       if(p>phi)  postdiscovery[subject,i,j] = p else postdiscovery[subject,i,j] = 0
     }
   }
   return(list(samples = idwt, postdiscovery = postdiscovery))
 }
 
 if(analysis == "single")
 {
   idwt = array( NA, dim=c(nsubject, grid, grid, nsample))
   postdiscovery = array( dim=c(nsubject, grid, grid) )
   
   for(subject in 1:nsubject)
   {
     C4 = hyperparam[subject,5]
     C5 = hyperparam[subject,6] 
   
     dwt = imwd(GLMCoeffStandardized[subject,,], type="wavelet", family=wave.family, filter.number=filter.number, bc=bc, RetFather=TRUE, verbose=FALSE)
      
     d.sim = array(NA,dim=c(grid^2-1, nsample))
     
     for(g in 1:nsample)
     {
       for(l in c(0:(log2(grid)-1)))
       for(j in c(1:(2^(2*l)*3)))
       { 
         if(l==0) d = WaveletCoefficientMatrix[subject, j]  else d = WaveletCoefficientMatrix[subject, sum(2^(2*(0:(l-1)))*3) + j]
         if(l==0) p = pklj.bar[subject, j] else p =  pklj.bar[subject, sum(2^(2*(0:(l-1)))*3) + j]   
         c_kl = C4 * 2^(-C5*l)
         U = rbinom(1, size=1, prob=min(p,1))
         if (l==0) {if (U==1) d.sim[j,g] = rnorm( 1, mean = c_kl/(1+c_kl) * d , sd = sqrt(c_kl/(1+c_kl)) )  else d.sim[j,g] = 0 } else {if (U==1) d.sim[sum(2^(2*(0:(l-1)))*3) + j,g] = rnorm( 1, mean = c_kl/(1+c_kl) * d , sd = sqrt(c_kl/(1+c_kl)) )  else d.sim[sum(2^(2*(0:(l-1)))*3) + j,g] = 0 }
       }

       w.test1 = dwt
       w.test1$w0L1 = d.sim[1,g]
       w.test1$w0L2 = d.sim[2,g]
       w.test1$w0L3 = d.sim[3,g]
       if(grid>2)
       {
       w.test1$w1L1 = d.sim[4:7,g]
       w.test1$w1L2 = d.sim[8:11,g]
       w.test1$w1L3 = d.sim[12:15,g]
       }
       if(grid>4)
       {
       w.test1$w2L1 = d.sim[16:31,g]
       w.test1$w2L2 = d.sim[32:47,g]
       w.test1$w2L3 = d.sim[48:63,g]
       }
       if(grid>8)
       {
       w.test1$w3L1 = d.sim[64:127,g]
       w.test1$w3L2 = d.sim[128:191,g]
       w.test1$w3L3 = d.sim[192:255,g]
       }
       if(grid>16)
       {
       w.test1$w4L1 = d.sim[256:511,g]
       w.test1$w4L2 = d.sim[512:767,g]
       w.test1$w4L3 = d.sim[768:1023,g]
       }
       if(grid>32)
       {
       w.test1$w5L1 = d.sim[1024:2047,g]
       w.test1$w5L2 = d.sim[2048:3071,g]
       w.test1$w5L3 = d.sim[3072:4095,g]
       }
       if(grid>64) 
       {
       w.test1$w6L1 = d.sim[4096:8191,g] 
       w.test1$w6L2 = d.sim[8192:12287,g] 
       w.test1$w6L3 = d.sim[12288:16383,g]
       }
       if(grid>128)
       {
       w.test1$w7L1 = d.sim[16384:32767,g]
       w.test1$w7L2 = d.sim[32768:49151,g]
       w.test1$w7L3 = d.sim[49152:65535,g]
       }
       if(grid>256)
       {
       w.test1$w8L1 = d.sim[65536:131071,g]
       w.test1$w8L2 = d.sim[131072:196607,g]
       w.test1$w8L3 = d.sim[196608:262143,g]
       }
       w.test1$w0Lconstant = dwt$w0Lconstant  
       
       idwt[subject,,,g] = imwr(w.test1) 
     }
     
     delta = 1
     phi = 0.999
     
     for(i in 1:grid)
     for(j in 1:grid)
     {
       v = abs(idwt[subject,i,j,])
       p =  length(v[v>delta])/length(v)
       if(p>phi)  postdiscovery[subject,i,j] = p else postdiscovery[subject,i,j] = 0
     }
   
   }
   return(list(samples = idwt, postdiscovery = postdiscovery))
 }

}