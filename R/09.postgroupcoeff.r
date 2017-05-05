# .......................... Group map (MSA) ..................................

postgroupcoeff = function( nsubject, grid, GLMCoeffStandardized, PostMeanWaveletCoeff, wave.family="DaubLeAsymm", filter.number=6, bc="periodic" )
{
 
 PostMeanWaveletCoeff = apply(PostMeanWaveletCoeff, 2, mean, na.rm=TRUE)
 
 scaling = c()
 for(subject in 1:nsubject)
 {
  dwt = imwd(GLMCoeffStandardized[subject,,], type="wavelet", family=wave.family, filter.number=filter.number, bc=bc, RetFather=TRUE, verbose=FALSE)
  scaling = c(scaling,dwt$w0Lconstant)
 }

 w.test1 = dwt
 w.test1$w0L1 = PostMeanWaveletCoeff[1]
 w.test1$w0L2 = PostMeanWaveletCoeff[2]
 w.test1$w0L3 = PostMeanWaveletCoeff[3]
 if(grid>2)
 {
 w.test1$w1L1 = PostMeanWaveletCoeff[4:7]
 w.test1$w1L2 = PostMeanWaveletCoeff[8:11]
 w.test1$w1L3 = PostMeanWaveletCoeff[12:15]
 }
 if(grid>4)
 {
 w.test1$w2L1 = PostMeanWaveletCoeff[16:31]
 w.test1$w2L2 = PostMeanWaveletCoeff[32:47]
 w.test1$w2L3 = PostMeanWaveletCoeff[48:63]
 }
 if(grid>8)
 {
 w.test1$w3L1 = PostMeanWaveletCoeff[64:127]
 w.test1$w3L2 = PostMeanWaveletCoeff[128:191]
 w.test1$w3L3 = PostMeanWaveletCoeff[192:255]
 }
 if(grid>16)
 {
 w.test1$w4L1 = PostMeanWaveletCoeff[256:511]
 w.test1$w4L2 = PostMeanWaveletCoeff[512:767]
 w.test1$w4L3 = PostMeanWaveletCoeff[768:1023]
 }
 if(grid>32) 
 {
 w.test1$w5L1 = PostMeanWaveletCoeff[1024:2047]
 w.test1$w5L2 = PostMeanWaveletCoeff[2048:3071]
 w.test1$w5L3 = PostMeanWaveletCoeff[3072:4095] 
 }
 if(grid>64) 
 {
 w.test1$w6L1 = PostMeanWaveletCoeff[4096:8191] 
 w.test1$w6L2 = PostMeanWaveletCoeff[8192:12287] 
 w.test1$w6L3 = PostMeanWaveletCoeff[12288:16383]
 }
 if(grid>128)
 {
 w.test1$w7L1 = PostMeanWaveletCoeff[16384:32767]
 w.test1$w7L2 = PostMeanWaveletCoeff[32768:49151]
 w.test1$w7L3 = PostMeanWaveletCoeff[49152:65535]
 }
 if(grid>256)
 {
 w.test1$w8L1 = PostMeanWaveletCoeff[65536:131071]
 w.test1$w8L2 = PostMeanWaveletCoeff[131072:196607]
 w.test1$w8L3 = PostMeanWaveletCoeff[196608:262143]
 }
 w.test1$w0Lconstant = mean(scaling)

 PostMeanReconsgroup = imwr(w.test1)
 
 return(list(groupcoeff=PostMeanReconsgroup))

}