read.fmridata <- function( directory, format, prefix, nimages, dim.image, nii=TRUE )
{
 
 fmridata <- array(NA, dim=c(dim.image,nimages))
 
 if(format=="Analyze")
 {
   for(i in 1:nimages)
   {
    if (i <= 9) aux <- paste0("000",i,".img") 
    if ((i > 9) && (i <= 99)) aux <- paste0("00",i,".img") 
    if ((i > 99) && (i <= 999)) aux <- paste0("0",i,".img") 
    if (i > 999) aux <- paste0(i,".img") 
    a <- f.read.analyze.volume(paste0(directory, "/", prefix, aux))   
    fmridata[,,,i] <- a[,,,1]
   }
   return(fmridata) 
 }
 
 if(format=="Nifti")
 { 
   for(i in 1:nimages)
   {
    if(nii==TRUE) ext <- ".nii" else ext <- ".img"
    if (i <= 9) aux <- paste0("000",i,ext) 
    if ((i > 9) && (i <= 99)) aux <- paste0("00",i,ext) 
    if ((i > 99) && (i <= 999)) aux <- paste0("0",i,ext) 
    if (i > 999) aux <- paste0(i,ext) 
    a <- f.read.nifti.volume(paste0(directory, "/", prefix, aux))   
    fmridata[,,,i] <- a
   }
   return(fmridata) 
 }
 
 if(format=="Afni")
 {
   a <- read.AFNI(paste0(directory, "/", prefix))   
   fmridata <- extract.data(a, what="data")
   return(fmridata) 
 }
 
}

