
pikljbar <- function(nsubject, grid, WaveletCoefficientMatrix, hyperparam, analysis)

{
 # ....................... piklj bar by Trapezoidal Rule ...................................

 if(analysis == "multi")
 {
   C0 <- hyperparam[1]
   C1 <- hyperparam[2]
   C2 <- hyperparam[3]
   C3 <- hyperparam[4]
   C4 <- hyperparam[5]
   C5 <- hyperparam[6]
   # ........................ Writing the loglikelihood ln(p(w_lj/y)) as a function of w_lj 

   LL <- function(w, l, j)
   {
     aux <- 0
     for(subject in 1:nsubject)
       { 
         if (l == 0) d <- WaveletCoefficientMatrix[subject,j]  else d <- WaveletCoefficientMatrix[subject, sum(2^(2*(0:(l-1)))*3) + j]       
         c.l <- C4 * 2^(-C5*l)   
         if (abs(d) <= 39.6) aux1 <- log(  exp( log(w *(1+c.l)^(-.5)) - 0.5*(d^2)/(1+c.l) ) + exp( log(1-w) - 0.5*(d^2) )  ) else aux1 <- log(  exp( log(w *(1+c.l)^(-.5)) - 0.5*(39.6^2)/(1+c.l) ) + exp( log(1-w) - 0.5*(39.6^2) )  )
         aux <- aux + aux1 + (C0*2^(-C1*l) -1)*log(w) + (C2*2^(-C3*l) -1) * log(1-w)
       }
     aux
   }

   # ...................... Writing piklj as a function of w_lj 

   pklj <- function(w,subject,l,j)
   { 
      if (l == 0) d <- WaveletCoefficientMatrix[subject,j]  else d <- WaveletCoefficientMatrix[subject, sum(2^(2*(0:(l-1)))*3) + j] 
      c.l <-  C4 * 2^(-C5*l)
      if(0.5 * c.l/(1+c.l) * d^2 <= 700) aux1 <- exp(0.5 * c.l/(1+c.l) * d^2 ) else aux1 <- exp(700)  # Note: Adjustment
      O <- (1+c.l)^(-0.5) * w/(1-w) * aux1
      p <- O/(1+O)
      p
   }

   # ........................ Using Trapezoidal rule to evaluate p_klj bar 

   N <- 1000
   a <- 0
   b <- 1

   pklj.bar <- matrix(nrow=nsubject, ncol=grid^2-1)
   for(subject in 1:nsubject)
   {
       
     for(l in c(0:(log2(grid)-1)))
     for(j in c(1:(2^(2*l)*3)))
     {   
       w.grid <- a+(0:N)*(b-a)/N
       w.density.unnorm <- exp(LL(w.grid,l,j) )
       for(id in 1:length(w.grid)) {if( is.nan(w.density.unnorm[id]) || (w.density.unnorm[id] == Inf) ) w.density.unnorm[id] <- 1}          
       sum.w <- sum(w.density.unnorm[2:N]) + 0.5 * (w.density.unnorm[1]+w.density.unnorm[N+1])
       w.density.norm <- w.density.unnorm / sum.w

       pklj.func.w <- pklj(w.grid,subject,l,j)
       pklj.func.w[N+1] <- 1

       prod.p.wd <- w.density.norm * pklj.func.w
       integral <-    sum(prod.p.wd[2:N])  + 0.5*(prod.p.wd[1]+prod.p.wd[N+1])   
  
       if (l==0)  pklj.bar[subject,j] <- integral   else pklj.bar[subject,sum(2^(2*(0:(l-1)))*3) + j] <- integral
     } 
     
   }
   return(list(pklj.bar=pklj.bar))
 }
 
 
 if(analysis == "single")
 {

   pklj.bar <- matrix(NA,nrow=nsubject, ncol=grid^2-1)
   
   for(subject in 1:nsubject)
   {
     C0 <- hyperparam[subject,1]
     C1 <- hyperparam[subject,2]
     C2 <- hyperparam[subject,3]
     C3 <- hyperparam[subject,4]
     C4 <- hyperparam[subject,5]
     C5 <- hyperparam[subject,6]
     
     # ........................ Writing the loglikelihood ln(p(w_lj/y)) as a function of w_lj 

     LL <- function(w, l, j)
     {
       aux <- 0
       if (l == 0) d <- WaveletCoefficientMatrix[subject,j]  else d <- WaveletCoefficientMatrix[subject, sum(2^(2*(0:(l-1)))*3) + j]       
       c.l <- C4 * 2^(-C5*l)   
       if (abs(d) <= 39.6) aux1 <- log(  exp( log(w *(1+c.l)^(-.5)) - 0.5*(d^2)/(1+c.l) ) + exp( log(1-w) - 0.5*(d^2) )  ) else aux1 <- log(  exp( log(w *(1+c.l)^(-.5)) - 0.5*(39.6^2)/(1+c.l) ) + exp( log(1-w) - 0.5*(39.6^2) )  )
       aux <- aux + aux1 + (C0*2^(-C1*l) -1)*log(w) + (C2*2^(-C3*l) -1) * log(1-w)
       aux
     }
     
     # ...................... Writing pklj as a function of w_lj 

     pklj <- function(w, subject, l, j)
     { 
        if (l == 0) d <- WaveletCoefficientMatrix[subject,j]  else d <- WaveletCoefficientMatrix[subject, sum(2^(2*(0:(l-1)))*3) + j] 
        c.l <-  C4 * 2^(-C5*l)
        if(0.5 * c.l/(1+c.l) * d^2 <= 700) aux1 <- exp(0.5 * c.l/(1+c.l) * d^2 ) else aux1 <- exp(700)  # Note: Adjustment
        O <- (1+c.l)^(-0.5) * w/(1-w) * aux1
        p <- O/(1+O)
        p
     }

     # ........................ Using Trapezoidal rule to evaluate p_klj bar 

     N <- 1000
     a <- 0
     b <- 1

     for(l in c(0:(log2(grid)-1)))
     for(j in c(1:(2^(2*l)*3)))
     {   
       w.grid <- a+(0:N)*(b-a)/N
       w.density.unnorm <- exp(LL(w.grid,l,j) )
       for(id in 1:length(w.grid)) {if( is.nan(w.density.unnorm[id]) || (w.density.unnorm[id] == Inf) ) w.density.unnorm[id] <- 1}          
       sum.w <- sum(w.density.unnorm[2:N]) + 0.5 * (w.density.unnorm[1]+w.density.unnorm[N+1])
       w.density.norm <- w.density.unnorm / sum.w

       pklj.func.w <- pklj(w.grid,subject,l,j)
       pklj.func.w[N+1] <- 1

       prod.p.wd <- w.density.norm * pklj.func.w
       integral <-    sum(prod.p.wd[2:N])  + 0.5*(prod.p.wd[1]+prod.p.wd[N+1])   
  
       if (l==0)  pklj.bar[subject,j] <- integral   else pklj.bar[subject,sum(2^(2*(0:(l-1)))*3) + j] <- integral
     }  
   }
   return(list(pklj.bar=pklj.bar))
 }
 
}