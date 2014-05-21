
hyperparamest <- function( nsubject, grid, WaveletCoefficientMatrix, analysis )

{

 #....................... Estimating C0, C1, C2, C3, C4 and C5..................
 
 if(analysis == "multi")
 {

   minus_loglikelihood_function <- function(C0,C1,C2,C3,C4,C5)
   {
     aux <- 0
     for(subject in c(1:nsubject))
     for(l in c(0:(log2(grid)-1)))
     for(j in c(1:(2^(2*l)*3)))
     { 
       if(l==0) d <- WaveletCoefficientMatrix[subject,j]  else d <- WaveletCoefficientMatrix[subject, sum(2^(2*(0:(l-1)))*3) + j]
       a.l <- C0 * 2^(-C1*l)
       b.l <- C2 * 2^(-C3*l)
       expect.w <- min(1, a.l/(a.l+b.l) )
       c.l <- C4 * 2^(-C5*l)       
       if (abs(d)<=39.6) aux1 <- log(  exp( log(expect.w *(1+c.l)^(-.5)) - 0.5*(d^2)/(1+c.l) ) + exp( log(1-expect.w) - 0.5*(d^2) )  )  else aux1 <- log(  exp( log(expect.w *(1+c.l)^(-.5)) - 0.5*(39.6^2)/(1+c.l) ) + exp( log(1-expect.w) - 0.5*(39.6^2) )  )
       aux <- aux + aux1
     }
     -aux
   }

   C0 <- 1
   C1 <- 1
   C2 <- 1
   C3 <- 1
   C4 <- 1
   C5 <- 1

   for(i in 1:50)
   {
     max_likelihood <- nlminb(start = C5, objective=minus_loglikelihood_function, lower=0, upper=Inf, C0=C0, C1=C1, C2=C2, C3=C3, C4=C4)
     C5 <- max_likelihood$par
     max_likelihood <- nlminb(start = C4, objective=minus_loglikelihood_function, lower=0, upper=Inf, C0=C0, C1=C1, C2=C2, C3=C3, C5=C5)
     C4 <- max_likelihood$par
     max_likelihood <- nlminb(start = C0, objective=minus_loglikelihood_function, lower=0, upper=Inf, C1=C1, C2=C2, C3=C3, C4=C4, C5=C5)
     C0 <- max_likelihood$par
     max_likelihood <- nlminb(start = C1, objective=minus_loglikelihood_function, lower=0, upper=Inf, C0=C0, C2=C2, C3=C3, C4=C4, C5=C5)
     C1 <- max_likelihood$par
     max_likelihood <- nlminb(start = C2, objective=minus_loglikelihood_function, lower=0, upper=Inf, C0=C0, C1=C1, C3=C3, C4=C4, C5=C5)
     C2 <- max_likelihood$par
     max_likelihood <- nlminb(start = C3, objective=minus_loglikelihood_function, lower=-Inf, upper=Inf, C0=C0, C1=C1, C2=C2, C4=C4, C5=C5)
     C3 <- max_likelihood$par
   }

   # ........... Computing MLE variance estimates................

   ll <- function(C0,C1,C2,C3,C4,C5)
   {
    - minus_loglikelihood_function(C0,C1,C2,C3,C4,C5)
   }

   IM <- matrix(0,nrow=6,ncol=6)

   h <- sqrt(.Machine$double.eps)*c(C0,C1,C2,C3,C4,C5)
   if(C0==0) h[1]<- 1.490116e-12
   if(C1==0) h[2]<- 1.490116e-12
   if(C2==0) h[3]<- 1.490116e-12
   if(C3==0) h[4]<- 1.490116e-12
   if(C4==0) h[5]<- 1.490116e-12
   if(C5==0) h[6]<- 1.490116e-12


   IM[1,1] <- -( ll(C0+h[1],C1,C2,C3,C4,C5) - ll(C0,C1,C2,C3,C4,C5) + ll(C0-h[1],C1,C2,C3,C4,C5) ) / (h[1]^2)
   IM[1,2] <- -( ll(C0+h[1],C1+h[2],C2,C3,C4,C5) - ll(C0+h[1],C1-h[2],C2,C3,C4,C5) - ll(C0-h[1],C1+h[2],C2,C3,C4,C5) + ll(C0-h[1],C1-h[2],C2,C3,C4,C5) ) / (4*h[1]*h[2])
   IM[1,3] <- -( ll(C0+h[1],C1,C2+h[3],C3,C4,C5) - ll(C0+h[1],C1,C2-h[3],C3,C4,C5) - ll(C0-h[1],C1,C2+h[3],C3,C4,C5) + ll(C0-h[1],C1,C2-h[3],C3,C4,C5) ) / (4*h[1]*h[3])
   IM[1,4] <- -( ll(C0+h[1],C1,C2,C3+h[4],C4,C5) - ll(C0+h[1],C1,C2,C3-h[4],C4,C5) - ll(C0-h[1],C1,C2,C3+h[4],C4,C5) + ll(C0-h[1],C1,C2,C3-h[4],C4,C5) ) / (4*h[1]*h[4])
   IM[1,5] <- -( ll(C0+h[1],C1,C2,C3,C4+h[5],C5) - ll(C0+h[1],C1,C2,C3,C4-h[5],C5) - ll(C0-h[1],C1,C2,C3,C4+h[5],C5) + ll(C0-h[1],C1,C2,C3,C4-h[5],C5) ) / (4*h[1]*h[5])
   IM[1,6] <- -( ll(C0+h[1],C1,C2,C3,C4,C5+h[6]) - ll(C0+h[1],C1,C2,C3,C4,C5-h[6]) - ll(C0-h[1],C1,C2,C3,C4,C5+h[6]) + ll(C0-h[1],C1,C2,C3,C4,C5-h[6]) ) / (4*h[1]*h[6])
   IM[2,2] <- -( ll(C0,C1+h[2],C2,C3,C4,C5) - ll(C0,C1,C2,C3,C4,C5) + ll(C0,C1-h[2],C2,C3,C4,C5) ) / (h[2]^2)
   IM[2,3] <- -( ll(C0,C1+h[2],C2+h[3],C3,C4,C5) - ll(C0,C1+h[2],C2-h[3],C3,C4,C5) - ll(C0,C1-h[2],C2+h[3],C3,C4,C5) + ll(C0,C1-h[2],C2-h[3],C3,C4,C5) ) / (4*h[2]*h[3])
   IM[2,4] <- -( ll(C0,C1+h[2],C2,C3+h[4],C4,C5) - ll(C0,C1+h[2],C2,C3-h[4],C4,C5) - ll(C0,C1-h[2],C2,C3+h[4],C4,C5) + ll(C0,C1-h[2],C2,C3-h[4],C4,C5) ) / (4*h[2]*h[4])
   IM[2,5] <- -( ll(C0,C1+h[2],C2,C3,C4+h[5],C5) - ll(C0,C1+h[2],C2,C3,C4-h[5],C5) - ll(C0,C1-h[2],C2,C3,C4+h[5],C5) + ll(C0,C1-h[2],C2,C3,C4-h[5],C5) ) / (4*h[2]*h[5])
   IM[2,6] <- -( ll(C0,C1+h[2],C2,C3,C4,C5+h[6]) - ll(C0,C1+h[2],C2,C3,C4,C5-h[6]) - ll(C0,C1-h[2],C2,C3,C4,C5+h[6]) + ll(C0,C1-h[2],C2,C3,C4,C5-h[6]) ) / (4*h[2]*h[6])
   IM[3,3] <- -( ll(C0,C1,C2+h[3],C3,C4,C5) - ll(C0,C1,C2,C3,C4,C5) + ll(C0,C1,C2-h[3],C3,C4,C5) ) / (h[3]^2)
   IM[3,4] <- -( ll(C0,C1,C2+h[3],C3+h[4],C4,C5) - ll(C0,C1,C2+h[3],C3-h[4],C4,C5) - ll(C0,C1,C2-h[3],C3+h[4],C4,C5) + ll(C0,C1,C2-h[3],C3-h[4],C4,C5) ) / (4*h[3]*h[4])
   IM[3,5] <- -( ll(C0,C1,C2+h[3],C3,C4+h[5],C5) - ll(C0,C1,C2+h[3],C3,C4-h[5],C5) - ll(C0,C1,C2-h[3],C3,C4+h[5],C5) + ll(C0,C1,C2-h[3],C3,C4-h[5],C5) ) / (4*h[3]*h[5])
   IM[3,6] <- -( ll(C0,C1,C2+h[3],C3,C4,C5+h[6]) - ll(C0,C1,C2+h[3],C3,C4,C5-h[6]) - ll(C0,C1,C2-h[3],C3,C4,C5+h[6]) + ll(C0,C1,C2-h[3],C3,C4,C5-h[6]) ) / (4*h[3]*h[6])
   IM[4,4] <- -( ll(C0,C1,C2,C3+h[4],C4,C5) - ll(C0,C1,C2,C3,C4,C5) + ll(C0,C1,C2,C3-h[4],C4,C5) ) / (h[4]^2)
   IM[4,5] <- -( ll(C0,C1,C2,C3+h[4],C4+h[5],C5) - ll(C0,C1,C2,C3+h[4],C4-h[5],C5) - ll(C0,C1,C2,C3-h[4],C4+h[5],C5) + ll(C0,C1,C2,C3-h[4],C4-h[5],C5) ) / (4*h[4]*h[5])
   IM[4,6] <- -( ll(C0,C1,C2,C3+h[4],C4,C5+h[6]) - ll(C0,C1,C2,C3+h[4],C4,C5-h[6]) - ll(C0,C1,C2,C3-h[4],C4,C5+h[6]) + ll(C0,C1,C2,C3-h[4],C4,C5-h[6]) ) / (4*h[4]*h[6])
   IM[5,5] <- -( ll(C0,C1,C2,C3,C4+h[5],C5) - ll(C0,C1,C2,C3,C4,C5) + ll(C0,C1,C2,C3,C4-h[5],C5) ) / (h[5]^2)
   IM[5,6] <- -( ll(C0,C1,C2,C3,C4+h[5],C5+h[6]) - ll(C0,C1,C2,C3,C4+h[5],C5-h[6]) - ll(C0,C1,C2,C3,C4-h[5],C5+h[6]) + ll(C0,C1,C2,C3,C4-h[5],C5-h[6]) ) / (4*h[5]*h[6])
   IM[6,6] <- -( ll(C0,C1,C2,C3,C4,C5+h[6]) - ll(C0,C1,C2,C3,C4,C5) + ll(C0,C1,C2,C3,C4,C5-h[6]) ) / (h[6]^2)

   IM[2,1] <- IM[1,2]
   IM[3,1] <- IM[1,3]
   IM[3,2] <- IM[2,3]
   IM[4,1] <- IM[1,4]
   IM[4,2] <- IM[2,4]
   IM[4,3] <- IM[3,4]
   IM[5,1] <- IM[1,5]
   IM[5,2] <- IM[2,5]
   IM[5,3] <- IM[3,5]
   IM[5,4] <- IM[4,5]
   IM[6,1] <- IM[1,6]
   IM[6,2] <- IM[2,6]
   IM[6,3] <- IM[3,6]
   IM[6,4] <- IM[4,6]
   IM[6,5] <- IM[5,6]
   
   VarMLE <- solve(IM, tol=1e-030)    # Note: taking tol=1e-030. Not setting appropriate tolerance level may show the "system is computationally singular" error.

   hyperparam <- c(C0,C1,C2,C3,C4,C5)
   return(list(hyperparam=hyperparam,hyperparamVar=VarMLE))
 }

 
 
 if(analysis == "single")
 {
   minus_loglikelihood_function <- function(C0,C1,C2,C3,C4,C5)
   {
     aux <- 0
     for(l in c(0:(log2(grid)-1)))
     for(j in c(1:(2^(2*l)*3)))
     { 
       if(l==0) d <- WaveletCoefficientMatrix[subject,j]  else d <- WaveletCoefficientMatrix[subject, sum(2^(2*(0:(l-1)))*3) + j]
       a.l <- C0*2^(-C1*l)
       b.l <- C2*2^(-C3*l)
       expect.w <- min(1, a.l/(a.l+b.l) )
       c.l <- C4*2^(-C5*l)       
       if (abs(d)<=39.6) aux1 <- log(  exp( log(expect.w *(1+c.l)^(-.5)) - 0.5*(d^2)/(1+c.l) ) + exp( log(1-expect.w) - 0.5*(d^2) )  )  else aux1 <- log(  exp( log(expect.w *(1+c.l)^(-.5)) - 0.5*(39.6^2)/(1+c.l) ) + exp( log(1-expect.w) - 0.5*(39.6^2) )  )
       aux <- aux + aux1
     }
     -aux
   }

   C0 <- 1
   C1 <- 1
   C2 <- 1
   C3 <- 1
   C4 <- 1
   C5 <- 1
   
   hyperparam <- matrix(NA,nrow=nsubject,ncol=6)
   hyperparamVar <- array(dim=c(nsubject,6,6))
   for(subject in 1:nsubject)
   {
     for(i in 1:50)
     {
       max_likelihood <- nlminb(start = C5, objective=minus_loglikelihood_function, lower=0, upper=Inf, C0=C0, C1=C1, C2=C2, C3=C3, C4=C4)
       C5 <- max_likelihood$par
       max_likelihood <- nlminb(start = C4, objective=minus_loglikelihood_function, lower=0, upper=Inf, C0=C0, C1=C1, C2=C2, C3=C3, C5=C5)
       C4 <- max_likelihood$par
       max_likelihood <- nlminb(start = C0, objective=minus_loglikelihood_function, lower=0, upper=Inf, C1=C1, C2=C2, C3=C3, C4=C4, C5=C5)
       C0 <- max_likelihood$par
       max_likelihood <- nlminb(start = C1, objective=minus_loglikelihood_function, lower=0, upper=Inf, C0=C0, C2=C2, C3=C3, C4=C4, C5=C5)
       C1 <- max_likelihood$par
       max_likelihood <- nlminb(start = C2, objective=minus_loglikelihood_function, lower=0, upper=Inf, C0=C0, C1=C1, C3=C3, C4=C4, C5=C5)
       C2 <- max_likelihood$par
       max_likelihood <- nlminb(start = C3, objective=minus_loglikelihood_function, lower=-Inf, upper=Inf, C0=C0, C1=C1, C2=C2, C4=C4, C5=C5)
       C3 <- max_likelihood$par
     }
   
   
   # ........... Computing MLE variance estimates ..........

   ll <- function(C0,C1,C2,C3,C4,C5)
   {
    - minus_loglikelihood_function(C0,C1,C2,C3,C4,C5)
   }
   
   IM <- matrix(0,nrow=6,ncol=6)

   h <- sqrt(.Machine$double.eps)*c(C0,C1,C2,C3,C4,C5)
   if(C0==0) h[1]<- 1.490116e-12
   if(C1==0) h[2]<- 1.490116e-12
   if(C2==0) h[3]<- 1.490116e-12
   if(C3==0) h[4]<- 1.490116e-12
   if(C4==0) h[5]<- 1.490116e-12
   if(C5==0) h[6]<- 1.490116e-12

   IM[1,1] <- -( ll(C0+h[1],C1,C2,C3,C4,C5) - ll(C0,C1,C2,C3,C4,C5) + ll(C0-h[1],C1,C2,C3,C4,C5) ) / (h[1]^2)
   IM[1,2] <- -( ll(C0+h[1],C1+h[2],C2,C3,C4,C5) - ll(C0+h[1],C1-h[2],C2,C3,C4,C5) - ll(C0-h[1],C1+h[2],C2,C3,C4,C5) + ll(C0-h[1],C1-h[2],C2,C3,C4,C5) ) / (4*h[1]*h[2])
   IM[1,3] <- -( ll(C0+h[1],C1,C2+h[3],C3,C4,C5) - ll(C0+h[1],C1,C2-h[3],C3,C4,C5) - ll(C0-h[1],C1,C2+h[3],C3,C4,C5) + ll(C0-h[1],C1,C2-h[3],C3,C4,C5) ) / (4*h[1]*h[3])
   IM[1,4] <- -( ll(C0+h[1],C1,C2,C3+h[4],C4,C5) - ll(C0+h[1],C1,C2,C3-h[4],C4,C5) - ll(C0-h[1],C1,C2,C3+h[4],C4,C5) + ll(C0-h[1],C1,C2,C3-h[4],C4,C5) ) / (4*h[1]*h[4])
   IM[1,5] <- -( ll(C0+h[1],C1,C2,C3,C4+h[5],C5) - ll(C0+h[1],C1,C2,C3,C4-h[5],C5) - ll(C0-h[1],C1,C2,C3,C4+h[5],C5) + ll(C0-h[1],C1,C2,C3,C4-h[5],C5) ) / (4*h[1]*h[5])
   IM[1,6] <- -( ll(C0+h[1],C1,C2,C3,C4,C5+h[6]) - ll(C0+h[1],C1,C2,C3,C4,C5-h[6]) - ll(C0-h[1],C1,C2,C3,C4,C5+h[6]) + ll(C0-h[1],C1,C2,C3,C4,C5-h[6]) ) / (4*h[1]*h[6])
   IM[2,2] <- -( ll(C0,C1+h[2],C2,C3,C4,C5) - ll(C0,C1,C2,C3,C4,C5) + ll(C0,C1-h[2],C2,C3,C4,C5) ) / (h[2]^2)
   IM[2,3] <- -( ll(C0,C1+h[2],C2+h[3],C3,C4,C5) - ll(C0,C1+h[2],C2-h[3],C3,C4,C5) - ll(C0,C1-h[2],C2+h[3],C3,C4,C5) + ll(C0,C1-h[2],C2-h[3],C3,C4,C5) ) / (4*h[2]*h[3])
   IM[2,4] <- -( ll(C0,C1+h[2],C2,C3+h[4],C4,C5) - ll(C0,C1+h[2],C2,C3-h[4],C4,C5) - ll(C0,C1-h[2],C2,C3+h[4],C4,C5) + ll(C0,C1-h[2],C2,C3-h[4],C4,C5) ) / (4*h[2]*h[4])
   IM[2,5] <- -( ll(C0,C1+h[2],C2,C3,C4+h[5],C5) - ll(C0,C1+h[2],C2,C3,C4-h[5],C5) - ll(C0,C1-h[2],C2,C3,C4+h[5],C5) + ll(C0,C1-h[2],C2,C3,C4-h[5],C5) ) / (4*h[2]*h[5])
   IM[2,6] <- -( ll(C0,C1+h[2],C2,C3,C4,C5+h[6]) - ll(C0,C1+h[2],C2,C3,C4,C5-h[6]) - ll(C0,C1-h[2],C2,C3,C4,C5+h[6]) + ll(C0,C1-h[2],C2,C3,C4,C5-h[6]) ) / (4*h[2]*h[6])
   IM[3,3] <- -( ll(C0,C1,C2+h[3],C3,C4,C5) - ll(C0,C1,C2,C3,C4,C5) + ll(C0,C1,C2-h[3],C3,C4,C5) ) / (h[3]^2)
   IM[3,4] <- -( ll(C0,C1,C2+h[3],C3+h[4],C4,C5) - ll(C0,C1,C2+h[3],C3-h[4],C4,C5) - ll(C0,C1,C2-h[3],C3+h[4],C4,C5) + ll(C0,C1,C2-h[3],C3-h[4],C4,C5) ) / (4*h[3]*h[4])
   IM[3,5] <- -( ll(C0,C1,C2+h[3],C3,C4+h[5],C5) - ll(C0,C1,C2+h[3],C3,C4-h[5],C5) - ll(C0,C1,C2-h[3],C3,C4+h[5],C5) + ll(C0,C1,C2-h[3],C3,C4-h[5],C5) ) / (4*h[3]*h[5])
   IM[3,6] <- -( ll(C0,C1,C2+h[3],C3,C4,C5+h[6]) - ll(C0,C1,C2+h[3],C3,C4,C5-h[6]) - ll(C0,C1,C2-h[3],C3,C4,C5+h[6]) + ll(C0,C1,C2-h[3],C3,C4,C5-h[6]) ) / (4*h[3]*h[6])
   IM[4,4] <- -( ll(C0,C1,C2,C3+h[4],C4,C5) - ll(C0,C1,C2,C3,C4,C5) + ll(C0,C1,C2,C3-h[4],C4,C5) ) / (h[4]^2)
   IM[4,5] <- -( ll(C0,C1,C2,C3+h[4],C4+h[5],C5) - ll(C0,C1,C2,C3+h[4],C4-h[5],C5) - ll(C0,C1,C2,C3-h[4],C4+h[5],C5) + ll(C0,C1,C2,C3-h[4],C4-h[5],C5) ) / (4*h[4]*h[5])
   IM[4,6] <- -( ll(C0,C1,C2,C3+h[4],C4,C5+h[6]) - ll(C0,C1,C2,C3+h[4],C4,C5-h[6]) - ll(C0,C1,C2,C3-h[4],C4,C5+h[6]) + ll(C0,C1,C2,C3-h[4],C4,C5-h[6]) ) / (4*h[4]*h[6])
   IM[5,5] <- -( ll(C0,C1,C2,C3,C4+h[5],C5) - ll(C0,C1,C2,C3,C4,C5) + ll(C0,C1,C2,C3,C4-h[5],C5) ) / (h[5]^2)
   IM[5,6] <- -( ll(C0,C1,C2,C3,C4+h[5],C5+h[6]) - ll(C0,C1,C2,C3,C4+h[5],C5-h[6]) - ll(C0,C1,C2,C3,C4-h[5],C5+h[6]) + ll(C0,C1,C2,C3,C4-h[5],C5-h[6]) ) / (4*h[5]*h[6])
   IM[6,6] <- -( ll(C0,C1,C2,C3,C4,C5+h[6]) - ll(C0,C1,C2,C3,C4,C5) + ll(C0,C1,C2,C3,C4,C5-h[6]) ) / (h[6]^2)

   IM[2,1] <- IM[1,2]
   IM[3,1] <- IM[1,3]
   IM[3,2] <- IM[2,3]
   IM[4,1] <- IM[1,4]
   IM[4,2] <- IM[2,4]
   IM[4,3] <- IM[3,4]
   IM[5,1] <- IM[1,5]
   IM[5,2] <- IM[2,5]
   IM[5,3] <- IM[3,5]
   IM[5,4] <- IM[4,5]
   IM[6,1] <- IM[1,6]
   IM[6,2] <- IM[2,6]
   IM[6,3] <- IM[3,6]
   IM[6,4] <- IM[4,6]
   IM[6,5] <- IM[5,6]
   
   VarMLE <- solve(IM, tol=1e-030)    # Note: taking tol=1e-030. Not setting appropriate tolerance level may show the "system is computationally singular" error.

   hyperparam[subject,] <- c(C0,C1,C2,C3,C4,C5)
   hyperparamVar[subject,,] <- VarMLE
 }
   
 return(list(hyperparam=hyperparam, hyperparamVar=hyperparamVar))
  
}
}
 

