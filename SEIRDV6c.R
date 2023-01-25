  SEIRD5V <- function( S0, E0, I0, RE0, RI0, D0, V0,
                     alpha1, beta1, betaI1, gamma1, 
                     rho1, rho1I, zeta1, n1,
                     mchpt1 ){
  # Here alpha1 must be one dimension higher than chpt1
  X1a <- mchpt1$X1a
  X1b <- mchpt1$X1b
  X1g <- mchpt1$X1g
  X1r <- mchpt1$X1r
  
  X1rI <- mchpt1$X1rI
  
  beta1v <- X1b*betaI1
  alpha1v <- X1a%*%alpha1
  
  rho1v <- X1r[,1]*rho1[1] + X1r[,2]*rho1[2]
  rho1Iv <- X1rI*rho1I
  
  gamma1v <- X1g%*%gamma1
  Out1a <- SEIRDV5cpp(S0, E0, I0, RE0, RI0, D0, V0, 
                      alpha1v, beta1, beta1v, gamma1v, 
                      rho1v, rho1Iv, rho1I, zeta1= zeta1, n1= n1 )
  Out1 <- data.frame( S = Out1a$S,
                      E = Out1a$E,
                      I = Out1a$I,
                      RE = Out1a$RE,
                      RI = Out1a$RI,
                      D = Out1a$D,
                      V = Out1a$V,
                      Rep0 = Out1a$R0 )

  return( Out1 )
}


LikeSEIRD5V2 <- function( data1,  S0, E0, I0, RE0, RI0, D0, V0,
                          alpha1, beta1, betaI1, gamma1, 
                          rho1, rho1I, zeta1, n1,
                          mchpt1 ){
  # Call the solver
  Out1 <- SEIRD5V(  S0, E0, I0, RE0, RI0, D0, V0,
                    alpha1, beta1, betaI1, gamma1, 
                    rho1, rho1I, zeta1, n1,
                    mchpt1 )
  # Evaluate the likelihood at each point.
  Lik1AI <- sum( dpois( data1$AdjInfect, Out1$I, log = TRUE ) )
  #Lik1ER <- sum( dpois( data1$Recovered, Out1$R, log = TRUE ) )
  Lik1IR <- sum( dpois( data1$Recovered, Out1$RI, log = TRUE ) )
  Lik1D <- sum( dpois( data1$Deaths, Out1$D, log = TRUE ) )
  Lik1V <- sum( dpois( data1$Vaccinated, Out1$V, log = TRUE ) )
  res1 <- Lik1AI + Lik1IR + Lik1D + Lik1V
  return( res1 )
}



PostSEIRD5V1 <- function( data1, S0, E0, I0, RE0, RI0, D0, V0,
                          alpha1, beta1, betaI1, gamma1, 
                          rho1, rho1I, zeta1, n1,
                          mchpt1,
                          prior1am, prior1b, prior1g, prior1z ,
                          prior1r){
  alpha1F <- c( alpha1 )
  betatest1 <- cumsum( c(beta1, betaI1) )
  if( min( alpha1F)  > 0 & 
      min( gamma1) > 0  ){
    Post1 <- dexp( zeta1, prior1z, log = TRUE ) 
    Post1 <- Post1 + sum( dexp( rho1, prior1r, log = TRUE ) )
    Post1 <- Post1 + sum( dexp( rho1I, prior1r, log = TRUE ) )
    Post1 <- Post1 + sum( dexp( beta1, prior1b, log = TRUE ) )
    Post1 <- Post1 + sum( dexp( betaI1, prior1b, log = TRUE ) )
    Post1 <- Post1 + sum( dexp( gamma1, prior1g, log = TRUE ) )
    Post1 <- Post1 + sum( dexp( alpha1F, prior1am, log = TRUE) )
    Lik1 <- LikeSEIRD5V2( data1, S0, E0, I0, RE0, RI0, D0, V0,
                          alpha1, beta1, betaI1, gamma1, 
                          rho1, rho1I, zeta1, n1,
                          mchpt1 )
    res1 <- Post1 + Lik1
  }else{
    res1 <- -Inf
  }
  return( res1 )
}



# Build the matrix for change points.
MatrixBuild1 <- function( chpt1, n1 ){
  l1 <- length( chpt1 )
  res1 <- matrix( 0, nrow = n1, ncol = l1 )
  for( i in 1:(l1-1)){
    res1[(chpt1[i]:(chpt1[i+1]-1)),i] <- 1 
  }
  res1[ chpt1[l1]:n1, l1 ] <- 1
  return(res1)
}



MCMCSEIRD5VPred <- function( data1, S0, E0, I0, RE0, RI0, D0, V0, 
                             alpha1, beta1, betaI1, gamma1, zeta1, rho1, rho1I,
                             n1, chpt1, rchpt1, gchpt1, ImpI1, ImpRI1,
                             nsamp1,
                             alpha1Step,
                             beta1Step,
                             betaI1Step,
                             gamma1Step,
                             zeta1Step, 
                             rho1Step, 
                             rho1IStep,
                             prior1am, prior1as,
                             prior1b, prior1g, prior1z, prior1r){
  X1a <- MatrixBuild1( chpt1, n1 )
  X1g <- MatrixBuild1( gchpt1, n1 )
  X1r <- MatrixBuild1( rchpt1, n1 )
  X1b <- rep(0, n1)
  X1b[ImpI1] <- 1
  X1rI <- rep(0, n1)
  X1rI[ImpRI1 ] <- 1
  mchpt1 <- list( X1a = X1a, X1g = X1g, X1r = X1r, X1b = X1b, X1rI= X1rI )
  Post1 <- PostSEIRD5V1( data1,  S0, E0, I0, RE0, RI0, D0, V0,
                         alpha1, beta1, betaI1, gamma1, 
                         rho1, rho1I, zeta1, n1,
                         mchpt1,
                         prior1am, prior1b, prior1g, prior1z ,
                         prior1r)
  alpha1Out <- matrix( 0, nrow = nsamp1, ncol = length( alpha1 ) )
  beta1Out <- rep( 0, nsamp1)
  betaI1Out <- rep( 0, nsamp1)
  gamma1Out <- matrix( 0, nrow = nsamp1, ncol = length( gamma1 ) )
  zeta1Out <- rep( 0, nsamp1 )
  rho1Out <- matrix( 0, nrow = nsamp1, ncol = length( rho1 ) )
  rho1IOut <- rep( 0, nsamp1 )
  Post1Out <- rep(0, nsamp1 )
  IPred1 <- matrix( 0, nrow = nsamp1, ncol = n1 )  
  REPred1 <- matrix( 0, nrow = nsamp1, ncol = n1 )
  RIPred1 <- matrix( 0, nrow = nsamp1, ncol = n1 )
  DPred1 <- matrix( 0, nrow = nsamp1, ncol = n1 )
  VPred1 <- matrix( 0, nrow = nsamp1, ncol = n1 )
  for( i in 1:nsamp1 ){
    for( j in 1:length( alpha1 ) ){
      alpha1t <- alpha1
      alpha1t[j] <- alpha1[j] + alpha1Step[j]*rnorm( 1, 0, 1 )
      if( alpha1t[j] > 0 ){
        Post1t <- PostSEIRD5V1( data1, S0, E0, I0, RE0, RI0, D0, V0,
                                alpha1t, beta1, betaI1, gamma1, 
                                rho1, rho1I, zeta1, n1,
                                mchpt1, prior1am, prior1b, 
                                prior1g, prior1z,
                                prior1r )
        diff1 <- Post1t - Post1
        U1 <- log( runif( 1, 0, 1) )
        if( diff1 > U1 ){
          Post1 <- Post1t
          alpha1 <- alpha1t 
        }
      }
    }
    
    beta1t <- beta1 + beta1Step*rnorm( 1, 0, 1 )
    if( beta1t > 0 ){
      Post1t <- PostSEIRD5V1( data1, S0, E0, I0, RE0, RI0, D0, V0,
                              alpha1, beta1t, betaI1, gamma1, 
                              rho1, rho1I, zeta1, n1,
                              mchpt1,
                              prior1am, prior1b, prior1g, prior1z,
                              prior1r)
      diff1 <- Post1t - Post1
      U1 <- log( runif( 1, 0, 1) )
      if( diff1 > U1 ){
        Post1 <- Post1t
        beta1 <- beta1t 
      }
    }
    
    betaI1t <- betaI1 + betaI1Step*rnorm( 1, 0, 1 )
    if( betaI1t > 0  ){
      Post1t <- PostSEIRD5V1( data1, S0, E0, I0, RE0, RI0, D0, V0,
                              alpha1, beta1, betaI1t, gamma1, 
                              rho1, rho1I, zeta1, n1,
                              mchpt1,
                              prior1am, prior1b, prior1g, prior1z,
                              prior1r ) 
      diff1 <- Post1t - Post1
      U1 <- log( runif( 1, 0, 1) )
      if( diff1 > U1 ){
        Post1 <- Post1t
        betaI1 <- betaI1t 
      }
    }
    for( j in 1:length( gamma1 ) ){
      gamma1t <- gamma1
      gamma1t[j] <- gamma1[j] + gamma1Step[j]*rnorm( 1, 0, 1 ) 
      if( gamma1t[j] > 0) {
        Post1t <- PostSEIRD5V1( data1, S0, E0, I0, RE0, RI0, D0, V0,
                                alpha1, beta1, betaI1, gamma1t, 
                                rho1, rho1I, zeta1, n1,
                                mchpt1,
                                prior1am, prior1b, prior1g, prior1z,
                                prior1r)
        diff1 <- Post1t - Post1
        U1 <- log( runif( 1, 0, 1) )
        if( diff1 > U1 ){
          Post1 <- Post1t
          gamma1 <- gamma1t 
        }
      }
    }
    zeta1t <- zeta1 + zeta1Step*rnorm( 1, 0, 1 )
    if( zeta1t > 0 ){
      Post1t <- PostSEIRD5V1( data1, S0, E0, I0, RE0, RI0, D0, V0,
                              alpha1, beta1, betaI1, gamma1, 
                              rho1, rho1I, zeta1t, n1,
                              mchpt1,
                              prior1am, prior1b, prior1g, prior1z,
                              prior1r)
      diff1 <- Post1t - Post1
      U1 <- log( runif( 1, 0, 1) )
      if( diff1 > U1 ){
        Post1 <- Post1t
        zeta1 <- zeta1t 
      }
    }
    for(j in 2:length(rho1)){
      rho1t <- rho1
      rho1t[j] <- rho1[j] + rho1Step[j]*rnorm( 1, 0, 1 )
      if( rho1t[j] > 0 ){
        Post1t <- PostSEIRD5V1( data1, S0, E0, I0, RE0, RI0, D0, V0,
                                alpha1, beta1, betaI1, gamma1, 
                                rho1t, rho1I, zeta1, n1,
                                mchpt1,
                                prior1am, prior1b, prior1g, prior1z,
                                prior1r )
        diff1 <- Post1t - Post1
        U1 <- log( runif( 1, 0, 1) )
        if( diff1 > U1 ){
          Post1 <- Post1t
          rho1 <- rho1t 
        }
      }
    }
    
    rho1It <- rho1I
    rho1It <- rho1I + rho1IStep*rnorm( 1, 0, 1 )
    if( rho1It > 0 ){
      Post1t <- PostSEIRD5V1( data1, S0, E0, I0, RE0, RI0, D0, V0,
                              alpha1, beta1, betaI1, gamma1, 
                              rho1, rho1It, zeta1, n1,
                              mchpt1,
                              prior1am, prior1b, prior1g, prior1z,
                              prior1r)
      diff1 <- Post1t - Post1
      U1 <- log( runif( 1, 0, 1) )
      if( diff1 > U1 ){
        Post1 <- Post1t
        rho1I <- rho1It 
      }
    }
    
  
    Fit1Q <- SEIRD5V( S0, E0, I0, RE0, RI0, D0, V0,
                      alpha1, beta1, betaI1, gamma1, 
                      rho1, rho1I, zeta1, n1,
                      mchpt1 )
    IPred1[i,] <- ifelse( Fit1Q$I < 1000, rpois( n1, Fit1Q$I), round(rnorm(n1, Fit1Q$I, sqrt(Fit1Q$I)),0))
    REPred1[i,] <- ifelse( Fit1Q$RE < 1000, rpois( n1, Fit1Q$RE), round(rnorm(n1, Fit1Q$RE, sqrt(Fit1Q$RE)),0))
    RIPred1[i,] <- ifelse( Fit1Q$RI < 1000, rpois( n1, Fit1Q$RI), round(rnorm(n1, Fit1Q$RI, sqrt(Fit1Q$RI)),0))
    DPred1[i,] <- ifelse( Fit1Q$D < 1000, rpois( n1, Fit1Q$D), round(rnorm(n1, Fit1Q$D, sqrt(Fit1Q$D)),0))
    VPred1[i,] <- ifelse( Fit1Q$V < 1000, rpois( n1, Fit1Q$V), round(rnorm(n1, Fit1Q$V, sqrt(Fit1Q$V)),0))
    
    
    alpha1Out[i,] <- alpha1
    beta1Out[i] <- beta1
    betaI1Out[i] <- betaI1
    gamma1Out[i,] <- gamma1
    zeta1Out[i] <- zeta1
    rho1Out[i,] <- rho1
    rho1IOut[i] <- rho1I
    Post1Out[i] <- Post1
  }
  res1 <- list( alpha1 = alpha1Out,
                beta1 = beta1Out,
                betaI1 = betaI1Out,
                gamma1 = gamma1Out,
                zeta1 = zeta1Out,
                rho1 = rho1Out,
                rho1I = rho1IOut,
                Post1 = Post1Out,
                IPred1 = IPred1,
                REPred1 = REPred1,
                RIPred1 = RIPred1,
                DPred1 = DPred1,
                VPred1 = VPred1,
                Rep0 = Rep0)
  
  return( res1 )
}






