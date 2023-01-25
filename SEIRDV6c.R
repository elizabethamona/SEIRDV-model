  SEIRD5V <- function( S0, E0, I0, RE0, RI0, D0, V0,
                     alpha1, beta1, betaI1, gamma1, 
                     rho1, rho1I, zeta1, n1,
                     mchpt1 ){
  # Here alpha1 must be one dimension higher than chpt1
  X1a <- mchpt1$X1a
  X1b <- mchpt1$X1b
  X1g <- mchpt1$X1g
  X1r <- mchpt1$X1r
  #X1l <- mchpt1$X1l
  #X1z2 <- mchpt1$X1z2
  X1rI <- mchpt1$X1rI
  #betaI1v <- X1b%*%betaI1
  beta1v <- X1b*betaI1
  alpha1v <- X1a%*%alpha1
  #lambda1v <- X1l%*%lambda1
  rho1v <- X1r[,1]*rho1[1] + X1r[,2]*rho1[2]
  rho1Iv <- X1rI*rho1I
  #zeta21v <- X1z2%*%zeta21
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
    #Post1 <- Post1 + sum( dbeta( kappa1, prior1k, 1, log = TRUE ) )        # This needs fixed
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
 # X1k <- MatrixBuild1( kchpt1, n1 )
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
  #kappa1Out <- matrix( 0, nrow = nsamp1, ncol = length( kappa1 ) )
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
    
    #for( j in 2:length(kappa1)){
      #kappa1t <- kappa1
      #kappa1t[j] <- kappa1[j] + kappa1Step[j]*rnorm( 1, 0, 1 )
      #if( kappa1t[j] > 0 & kappa1t[j] <= 1 ){
       # Post1t <- PostSEIRD5V1( data1, S0, E0, I0, RE0, RI0, D0, V0,
        #                        alpha1, beta1, betaI1, gamma1, eta1, rho1, rho1I, kappa1t, n1, mchpt1,
        #                        prior1am, prior1b, prior1g, prior1z,
        #                        prior1r, prior1k ) 
        #diff1 <- Post1t - Post1
        #U1 <- log( runif( 1, 0, 1) )
        #if( diff1 > U1 ){
          #Post1 <- Post1t
        #%kappa1 <- kappa1t 
       # }
     # }
 # }
  
    Fit1Q <- SEIRD5V( S0, E0, I0, RE0, RI0, D0, V0,
                      alpha1, beta1, betaI1, gamma1, 
                      rho1, rho1I, zeta1, n1,
                      mchpt1 )
    IPred1[i,] <- ifelse( Fit1Q$I < 1000, rpois( n1, Fit1Q$I), round(rnorm(n1, Fit1Q$I, sqrt(Fit1Q$I)),0))
    REPred1[i,] <- ifelse( Fit1Q$RE < 1000, rpois( n1, Fit1Q$RE), round(rnorm(n1, Fit1Q$RE, sqrt(Fit1Q$RE)),0))
    RIPred1[i,] <- ifelse( Fit1Q$RI < 1000, rpois( n1, Fit1Q$RI), round(rnorm(n1, Fit1Q$RI, sqrt(Fit1Q$RI)),0))
    DPred1[i,] <- ifelse( Fit1Q$D < 1000, rpois( n1, Fit1Q$D), round(rnorm(n1, Fit1Q$D, sqrt(Fit1Q$D)),0))
    VPred1[i,] <- ifelse( Fit1Q$V < 1000, rpois( n1, Fit1Q$V), round(rnorm(n1, Fit1Q$V, sqrt(Fit1Q$V)),0))
    #Out1IN[i,] <- diff( Out1I[i,] )
    #Out1Max1[i] <- which( Fit1Q$I == max(Fit1Q$I) )
    
    
    alpha1Out[i,] <- alpha1
    beta1Out[i] <- beta1
    betaI1Out[i] <- betaI1
    gamma1Out[i,] <- gamma1
    zeta1Out[i] <- zeta1
    rho1Out[i,] <- rho1
    #kappa1Out[i,] <- kappa1
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
               # kappa1 = kappa1Out,
                Post1 = Post1Out,
                IPred1 = IPred1,
                REPred1 = RIPred1,
                RIPred1 = RIPred1,
                DPred1 = DPred1,
                VPred1 = VPred1,
                Rep0 = Rep0)
  
  return( res1 )
}





##########################################################################
######################################################################
########################################################################

# SEIRD4V <- function( S0, E0, I0, RE0, RI0, D0, V0,
#                      alpha1, beta1, gamma1, eta1, rho1, kappa1, n1, chpt1,
#                      rchpt1, kchpt1, gchpt1, betaI1, ImpI1 ){
#   # Here alpha1 must be one dimension higher than chpt1
#   Out1 <- data.frame( S = rep( S0, n1 ),
#                       E = rep( E0, n1 ),
#                       I = rep( I0, n1 ),
#                       R = rep( RE0, RI0, n1 ),
#                       D = rep( D0, n1 ),
#                       V = rep( V0, n1 ))
#   X1a <- rep(1,n1)
#   X1b <- rep(1,n1)
#   X1g <- rep(1,n1)
#   X1r <- rep(1,n1)
#   X1k <- rep(1,n1)
#   for( i in 2:length( alpha1 ) ){
#     X1a <- cbind( X1a, ifelse( (1:n1) > chpt1[ i-1 ], 1, 0 ) )
#   }
#   for( i in 1:length( betaI1 )){
#     X1b <- cbind( X1b, ifelse( (1:n1) == ImpI1[i] , 1, 0 ) )
#   }
#   for( i in 2:length( rho1 ) ){
#     X1r <- cbind( X1r, ifelse( (1:n1) > rchpt1[ i-1 ], 1, 0 ) )
#   }
#   for( i in 2:length( kappa1 ) ){
#     X1k <- cbind( X1k, ifelse( (1:n1) > kchpt1[ i-1 ], 1, 0 ) )
#   }  
#   for( i in 2:length( gamma1 ) ){
#     X1g <- cbind( X1g, ifelse( (1:n1) > gchpt1[ i-1 ], 1, 0 ) )
#   }  
#   beta1v <- X1b%*%c( beta1, betaI1 )
#   alpha1v <- X1a%*%alpha1
#   rho1v <- X1r%*%rho1
#   kappa1v <- X1k%*%kappa1  
#   gamma1v <- X1g%*%gamma1
#   for( i in 2:n1 ){
#     S0n <- S0
#     E0n <- E0
#     I0n <- I0
#     RE0n <- RE0
#     RI0n <- RI0
#     D0n <- D0
#     V0n <- V0
#     betaE0n <- min( E0n, beta1v[i]*E0n)
#     S0 <- max(0, S0n - alpha1v[i]*S0n*E0n -rho1v[i]*S0n )
#     E0 <- max(0, E0n + alpha1v[i]*S0n*E0n - beta1*E0n + (1-kappa1v[i])*alpha1v[i]*V0n*E0n-gamma1v[i]*E0n)
#     I0 <- max(0, I0n + beta1*E0n - gamma1v[i]*I0n - eta1*I0n )
#     RE0 <- max(0,RE0n + gamma1v[i]*I0n +gamma1v[i]*E0n )
#     RI0 <- max(0, RI0, RI0n + gamma1v[i]*I0n +gamma1v[i]*E0n )
#     D0 <- max(0, D0n + eta1*I0n )
#     V0 <- max(0, V0n + rho1v[i]*S0n - (1-kappa1v[i])*alpha1v[i]*V0n*E0n)
#     Out1[ i, ] <- c( S0, E0, I0, RE0, RI0, D0, V0 )
#   }
#   return( Out1 )
# }
# 
# 
# 
# 
# 
# # Function for SEIRD1  Base function.
# 
# SEIRD1 <- function( S0, E0, I0, RE0, RI0, D0, V0,
#                     alpha1, alpha2, beta1, gamma1, eta1, rho1, kappa1, n1, chpt1 = chpt1 ){
#   Out1 <- data.frame( S = rep( S0, n1 ),
#                       E = rep( E0, n1 ),
#                       I = rep( I0, n1 ),
#                       D = rep( D0, n1 ) )
#   for( i in 2:n1 ){
#     S0n <- S0
#     E0n <- E0
#     I0n <- I0
#     RE0, RI0n <- RE0, RI0
#     D0n <- D0
#     V0n <- V0
#     
#     Intervention1 <- ifelse( i < chpt1, 0, 1 )
#     S0 <- max(0, S0n - alpha1v[i]*S0n*E0n -rho1*S0n )
#     E0 <- max(0, E0n + alpha1v[i]*S0n*E0n - beta1*E0n + (1-kappa1)*alpha1v[i]*V0n*I0n-gamma1*E0n)
#     I0 <- max(0, I0n + beta1*E0n - gamma1*I0n - eta1*I0n )
#     RE0, RI0 <- max(0, RE0, RI0n + gamma1*I0n +gamma1*E0n )
#     D0 <- max(0, D0n + eta1*I0n )
#     V0 <- max(0, V0n + rho1*S0n -(1-kappa1)*alpha1v[i]*V0n*I0n)
#     Out1[ i, ] <- c( S0, E0, I0, RE0, RI0, D0, V0 )
#   }
#   return( Out1 )
# }
# 
# 
# SEIRD2 <- function( S0, E0, I0, RE0, RI0, D0,V0,
#                     alpha1, beta1, gamma1, eta1,rho1, kappa1, n1, chpt1  ){
#   # Here alpha1 must be one dimension higher than chpt1
#   Out1 <- data.frame( S = rep( S0, n1 ),
#                       E = rep( E0, n1 ),
#                       I = rep( I0, n1 ),
#                       R = rep( RE0, RI0, n1 ),
#                       D = rep( D0, n1 ),
#                       V = rep( V0, n1 ))
#   X1 <- rep(1,n1)
#   for( i in 2:length( alpha1 ) ){
#     X1 <- cbind( X1, ifelse( (1:n1) > chpt1[ i-1 ], 1, 0 ) )
#   }
#   alpha1v <- X1%*%alpha1
#   for( i in 2:n1 ){
#     S0n <- S0
#     E0n <- E0
#     I0n <- I0
#     RE0, RI0n <- RE0, RI0
#     D0n <- D0
#     V0n <- V0
#     
#     S0 <- max(0, S0n - alpha1v[i]*S0n*E0n -rho1*S0n )
#     E0 <- max(0, E0n + alpha1v[i]*S0n*E0n - beta1*E0n + (1-kappa1)*alpha1v[i]*V0n*I0n-gamma1*E0n)
#     I0 <- max(0, I0n + beta1*E0n - gamma1*I0n - eta1*I0n )
#     RE0, RI0 <- max(0, RE0, RI0n + gamma1*I0n +gamma1*E0n )
#     D0 <- max(0, D0n + eta1*I0n )
#     V0 <- max(0, V0n + rho1*S0n -(1-kappa1)*alpha1v[i]*V0n*I0n)
#     Out1[ i, ] <- c( S0, E0, I0, RE0, RI0, D0, V0 )
#   }
#   return( Out1 )
# }
# 
# LikeSEIRD2 <- function( data1, S0, E0, I0, RE0, RI0, D0, V0,
#                         alpha1, beta1, gamma1, eta1, rho1, kappa1, n1, chpt1 ){
#   # Call the solver
#   Out1 <- SEIRD2( S0, E0, I0, RE0, RI0, D0, V0,
#                   alpha1, beta1, gamma1, eta1, rho1, kappa1, n1, chpt1 )
#   # Evaluate the likelihood at each point.
#   Lik1AI <- sum( dpois( data1$AdjInfect, Out1$I, log = TRUE ) )
#   Lik1R <- sum( dpois( data1$Recovered, Out1$R, log = TRUE ) )
#   Lik1D <- sum( dpois( data1$Deaths, Out1$D, log = TRUE ) )
#   Lik1V <- sum( dpois( data1$Vaccinated, Out1$V, log = TRUE ) )
#   res1 <- Lik1AI + Lik1R + Lik1D + Lik1V
#   return( res1 )
# }
# 
# 
# 
# SEIRD3 <- function( S0, E0, I0, RE0, RI0, D0, V0,
#                     alpha1, beta1, gamma1, eta1, rho1, kappa1, n1, chpt1,
#                     betaI1, ImpI1 ){
#   # Here alpha1 must be one dimension higher than chpt1
#   Out1 <- data.frame( S = rep( S0, n1 ),
#                       E = rep( E0, n1 ),
#                       I = rep( I0, n1 ),
#                       R = rep( RE0, RI0, n1 ),
#                       D = rep( D0, n1 ),
#                       V = rep( V0, n1 ) )
#   X1a <- rep(1,n1)
#   X1b <- rep(1,n1)
#   for( i in 2:length( alpha1 ) ){
#     X1a <- cbind( X1a, ifelse( (1:n1) > chpt1[ i-1 ], 1, 0 ) )
#   }
#   for( i in 1:length( betaI1 )){
#     X1b <- cbind( X1b, ifelse( (1:n1) == ImpI1[i] , 1, 0 ) )
#   }
#   beta1v <- X1b%*%c( beta1, betaI1 )
#   alpha1v <- X1a%*%alpha1
#   for( i in 2:n1 ){
#     S0n <- S0
#     E0n <- E0
#     I0n <- I0
#     RE0, RI0n <- RE0, RI0
#     D0n <- D0
#     V0n <- V0
#     S0 <- max(0, S0n - alpha1v[i]*S0n*E0n -rho1*S0n )
#     E0 <- max(0, E0n + alpha1v[i]*S0n*E0n - beta1*E0n + (1-kappa1)*alpha1v[i]*V0n*I0n-gamma1*E0n)
#     I0 <- max(0, I0n + beta1*E0n - gamma1*I0n - eta1*I0n )
#     RE0, RI0 <- max(0, RE0, RI0n + gamma1*I0n +gamma1*E0n )
#     D0 <- max(0, D0n + eta1*I0n )
#     V0 <- max(0, V0n + rho1*S0n -(1-kappa1)*alpha1v[i]*V0n*I0n)
#     Out1[ i, ] <- c( S0, E0, I0, RE0, RI0, D0, V0 )
#   }
#   return( Out1 )
# }
# 
# 
# SEIRD4 <- function( S0, E0, I0, RE0, RI0, D0, V0,
#                     alpha1, beta1, gamma1, eta1, rho1, kappa1, n1, chpt1,
#                     betaI1, ImpI1 ){
#   # Here alpha1 must be one dimension higher than chpt1
#   Out1 <- data.frame( S = rep( S0, n1 ),
#                       E = rep( E0, n1 ),
#                       I = rep( I0, n1 ),
#                       R = rep( RE0, RI0, n1 ),
#                       D = rep( D0, n1 ),
#                       V = rep( V0, n1 ))
#   X1a <- rep(1,n1)
#   X1b <- rep(1,n1)
#   for( i in 2:length( alpha1 ) ){
#     X1a <- cbind( X1a, ifelse( (1:n1) > chpt1[ i-1 ], 1, 0 ) )
#   }
#   for( i in 1:length( betaI1 )){
#     X1b <- cbind( X1b, ifelse( (1:n1) == ImpI1[i] , 1, 0 ) )
#   }
#   beta1v <- X1b%*%c( beta1, betaI1 )
#   alpha1v <- X1a%*%alpha1
#   for( i in 2:n1 ){
#     S0n <- S0
#     E0n <- E0
#     I0n <- I0
#     RE0, RI0n <- RE0, RI0
#     D0n <- D0
#     V0n <- V0
#     betaE0n <- min( E0n, beta1v[i]*E0n)
#     S0 <- max(0, S0n - alpha1v[i]*S0n*E0n -rho1*S0n )
#     E0 <- max(0, E0n + alpha1v[i]*S0n*E0n - beta1*E0n + (1-kappa1)*alpha1v[i]*V0n*I0n-gamma1*E0n)
#     I0 <- max(0, I0n + beta1*E0n - gamma1*I0n - eta1*I0n )
#     RE0, RI0 <- max(0, RE0, RI0n + gamma1*I0n +gamma1*E0n )
#     D0 <- max(0, D0n + eta1*I0n )
#     V0 <- max(0, V0n + rho1*S0n - (1-kappa1)*alpha1v[i]*V0n*I0n)
#     Out1[ i, ] <- c( S0, E0, I0, RE0, RI0, D0, V0 )
#   }
#   return( Out1 )
# }
# 
# 
# SEIRD4V <- function( S0, E0, I0, RE0, RI0, D0, V0,
#                      alpha1, beta1, gamma1, eta1, rho1, kappa1, n1, chpt1,
#                      rchpt1, kchpt1, gchpt1, betaI1, ImpI1 ){
#   # Here alpha1 must be one dimension higher than chpt1
#   Out1 <- data.frame( S = rep( S0, n1 ),
#                       E = rep( E0, n1 ),
#                       I = rep( I0, n1 ),
#                       R = rep( RE0, RI0, n1 ),
#                       D = rep( D0, n1 ),
#                       V = rep( V0, n1 ))
#   X1a <- rep(1,n1)
#   X1b <- rep(1,n1)
#   X1g <- rep(1,n1)
#   X1r <- rep(1,n1)
#   X1k <- rep(1,n1)
#   for( i in 2:length( alpha1 ) ){
#     X1a <- cbind( X1a, ifelse( (1:n1) > chpt1[ i-1 ], 1, 0 ) )
#   }
#   for( i in 1:length( betaI1 )){
#     X1b <- cbind( X1b, ifelse( (1:n1) == ImpI1[i] , 1, 0 ) )
#   }
#   for( i in 2:length( rho1 ) ){
#     X1r <- cbind( X1r, ifelse( (1:n1) > rchpt1[ i-1 ], 1, 0 ) )
#   }
#   for( i in 2:length( kappa1 ) ){
#     X1k <- cbind( X1k, ifelse( (1:n1) > kchpt1[ i-1 ], 1, 0 ) )
#   }  
#   for( i in 2:length( gamma1 ) ){
#     X1g <- cbind( X1g, ifelse( (1:n1) > gchpt1[ i-1 ], 1, 0 ) )
#   }  
#   beta1v <- X1b%*%c( beta1, betaI1 )
#   alpha1v <- X1a%*%alpha1
#   rho1v <- X1r%*%rho1
#   kappa1v <- X1k%*%kappa1  
#   gamma1v <- X1g%*%gamma1
#   for( i in 2:n1 ){
#     S0n <- S0
#     E0n <- E0
#     I0n <- I0
#     RE0, RI0n <- RE0, RI0
#     D0n <- D0
#     V0n <- V0
#     betaE0n <- min( E0n, beta1v[i]*E0n)
#     S0 <- max(0, S0n - alpha1v[i]*S0n*E0n -rho1v[i]*S0n )
#     E0 <- max(0, E0n + alpha1v[i]*S0n*E0n - beta1*E0n + (1-kappa1v[i])*alpha1v[i]*V0n*E0n-gamma1v[i]*E0n)
#     I0 <- max(0, I0n + beta1*E0n - gamma1v[i]*I0n - eta1*I0n )
#     RE0, RI0 <- max(0, RE0, RI0n + gamma1v[i]*I0n +gamma1v[i]*E0n )
#     D0 <- max(0, D0n + eta1*I0n )
#     V0 <- max(0, V0n + rho1v[i]*S0n - (1-kappa1v[i])*alpha1v[i]*V0n*E0n)
#     Out1[ i, ] <- c( S0, E0, I0, RE0, RI0, D0, V0 )
#   }
#   return( Out1 )
# }
# 
# 
# LikeSEIRD4V2 <- function( data1, S0, E0, I0, RE0, RI0, D0, V0, 
#                           alpha1, beta1, gamma1, eta1, rho1, kappa1, n1, chpt1,
#                           rchpt1, kchpt1, gchpt1, betaI1, ImpI1  ){
#   # Call the solver
#   Out1 <- SEIRD4V( S0, E0, I0, RE0, RI0, D0, V0,
#                    alpha1, beta1, gamma1, eta1, rho1, kappa1, n1, chpt1,
#                    rchpt1, kchpt1, gchpt1, betaI1, ImpI1 )
#   # Evaluate the likelihood at each point.
#   Lik1AI <- sum( dpois( data1$AdjInfect, Out1$I, log = TRUE ) )
#   Lik1R <- sum( dpois( data1$Recovered, Out1$R, log = TRUE ) )
#   Lik1D <- sum( dpois( data1$Deaths, Out1$D, log = TRUE ) )
#   Lik1V <- sum( dpois( data1$Vaccinated, Out1$V, log = TRUE ) )
#   res1 <- Lik1AI + Lik1R + Lik1D + Lik1V
#   return( res1 )
# }
# 
# 
# PostSEIRD4V1 <- function( data1, S0, E0, I0, RE0, RI0, D0, V0,
#                           alpha1, beta1, gamma1, eta1, rho1, kappa1, n1, chpt1,
#                           rchpt1, kchpt1, gchpt1, betaI1, ImpI1, 
#                           prior1am, prior1as, prior1b = 1, prior1g = 1, prior1z = 1,
#                           prior1r =1, prior1k = 1000 ){
#   alpha1F <- c( alpha1 )
#   betatest1 <- cumsum( c(beta1, betaI1) )
#   if( min( cumsum( alpha1F) ) > 0 & 
#       min( betatest1) > 0 & 
#       max( betatest1 ) < 1 ){
#     Post1 <- sum( dexp( eta1, prior1z, log = TRUE ) )
#     Post1 <- Post1 + sum( dexp( rho1, prior1r, log = TRUE ) )
#     Post1 <- Post1 + sum( dbeta( cumsum(kappa1), prior1k, 1, log = TRUE ) )        # This needs fixed
#     Post1 <- Post1 + sum( dexp( beta1, prior1b, log = TRUE ) )
#     Post1 <- Post1 + sum( dexp( betaI1, prior1b, log = TRUE ) )
#     Post1 <- Post1 + sum( dexp( cumsum(gamma1), prior1g, log = TRUE ) )
#     Post1 <- Post1 + sum( dnorm( alpha1F, prior1am, prior1as, log = TRUE) )
#     Lik1 <- LikeSEIRD4V2( data1, S0, E0, I0, RE0, RI0, D0, V0, 
#                           alpha1, beta1, gamma1, eta1, rho1, kappa1, n1, chpt1,
#                           rchpt1, kchpt1, gchpt1, betaI1, ImpI1 )
#     res1 <- Post1 + Lik1
#   }else{
#     res1 <- -Inf
#   }
#   return( res1 )
# }
# 
# 
# MCMCSEIRD4V1 <- function( data1, S0, E0, I0, RE0, RI0, D0, V0, 
#                           alpha1, beta1, gamma1, eta1, rho1, kappa1, n1, chpt1,
#                           rchpt1, kchpt1, gchpt1, betaI1, ImpI1, 
#                           nsamp1,
#                           alpha1Step,
#                           beta1Step,
#                           betaI1Step,
#                           gamma1Step,
#                           eta1Step, 
#                           rho1Step, 
#                           kappa1Step, prior1am, prior1as,
#                           prior1b, prior1g, prior1z, prior1r, prior1k){
#   Post1 <- PostSEIRD4V1( data1, S0, E0, I0, RE0, RI0, D0, V0, 
#                          alpha1, beta1, gamma1, eta1, rho1, kappa1, n1, chpt1, 
#                          rchpt1, kchpt1, gchpt1, betaI1, ImpI1,
#                          prior1am, prior1as, prior1b, prior1g, prior1z, prior1r, prior1k)
#   alpha1Out <- matrix( 0, nrow = nsamp1, ncol = length( alpha1 ) )
#   beta1Out <- rep( 0, nsamp1)
#   betaI1Out <- rep( 0, nsamp1)
#   gamma1Out <- matrix( 0, nrow = nsamp1, ncol = length( gamma1 ) )
#   eta1Out <- rep( 0, nsamp1 )
#   rho1Out <- matrix( 0, nrow = nsamp1, ncol = length( rho1 ) )
#   kappa1Out <- matrix( 0, nrow = nsamp1, ncol = length( kappa1 ) )
#   Post1Out <- rep(0, nsamp1 )
#   for( i in 1:nsamp1 ){
#     for( j in 1:length( alpha1 ) ){
#       alpha1t <- alpha1
#       alpha1t[j] <- alpha1[j] + alpha1Step[j]*rnorm( 1, 0, 1 )
#       if( min(cumsum(alpha1t)) > 0 ){
#         Post1t <- PostSEIRD4V1( data1, S0, E0, I0, RE0, RI0, D0, V0,
#                                 alpha1t, beta1, gamma1, eta1, rho1, kappa1, n1, chpt1,
#                                 rchpt1, kchpt1, gchpt1, betaI1, ImpI1,
#                                 prior1am, prior1as, prior1b, prior1g, prior1z, prior1r, prior1k) 
#         diff1 <- Post1t - Post1
#         U1 <- log( runif( 1, 0, 1) )
#         if( diff1 > U1 ){
#           Post1 <- Post1t
#           alpha1 <- alpha1t 
#         }
#       }
#     }
#     beta1t <- beta1 + beta1Step*rnorm( 1, 0, 1 )
#     Post1t <- PostSEIRD4V1( data1, S0, E0, I0, RE0, RI0, D0, V0,
#                             alpha1, beta1t, gamma1, eta1, rho1, kappa1, n1, chpt1,
#                             rchpt1, kchpt1, gchpt1, betaI1, ImpI1,
#                             prior1am, prior1as, prior1b, prior1g, prior1z, prior1r, prior1k) 
#     diff1 <- Post1t - Post1
#     U1 <- log( runif( 1, 0, 1) )
#     if( diff1 > U1 ){
#       Post1 <- Post1t
#       beta1 <- beta1t 
#     }
#     betaI1t <- betaI1 + betaI1Step*rnorm( 1, 0, 1 )
#     Post1t <- PostSEIRD4V1( data1, S0, E0, I0, RE0, RI0, D0, V0,
#                             alpha1, beta1, gamma1, eta1, rho1, kappa1, n1, chpt1,
#                             rchpt1, kchpt1, gchpt1, betaI1, ImpI1,
#                             prior1am, prior1as, prior1b, prior1g, prior1z, prior1r, prior1k) 
#     diff1 <- Post1t - Post1
#     U1 <- log( runif( 1, 0, 1) )
#     if( diff1 > U1 ){
#       Post1 <- Post1t
#       betaI1 <- betaI1t 
#     }
#     for( j in 1:length( gamma1 ) ){
#       gamma1t <- gamma1 
#       gamma1t[j] <- gamma1[j] + gamma1Step[j]*rnorm( 1, 0, 1 )
#       Post1t <- PostSEIRD4V1( data1, S0, E0, I0, RE0, RI0, D0, V0,
#                               alpha1, beta1, gamma1t, eta1, rho1, kappa1, n1, chpt1,
#                               rchpt1, kchpt1, gchpt1, betaI1, ImpI1,
#                               prior1am, prior1as, prior1b, prior1g, prior1z, prior1r, prior1k) 
#       diff1 <- Post1t - Post1
#       U1 <- log( runif( 1, 0, 1) )
#       if( diff1 > U1 ){
#         Post1 <- Post1t
#         gamma1 <- gamma1t 
#       }
#     }
#     eta1t <- eta1 + eta1Step*rnorm( 1, 0, 1 )
#     Post1t <- PostSEIRD4V1( data1, S0, E0, I0, RE0, RI0, D0, V0,
#                             alpha1, beta1, gamma1, eta1t, rho1, kappa1, n1, chpt1,
#                             rchpt1, kchpt1, gchpt1, betaI1, ImpI1,
#                             prior1am, prior1as, prior1b, prior1g, prior1z, prior1r, prior1k) 
#     diff1 <- Post1t - Post1
#     U1 <- log( runif( 1, 0, 1) )
#     if( diff1 > U1 ){
#       Post1 <- Post1t
#       eta1 <- eta1t 
#     }
#     for(j in 2:length(rho1)){
#       rho1t <- rho1
#       rho1t[j] <- rho1[j] + rho1Step[j]*rnorm( 1, 0, 1 )
#       Post1t <- PostSEIRD4V1( data1, S0, E0, I0, RE0, RI0, D0, V0,
#                               alpha1, beta1, gamma1, eta1, rho1t, kappa1, n1, chpt1,
#                               rchpt1, kchpt1, gchpt1, betaI1, ImpI1,
#                               prior1am, prior1as, prior1b, prior1g, prior1z, prior1r, prior1k) 
#       diff1 <- Post1t - Post1
#       U1 <- log( runif( 1, 0, 1) )
#       if( diff1 > U1 ){
#         Post1 <- Post1t
#         rho1 <- rho1t 
#       }
#     }
#     
#     for( j in 2:length(kappa1)){
#       kappa1t <- kappa1
#       kappa1t[j] <- kappa1[j] + kappa1Step[j]*rnorm( 1, 0, 1 )
#       Post1t <- PostSEIRD4V1( data1, S0, E0, I0, RE0, RI0, D0, V0,
#                               alpha1, beta1, gamma1, eta1, rho1, kappa1t, n1, chpt1,
#                               rchpt1, kchpt1, gchpt1, betaI1, ImpI1,
#                               prior1am, prior1as, prior1b, prior1g, prior1z, prior1r, prior1k) 
#       diff1 <- Post1t - Post1
#       U1 <- log( runif( 1, 0, 1) )
#       if( diff1 > U1 ){
#         Post1 <- Post1t
#         kappa1 <- kappa1t 
#       }
#     }
#     alpha1Out[i,] <- alpha1
#     beta1Out[i] <- beta1
#     betaI1Out[i] <- betaI1
#     gamma1Out[i,] <- gamma1
#     eta1Out[i] <- eta1
#     rho1Out[i,] <- rho1
#     kappa1Out[i,] <- kappa1
#     Post1Out[i] <- Post1
#   }
#   res1 <- list( alpha1 = alpha1Out,
#                 beta1 = beta1Out,
#                 betaI1 = betaI1Out,
#                 gamma1 = gamma1Out,
#                 eta1 = eta1Out,
#                 rho1 = rho1Out,
#                 kappa1 = kappa1Out,
#                 Post1 = Post1Out)
#   return( res1 )
# }
# 
# 
# 
# 
# 
# PostPred4V1 <- function( MCMCPar1, S0, E0, I0, RE0, RI0, D0, V0, n1Q, chpt1, ImpI1,rchpt1, kchpt1, 
#                          prior1am, prior1as, prior1b, prior1g, prior1z, prior1r, prior1k ){
#   n1 <- nrow( MCMCPar1$alpha1 )
#   Out1I <- matrix( 0, nrow = n1, ncol = n1Q )
#   Out1R <- matrix( 0, nrow = n1, ncol = n1Q )
#   Out1D <- matrix( 0, nrow = n1, ncol = n1Q )
#   Out1V <- matrix( 0, nrow = n1, ncol = n1Q )
#   Out1IN <- matrix( 0, nrow = n1, ncol = (n1Q-1) )
#   Out1Max1 <- rep( 0, n1 )
#   for( i in 1:n1 ){
#     alpha1 <- MCMCPar1$alpha1[i,]
#     beta1 <- MCMCPar1$beta1[i]
#     betaI1 <- MCMCPar1$betaI1[i]
#     gamma1 <- MCMCPar1$gamma1[i]
#     eta1 <- MCMCPar1$eta1[i]
#     rho1 <- MCMCPar1$rho1[i,]
#     kappa1 <- MCMCPar1$kappa1[i,]
#     Fit1Q <- SEIRD4V( S0, E0, I0, RE0, RI0, D0, V0,
#                       alpha1 = alpha1, 
#                       beta1Q, gamma1Q, eta1Q, rho1, kappa1, n1Q, 
#                       chpt1,
#                       rchpt1,
#                       kchpt1,
#                       betaI1, ImpI1 )
#     Out1I[i,] <- ifelse( Fit1Q$I < 1000, rpois( n1Q, Fit1Q$I), round(rnorm(n1Q, Fit1Q$I, sqrt(Fit1Q$I)),0))
#     Out1R[i,] <- ifelse( Fit1Q$R < 1000, rpois( n1Q, Fit1Q$R), round(rnorm(n1Q, Fit1Q$R, sqrt(Fit1Q$R)),0))
#     Out1D[i,] <- ifelse( Fit1Q$D < 1000, rpois( n1Q, Fit1Q$D), round(rnorm(n1Q, Fit1Q$D, sqrt(Fit1Q$D)),0))
#     Out1V[i,] <- ifelse( Fit1Q$V < 1000, rpois( n1Q, Fit1Q$V), round(rnorm(n1Q, Fit1Q$V, sqrt(Fit1Q$V)),0))
#     Out1IN[i,] <- diff( Out1I[i,] )
#     Out1Max1[i] <- which( Fit1Q$I == max(Fit1Q$I) )
#   }
#   
#   Out2I <- apply( Out1I, 2, quantile, c(0.025,0.5,0.975) )
#   Out2R <- apply( Out1R, 2, quantile, c(0.025,0.5,0.975) )
#   Out2D <- apply( Out1D, 2, quantile, c(0.025,0.5,0.975) )
#   Out2V <- apply( Out1V, 2, quantile, c(0.025,0.5,0.975) )
#   Out2IN <- apply( Out1IN, 2, quantile, c(0.025,0.5,0.975) )
#   Out2Max <- quantile( Out1Max1, c(0.025,0.5,0.975) )
#   res1 <- list( I = Out2I,
#                 R = Out2R,
#                 D = Out2D,
#                 V = Out2V,
#                 IN = Out2IN,
#                 Max = Out2Max)
#   return( res1 )
# }
# 
# 

