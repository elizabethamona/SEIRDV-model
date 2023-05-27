######################################################



library( Rcpp )
sourceCpp( "SEIRDV6cc.cpp" )

# Read in the functions we will use.
source("SEIRDV6ccc.R")


# Read in the data
Cdata1 <- read.csv( "QatarIRDV_ed.csv",header = TRUE)
#View(Cdata1)
# re-naming the dataset
Q1 <- as.data.frame( Cdata1[39:463,-9])    #The first 38 days are zero so there is no need for them.

dim(Q1)
#View(Q1)
n1Q <- nrow(Q1)
# Adjust Infected to remove recovered and remove deaths
Q1$AdjInfect <- Q1$Confirmed - Q1$Recovered - Q1$Deaths

#Cdata1 <- Q1
data1 <- Q1
n1 <- n1Q
#View(data1)


###########################################################################
#
# Model for Qatar
#
###########################################################################
# Starting Values
# Data for Qatar
S0 <- 2782000
E0 <- 5
I0 <- 1
RE0 <- 0
RI0 <- 0
D0 <- 0
#V0 <- 0
Rep0 <- 0



alpha0Q <- 1/2782000
alpha1Q <- 1/2782000
alpha2Q <- 1/2782000
alpha3Q <- 1/2782000
alpha4Q <- 1/2782000
alpha5Q <- 1/2782000
alpha6Q <- 1/2782000
alpha7Q <- 1/2782000
alpha8Q <- 1/2782000
alpha9Q <- 1/2782000
alpha10Q <- 1/2782000
alpha11Q <- 1/2782000
alpha12Q <- 1/2782000
alpha13Q <- 1/2782000
alpha14Q <- 1/2782000


betaI1 <- 0.000023
beta1 <- 0.07640


gamma0Q <- 1.88e-12/2
gamma1Q <- 1.88e-12/2
gamma2Q <- 1.00e-12/2
gamma3Q <- 1.00e-12/2
gamma4Q <- 1.00e-12/2
gamma5Q <- 1.00e-12/2
gamma6Q <- 1.00e-12/2




zeta1 <- 0.000013

ImpI1 <- 12    # Time at which an impulse jump occurs.


#################################################################################
#
#
# Infection rate intervention points
chpt1 <- c(1, 12, 35, 48, 60, 71, 80, 87, 95, 115, 300)
alpha1 <- rep( 1/2782000,  length( chpt1 ) )

# Recovery rate changepoints
gchpt1 <- c(1, 60, 89, 130, 190, 300)
gamma1 <- rep( 1.00e-12, length( gchpt1 ) )



gamma1Step <- rep(1.00e-15, length(gamma1) )
beta1Step <- 0.0001
zeta1Step <- 0.000001
alpha1Step <- rep( 1.00e-8, length( alpha1) ) # c(6.369427e-11, 8.000000e-11, 1.000000e-11, 5.000000e-11,5.000000e-11, 5.000000e-11, 5.000000e-11 )
betaI1Step <- 0.001



prior1am = 1
prior1b = 1
prior1g = 1
prior1z = 1


####### Read in the start step values from a file...
filename2 <- "StartStepValuesM1V5E3.RData"
load( filename2 )    # Everything will be named correctly since it was saved from R... below....
# Comment the line above if you don't have starting start step values in a file.


tic <- Sys.time()
nMCMC1 <- 1000
Q1MCMC <- MCMCSEIRD5VPred(data1, S0, E0, I0, RE0, RI0, D0, 
                          alpha1, beta1, betaI1, gamma1, zeta1,
                          n1, chpt1, gchpt1, ImpI1,
                          nMCMC1,
                          alpha1Step,
                          beta1Step,
                          betaI1Step,
                          gamma1Step,
                          zeta1Step, 
                          prior1am, prior1as,
                          prior1b, prior1g, prior1z)
Sys.time() - tic



#######################################################################
#
# Trace plots
#
#######################################################################
# for( i in 1:length(alpha1) ){
#   plot( Q1MCMC$alpha1[,i], type = "l" )
# }
# plot( Q1MCMC$beta1, type = "l" )
# plot( Q1MCMC$betaI1, type = "l" )
# for( i in 1:length( gamma1 ) ){
#   plot( Q1MCMC$gamma1[,i], type = "l" )
# }
# plot( Q1MCMC$zeta1, type = "l" )
# plot( Q1MCMC$zeta1, type = "l" )

# Create files so you can read these in 
alpha1Step <- apply( Q1MCMC$alpha1, 2, sd )/4
gamma1Step <- apply( Q1MCMC$gamma1, 2, sd )/4
betaI1Step <- sd( Q1MCMC$betaI1)/4
beta1Step <- sd( Q1MCMC$beta1 )/4
zeta1Step <- sd( Q1MCMC$zeta1 )/4


# Step #So we can append other generated values to the above list 

# Create files so you can read these in
alpha1 <- Q1MCMC$alpha1[nMCMC1,]
beta1 <- Q1MCMC$beta1[nMCMC1]
betaI1 <- Q1MCMC$betaI1[nMCMC1]
gamma1 <- Q1MCMC$gamma1[nMCMC1,]
zeta1 <- Q1MCMC$zeta1[nMCMC1]



#########  Why not just use a save function?  This will create a R database and save it.  When you load it everything
#########  will be there.
save(alpha1, beta1, gamma1, zeta1, n1Q, chpt1,
     gchpt1, betaI1, ImpI1,
     alpha1Step, beta1Step, betaI1Step, gamma1Step, zeta1Step, 
     file = filename2 )
########  You can add a load statement above and it will pull all these in so you don't need to start all over.

Q1MCMC20_p <- Q1MCMC
save( Q1MCMC20_p, file = "Q1MCMC20_p.Rdat" )

#############################################################################

# Load the data
library( "Rcpp" )
sourceCpp("SEIRDV6cc.cpp")
source( "SEIRDV6ccc.R")
load( "Q1MCMC20_p.Rdat" )
load( "StartStepValuesM1V5E2.RData")

S0 <- 2782000
E0 <- 5
I0 <- 1
RE0 <- 0
RI0 <- 0
D0 <- 0
#V0 <- 0
Rep0 <- 0

n1 <- 425
# Deaths at time 425
D425 <- Cdata1[n1,]$Deaths

X1a <- MatrixBuild1( chpt1, n1 )
X1g <- MatrixBuild1( gchpt1, n1 )
X1b <- rep(0, n1)
X1b[ImpI1] <- 1
mchpt1 <- list( X1a = X1a, X1g = X1g, X1b = X1b)


n1P <- rep( 1, length(426:631) )
Dout1 <- matrix(0, nrow = 1000, ncol = length( n1P ) )
for( i in 1:1000){
  alpha1t <- Q1MCMC20_p$alpha1[i,]
  beta1t <- Q1MCMC20_p$beta1[i]
  betaI1t <- Q1MCMC20_p$betaI1[i]
  gamma1t <- Q1MCMC20_p$gamma1[i,]
  zeta1t <- Q1MCMC20_p$zeta1[i]

  Fit1Q <- SEIRD5V( S0, E0, I0, RE0, RI0, D0,
                    alpha1t, beta1t, betaI1t, gamma1t,
                    zeta1t, n1,
                    mchpt1 )

  Pred1Q <- SEIRDV5cpp(Fit1Q$S[n1],
                       Fit1Q$E[n1],
                       Fit1Q$I[n1],
                       Fit1Q$RE[n1],
                       Fit1Q$RI[n1],
                       Fit1Q$D[n1],
                       alpha1[11]*n1P,
                       beta1,
                       beta1v = beta1*n1P,
                       gamma1[7]*n1P,
                       zeta1 = zeta1t,
                       n1 = length(n1P))
  IPred1 <- ifelse( Pred1Q$I < 1000, rpois( n1, Pred1Q$I), round(rnorm(n1, Pred1Q$I, sqrt(Pred1Q$I)),0))
  REPred1 <- ifelse( Pred1Q$RE < 1000, rpois( n1, Pred1Q$RE), round(rnorm(n1, Pred1Q$RE, sqrt(Pred1Q$RE)),0))
  RIPred1 <- ifelse( Pred1Q$RI < 1000, rpois( n1, Pred1Q$RI), round(rnorm(n1, Pred1Q$RI, sqrt(Pred1Q$RI)),0))
  Ddiff1 <- diff(Pred1Q$D)
  DPois1 <- rpois( length(Ddiff1), Ddiff1)
  DPred1 <- D425 + c(0,cumsum(DPois1))
  Dout1[i,] <- DPred1
}

DQuant1 <- apply( Dout1, 2, quantile, c(0.5,0.025, 0.975))
plot( 426:631, DQuant1[1,], ylim=c(min(DQuant1[2,]),max(DQuant1[3,])), type = "l", col = "darkgrey", lty = 3,
      ylab = "Deaths",
      xlab = "Days", main="Death-Post without vaccine" )
polygon( c(426:631,631:426), c(DQuant1[2,],rev(DQuant1[3,])),col = "azure4", border = "azure4")
lines( 426:631, DQuant1[1,], col = "black", lwd = 2 )
# points( 426:631, Q1$Deaths[426:631], col = "black" )
# lines( 426:631, DQuant1[1,], type = "l", col = "black")
# lines( 426:631, DQuant1[2,], type = "l", col = "black", lty = 3)
# lines( 426:631, DQuant1[3,], type = "l", col = "black", lty = 3)




