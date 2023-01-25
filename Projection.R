######################################################



library( Rcpp )
sourceCpp( "SEIRDV6cc.cpp" )

# Read in the functions we will use.
source("SEIRDV6ccc.R")


# Read in the data
Cdata1 <- read.csv( "QatarIRDV_ed.csv",header = TRUE)
View(Cdata1)
# re-naming the dataset
Q1 <- as.data.frame( Cdata1[39:463,-9])    #The first 38 days are zero so there is no need for them.

dim(Q1)
View(Q1)
n1Q <- nrow(Q1)
# Adjust Infected to remove recovered and remove deaths
Q1$AdjInfect <- Q1$Confirmed - Q1$Recovered - Q1$Deaths

#Cdata1 <- Q1
data1 <- Q1
n1 <- n1Q
View(data1)


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
#rho0Q <- 0      # For when there is no vaccine
#rho1Q <- 0.000391 ### assumed value for vaccine rate
#kappa0Q <- 1    # For when there is no vacccine
#kappa1Q <- -0.9 #0.9 ### assumed value


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

# Vaccine intervention points.
#rchpt1 <- c( 1, 426 )  # Time at which vaccine was deployed
#rho1 <- c(0, 0.000391)

# When the vaccine is deployed
#ImpRI1 <- 426
# ( Q1$Vaccinated[426])/( S0 - Q1$Confirmed[425] - Q1$Deaths[425] - Q1$Recovered[425] - Q1$Vaccinated[425])
#rho1I <- 0.21



gamma1Step <- rep(1.00e-15, length(gamma1) )
beta1Step <- 0.0001
zeta1Step <- 0.000001
#rho1Step <- c(0,0.00001)
alpha1Step <- rep( 1.00e-8, length( alpha1) ) # c(6.369427e-11, 8.000000e-11, 1.000000e-11, 5.000000e-11,5.000000e-11, 5.000000e-11, 5.000000e-11 )
betaI1Step <- 0.001
#rho1IStep <- 0.00001


prior1am = 1
prior1b = 1
prior1g = 1
prior1z = 1
#prior1r = 1


####### Read in the start step values from a file...
filename2 <- "StartStepValuesM1V5E3.RData"
load( filename2 )    # Everything will be named correctly since it was saved from R... below....
# Comment the line above if you don't have starting start step values in a file.


tic <- Sys.time()
nMCMC1 <- 10000
Q1MCMC <- MCMCSEIRD5VPred(data1, S0, E0, I0, RE0, RI0, D0, 
                          alpha1, beta1, betaI1, gamma1, zeta1,
                          n1, chpt1, gchpt1, ImpI1,
                          nMCMC1,
                          alpha1Step,
                          beta1Step,
                          betaI1Step,
                          gamma1Step,
                          zeta1Step, 
                          #rho1Step, 
                          #rho1IStep,
                          prior1am, prior1as,
                          prior1b, prior1g, prior1z)
Sys.time() - tic



#######################################################################
#
# Trace plots
#
#######################################################################
for( i in 1:length(alpha1) ){
  plot( Q1MCMC$alpha1[,i], type = "l" )
}
plot( Q1MCMC$beta1, type = "l" )
plot( Q1MCMC$betaI1, type = "l" )
for( i in 1:length( gamma1 ) ){
  plot( Q1MCMC$gamma1[,i], type = "l" )
}
plot( Q1MCMC$zeta1, type = "l" )
# plot( Q1MCMC$rho1[,1], type = "l" )
# plot( Q1MCMC$rho1[,2], type = "l" )
#plot( Q1MCMC$kappa1[,1], type = "l" )
#plot( Q1MCMC$kappa1[,2], type = "l" )
#plot( Q1MCMC$rho1I, type = "l" )
plot( Q1MCMC$zeta1, type = "l" )

# Create files so you can read these in 
alpha1Step <- apply( Q1MCMC$alpha1, 2, sd )/4
gamma1Step <- apply( Q1MCMC$gamma1, 2, sd )/4
betaI1Step <- sd( Q1MCMC$betaI1)/4
beta1Step <- sd( Q1MCMC$beta1 )/4
zeta1Step <- sd( Q1MCMC$zeta1 )/4
#rho1Step <- apply( Q1MCMC$rho1, 2, sd )/4
#kappa1Step <- apply( Q1MCMC$kappa1, 2, sd )/4
#rho1IStep <- sd( Q1MCMC$rho1I)/4

# Step #So we can append other generated values to the above list 

# Create files so you can read these in
alpha1 <- Q1MCMC$alpha1[nMCMC1,]
beta1 <- Q1MCMC$beta1[nMCMC1]
betaI1 <- Q1MCMC$betaI1[nMCMC1]
gamma1 <- Q1MCMC$gamma1[nMCMC1,]
zeta1 <- Q1MCMC$zeta1[nMCMC1]
#rho1 <- Q1MCMC$rho1[nMCMC1,]
#kappa1 <- Q1MCMC$kappa1[nMCMC1,]
#rho1I <- Q1MCMC$rho1I[nMCMC1]


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
load( "StartStepValuesM1V5E3.RData")

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
# X1r <- MatrixBuild1( rchpt1, n1 )
# X1k <- MatrixBuild1( kchpt1, n1 )
X1b <- rep(0, n1)
X1b[ImpI1] <- 1
#X1rI <- rep(0, n1)
#X1rI[ImpRI1 ] <- 1
mchpt1 <- list( X1a = X1a, X1g = X1g, X1b = X1b)


n1P <- rep( 1, length(426:631) )
Dout1 <- matrix(0, nrow = 10000, ncol = length( n1P ) )
for( i in 1:10000){
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


hist( Dout1[,206]-D425)

quantile( Dout1[,206]-D425, c(0.025, 0.5, 0.975) )


# 
# 
# #### Plotting the projection
# 
# 
#  plot( 426:631, DQuant1[1,], ylim = c(350-D425, 450-D425 ), type = "l", col = "darkgrey", lty = 3,
#       ylab = "Deaths",
#       xlab = "Days", main="Death-Post without vaccine")
# polygon( x = c(426:631,631:426), c(DQuant1[2,], sort(DQuant1[3,], decreasing=TRUE) ) ,
#          col = "azure4", border = "azure4") 
# lines( 426:631,  DQuant1[1,], col = "black", lwd = 2 )
# points( 426:631, Q1$Deaths[426:631], col = "black" )
# 
# 
# plot( 426:631, DQuant1[1,], ylim = c(350-D425, 450-D425 ), type = "l", col = "darkgrey", lty = 3,
#       ylab = "Deaths",
#       xlab = "Days", main="Death-Post without vaccine")
# polygon( x = c(426:631,631:426), c(DQuant1[2,], sort(DQuant1[3,], decreasing=TRUE) ) ,
#          col = "azure4", border = "azure4") 
# lines( 426:631,  DQuant1[1,], col = "black", lwd = 2 )
# points( 426:631, Q1$Deaths[426:631], col = "black" )
# 
# ##############################################################
# 
# #### Take the difference between before and after
# 
# #### Take the Pois(Mean)
# 
# x <- 426:631
# diff <- D425-DPred1
# lambda <- mean(diff)
# Distn <- dpois(x,           # X-axis values (x = 0, 1, 2, ...)
#       lambda,      # Mean number of events that occur on the interval
#       log = FALSE) # If TRUE, probabilities are given as log
# 
#  P.death <- cumsum(Distn)
# 
#  
#  #### ploting the distribution
#  
#  plot(dpois(x, lambda), type = "h",lwd = 2,
#       main = "Poisson probability function",
#       ylab = "P(X = x)", xlab = "Number of Deaths")
#  
#  

#################################################################################

#################################################################################

IPred1m <- apply( Q1MCMC$IPred1, 2, quantile, c(0.5,0.025,0.975) )
REPred1m <- apply( Q1MCMC$REPred1, 2, quantile, c(0.5,0.025,0.975) )
RIPred1m <- apply( Q1MCMC$RIPred1, 2, quantile, c(0.5,0.025,0.975) )
DPred1m <- apply( Q1MCMC$DPred1, 2, quantile, c(0.5,0.025,0.975) )
#VPred1m <- apply( Q1MCMC$VPred1, 2, quantile, c(0.5,0.25,0.975) )


n1Q <- n1
plot( 1:n1Q, Q1$AdjInfect, col = "orange", type = "b", main='Infected', xlab = "Days", ylab = "Infected")
lines( 1:n1Q, IPred1m[1,], col ="orange" )
lines( 1:n1Q, IPred1m[2,], col = "orange", lty = 2 )
lines( 1:n1Q, IPred1m[3,], col = "orange", lty = 2 )


plot(1:n1Q, Q1$Recovered, col = "blue", type = "b", main= "Recovered", xlab = "Days", ylab = "Recovered")
lines( 1:n1Q, RIPred1m[1,], col = "blue")
lines( 1:n1Q, RIPred1m[2,], col = "blue", lty = 2)
lines( 1:n1Q, RIPred1m[3,], col = "blue", lty = 2)

 

# plot( 1:n1Q, Q1$Vaccinated, col = "purple", type = "b", main='Vaccinated', xlab = "Days", ylab = "Vaccinated")
# lines( 1:n1Q, VPred1m[1,], col = "purple")
# lines( 1:n1Q, VPred1m[2,], col = "purple", lty = 2)
# lines( 1:n1Q, VPred1m[3,], col = "purple", lty = 2)


plot( 1:n1Q, Q1$Deaths, col = "black", type = "b", main='Death', xlab = "Days", ylab = "Death")
lines( 1:n1Q, DPred1m[1,], col = "black")
lines( 1:n1Q, DPred1m[2,], col = "black", lty = 2)
lines( 1:n1Q, DPred1m[3,], col = "black", lty = 2)



# Basic Model Fit Statistics
AIError <- sum( (Q1$AdjInfect - IPred1m[1,])^2 );AIError
RIError <- sum( (Q1$Recovered - RIPred1m[1,])^2 );RIError
DError <- sum( (Q1$Deaths - DPred1m[1,])^2 );DError
#VError <- sum( (Q1$Vaccinated - VPred1m[1,])^2 );VError
TotError1 <- (n1Q-1)*(sd( Q1$AdjInfect)^2 + sd( Q1$Recovered )^2 + sd( Q1$Deaths )^2 + sd( Q1$Vaccinated )^2 )
TotError1
# Pseudo R2
PsuedoR2 <- 1 - (AIError + RIError + DError)/TotError1

PsuedoR2


## The R2 looks great. Now, plotting the posterior distribution plots:

plot( 1:n1Q, IPred1m[3,], type = "l", col = "red3", lty = 3,
      ylab = "Infections",
      xlab = "Days since 29Feb2020", main = "Infected-Post")
polygon( x = c( 1:n1Q,rev(1:n1Q)), y =c( IPred1m[3,], rev( IPred1m[2,] ) ),
         col = "darksalmon", border = "darksalmon") 
lines( 1:n1Q, IPred1m[1,], col = "red3", lwd = 2 )
points( 1:n1Q, Q1$AdjInfect, col = "red3" )


plot( 1:n1Q, RIPred1m[3,], type = "l", col = "aquamarine4", lty = 3,
      ylab = "Recovered",
      xlab = "Days since 29Feb2020", main="Recovered-Post")
polygon( x = c( 1:n1Q,rev(1:n1Q)), y =c( RIPred1m[3,], rev( RIPred1m[2,] ) ),
         col = "aquamarine3", border = "aquamarine3") 
lines( 1:n1Q, RIPred1m[1,], col = "aquamarine4", lwd = 2 )
points( 1:n1Q, Q1$Recovered, col = "aquamarine4" )


plot( 421:n1Q, DPred1m[1,][421:512], ylim=c(min(DPred1m[2,][421:512]),max(DPred1m[3,][421:512])), type = "l", col = "darkgrey", lty = 3,
      ylab = "Deaths",
      xlab = "Days", main="Deaths-Post with Vaccine")
polygon( x = c( 421:n1Q,rev(421:n1Q)), y =c( DPred1m[2,][421:512], rev( DPred1m[3,][421:512] ) ),
         col = "azure4", border = "azure4") 
lines( 421:n1Q, DPred1m[1,][421:512], col = "black", lwd = 2 )
#points( 421:n1Q, Q1$Deaths[421:n1Q], col = "black" )


plot( 426:631, DQuant1[1,], ylim=c(min(DQuant1[2,]),max(DQuant1[3,])), type = "l", col = "darkgrey", lty = 3,
      ylab = "Deaths",
      xlab = "Days", main="Death-Post without vaccine" )
polygon( c(426:631,631:426), c(DQuant1[2,],rev(DQuant1[3,])),col = "azure4", border = "azure4")
lines( 426:631, DQuant1[1,], col = "black", lwd = 2 )
lines( 421:n1Q, DPred1m[1,][421:512], col = "black", lwd = 2 )

# DQuant1 <- apply( Dout1, 2, quantile, c(0.5,0.025, 0.975))
# plot( 426:631, DQuant1[1,], ylim=c(min(DQuant1[2,]),max(DQuant1[3,])), type = "l", col = "darkgrey", lty = 3,
#       ylab = "Deaths",
#       xlab = "Days", main="Death-Post without vaccine" )
# polygon( c(426:631,631:426), c(DQuant1[2,],rev(DQuant1[3,])),col = "azure4", border = "azure4")
# lines( 426:631, DQuant1[1,], col = "black", lwd = 2 )
# points( 426:631, Q1$Deaths[426:631], col = "black" )



# plot( 421:n1Q, DPred1m[1,][421:512], ylim=c(min(DPred1m[2,][421:512]),max(DPred1m[3,][421:512])), type = "l", col = "darkgrey", lty = 3,
#       ylab = "Deaths",
#       xlab = "Days", main="Deaths-Post with Vaccine")
# polygon( x = c( 421:n1Q,rev(421:n1Q)), y =c( DPred1m[2,][421:512], rev( DPred1m[3,][421:512] ) ),
#          col = "azure4", border = "azure4")
# lines( 421:n1Q, DPred1m[1,][421:512], col = "black", lwd = 2 )
# points( 421:n1Q, Q1$Deaths[421:n1Q], col = "black" )



# DQuant1 <- apply( Dout1, 2, quantile, c(0.5,0.025, 0.975))
# plot( 426:631, DQuant1[1,][426:631], type = "l", col = "darkgrey", lty = 3,
#       ylab = "Deaths",
#       xlab = "Days", main="Death-Post without vaccine" )
# polygon( c(426:631,631:426), c(DQuant1[2,],rev(DQuant1[3,])),col = "azure4", border = "azure4")
# lines( 426:631, DQuant1[1,], col = "black", lwd = 2 )

# 
# hist( Dout1[,206]-D425)
# 
# quantile( Dout1[,206]-D425, c(0.025, 0.5, 0.975) )

# plot( 426:631, DQuant1[1,], ylim=c(min(DQuant1[2,]),max(DQuant1[3,])), type = "l", col = "darkgrey", lty = 3,
#       #       ylab = "Deaths",
#       #       xlab = "Days", main="Death-Post without vaccine" )
#       # polygon( c(426:631,631:426), c(DQuant1[2,],rev(DQuant1[3,])),col = "azure4", border = "azure4")
#       # lines( 426:631, DQuant1[1,], col = "black", lwd = 2 )
#       

#################################################################
################################################################
# Reproduction number....

load( "Q1MCMC20_p.Rdat" )
ls()
n1 <- n1Q
length(Q1MCMC20_p )

X1a <- MatrixBuild1( chpt1, n1 )
X1g <- MatrixBuild1( gchpt1, n1 )
#X1r <- MatrixBuild1( rchpt1, n1 )
#X1k <- MatrixBuild1( kchpt1, n1 )
X1b <- rep(0, n1)
X1b[ImpI1] <- 1
#X1rI <- rep(0, n1)
#X1rI[ImpRI1 ] <- 1
mchpt1 <- list( X1a = X1a, X1g = X1g, X1b = X1b )



### Renaming the parameters:
alpha1 <- Q1MCMC20_p$alpha1
beta1 <- Q1MCMC20_p$beta1
betaI1 <- Q1MCMC20_p$betaI1
gamma1 <- Q1MCMC20_p$gamma1
rho1 <- Q1MCMC20_p$rho1
rho1I <- Q1MCMC20_p$rho1I
zeta1 <- Q1MCMC20_p$zeta1

##  Mean median, sd, and quantile for alpha1 ###########
mean_alpha1 <- apply(alpha1, 2, mean);mean_alpha1
sd_alpha1 <- apply(alpha1, 2, sd);sd_alpha1
median_alpha1 <- apply(alpha1, 2, median);median_alpha1
quantile_alpha1 <- apply(alpha1, 2, quantile, c(0.025,0.5,0.975));quantile_alpha1

###  Mean median, sd, and quantile for beta1 ######
mean(beta1)
sd(beta1)
median(beta1)
quantile(beta1, c(0.025,0.5,0.975))

###  Mean median, sd, and quantile for BetaI1 ######
mean(betaI1)
sd(betaI1)
median(betaI1)
quantile(betaI1, c(0.025,0.5,0.975))

mean(zeta1)
sd(zeta1)
median(zeta1)
quantile(zeta1, c(0.025,0.5,0.975))


###  Mean median, sd, and quantile for gamma1 ######
mean_gamma1 <- apply(gamma1, 2, mean);mean_gamma1
sd_gamma1 <- apply(gamma1, 2, sd);sd_gamma1
median_gamma1 <- apply(gamma1, 2, median);median_gamma1
quantile_gamma1 <- apply(gamma1, 2, quantile, c(0.025,0.5,0.975));quantile_gamma1




########################################################################
########################################################################









# Get the posterior predictive distribution
 #n1pred <- n1Q + 125


# PostPred <- MCMCSEIRD5VPred( Q1MCMC,S0, E0, I0, RE0, RI0, D0, 
#                              alpha1, beta1, betaI1, gamma1, 
#                              zeta1, n1=n1pred,
#                              prior1b, prior1g, prior1z)
# 
# 


# plot( 1:n1pred, DPred1m[3,], type = "l", col = "red", lty = 3,
#       ylab = "Death",
#       xlab = "Days since 22Jan2020")
# polygon( x = c( 1:n1pred,rev(1:n1pred)), y =c( DPred1m[3,], rev( DPred1m[1,] ) ),
#          col = "lightpink", border = "lightpink")
# lines( 1:n1pred, DPred1m[2,], col = "red", lwd = 2 )
# abline( v = (n1Q+.5), lty = 3)
#  




# png(paste( "PlotsProp/QInfectPred1-",Sys.Date(),".png", sep = "" ) )
# plot( 1:n1pred, Q2PostPred$I[3,], type = "l", col = "red", lty = 3,
#       ylab = "Infections",
#       xlab = "Days since 22Jan2020")
# polygon( x = c( 1:n1pred,rev(1:n1pred)), y =c( Q2PostPred$I[3,], rev( Q2PostPred$I[1,] ) ),
#          col = "lightpink", border = "lightpink")
# lines( 1:n1pred, Q2PostPred$I[2,], col = "red", lwd = 2 )
# points( 1:n1Q, Q1$AdjInfect, col = "red" )
# #points( (n1Qa-6):n1Qa, Q1p$AdjInfect, col = "red", pch = "+")
# abline( v = (n1Q+.5), lty = 3)
# 
# 
# png( paste( "PlotsProp/QRecoveredPred1-",Sys.Date(),".png" , sep = "" ) )
# plot( 1:n1pred, Q2PostPred$R[3,], type = "l", col = "seagreen", lty = 3,
#       ylab = "Recovered",
#       xlab = "Days since 22Jan2020")
# polygon( x = c( 1:n1pred,rev(1:n1pred)), y =c( Q2PostPred$R[3,], rev( Q2PostPred$R[1,] ) ),
#          col = "palegreen", border = "palegreen")
# lines( 1:n1pred, Q2PostPred$R[2,], col = "seagreen", lwd = 2 )
# points( 1:n1Q, Q1$Recovered, col = "seagreen" )
# points( (n1Qa-6):n1Qa, Q1p$Recovered, col = "seagreen", pch = "+")
# abline( v = (n1Q+.5), lty = 3)
# 
# 
# png( paste( "PlotsProp/QDeathsPred1-",Sys.Date(),".png" , sep = "" ) )
# plot( 1:n1pred, Q2PostPred$D[3,], type = "l", col = "black", lty = 3,
#       ylab = "Deaths",
#       xlab = "Days since 22Jan2020")
# polygon( x = c( 1:n1pred,rev(1:n1pred)), y =c( Q2PostPred$D[3,], rev( Q2PostPred$D[1,] ) ),
#          col = "grey80", border = "grey80")
# lines( 1:n1pred, Q2PostPred$D[2,], col = "black", lwd = 2 )
# points( 1:n1Q, Q1$Deaths, col = "black" )
# points( (n1Qa-6):n1Qa, Q1p$Deaths, col = "black", pch = "+")
# abline( v = (n1Q+.5), lty = 3)
# 
# 
# png( paste( "PlotsProp/QVacinatedPred1-",Sys.Date(),".png" , sep = "" ) )
# plot( 1:n1pred, Q2PostPred$V[3,], type = "l", col = "coral4", lty = 3,
#       ylab = "Vaccinated",
#       xlab = "Days since 29APR2021")
# polygon( x = c( 1:n1pred,rev(1:n1pred)), y =c( Q2PostPred$V[3,], rev( Q2PostPred$V[1,] ) ),
#          col = "darkcyan", border = "darkcyan")
# lines( 1:n1pred, Q2PostPred$V[2,], col = "coral4", lwd = 2 )
# points( 1:n1Q, Q1$Vaccinated, col = "coral4" )
# points( (n1Qa-6):n1Qa, Q1p$vaccinated, col = "coral4", pch = "+")
# abline( v = (n1Q+.5), lty = 3)
# 
# 
# 
# png( paste( "PlotsProp/QNewInfPred1-",Sys.Date(),".png" , sep = "" ) )
# plot( 2:(n1pred), Q2PostPred$IN[3,], type = "l", col = "red", lty = 3,
#       ylab = "New Infections",
#       xlab = "Days since 22Jan2020",
#       ylim = c(0,511))
# polygon( x = c( 2:(n1pred),rev(2:(n1pred))), y =c( Q2PostPred$IN[3,], rev( Q2PostPred$IN[1,] ) ),
#          col = "lightpink", border = "lightpink")
# lines( 2:(n1pred), Q2PostPred$IN[2,], col = "red", lwd = 2 )
# points( 2:(n1Q), diff(Q1$AdjInfect), col = "red" )
# points( (n1Qa-5):n1Qa, diff(Q1p$AdjInfect), col = "red", pch = "+")
# abline( v = (n1Q+.5), lty = 3)



