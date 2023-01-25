######################################################
#
# Lecture notes 4/8/2020

# Read in the Covid-19 data for 4/8/2020
# Under course documents in Blackboard
#s#etwd( "G:/My Drive/RGandED/Amona/RCode")
#Cdata1 <- read.csv( "GlobalCovid19-2020-05-01.csv",
#                    header = TRUE)
# setwd("/Volumes/GoogleDrive/My Drive/RGandED/Amona/RCode")

# Clear Everything
#rm( list= ls() )

Rep0=1 

library( Rcpp )
sourceCpp( "SEIRDV6.cpp" )

# Read in the functions we will use.
source("SEIRDV6c.R")


# Read in the data
Cdata1 <- read.csv( "QatarIRDV_ed.csv",header = TRUE)
dim(Cdata1)
#View(Cdata1)
# re-naming the dataset
Q1 <- Cdata1[39:550,]    #The first 38 days are zero so there is no need for them.
n1Q <- nrow(Q1)
# Adjust Infected to remove recovered and remove deaths
Q1$AdjInfect <- Q1$Confirmed - Q1$Recovered - Q1$Deaths

Cdata1 <- Q1
data1 <- Cdata1
n1 <- n1Q
#View(Cdata1)


###########################################################################
#
# Model for Qatar
#
###########################################################################
# Starting Values
# Data for Qatar
# S0 <- 2782000
# E0 <- 10
# I0 <- 1
# RE0 <- 0
# RI0 <- 0
# D0 <- 0
# V0 <- 0
# Rep0 <- 0

# Parameter value
# alpha0Q <-1.28e-08
# alpha1Q <- 2.74e-08
# alpha2Q <- 1.22e-08
# alpha3Q <- 6.11e-08
# alpha4Q <- 3.85e-08
# alpha5Q <- 4.01e-08
# alpha6Q <- 2.19e-08
# alpha7Q <- 4.00e-08
# alpha8Q <- 5.39e-08
# alpha9Q <- 2.74e-08
# alpha10Q <- 4.38e-08
# alpha11Q <- 5.49e-08
# alpha12Q <- 2.76e-08
# alpha13Q <- 3.61e-08
# alpha14Q <- 5.31e-08

# alpha0Q <- 1/27820000
# alpha1Q <- 1/27820000
# alpha2Q <- 1/27820000
# alpha3Q <- 1/27820000
# alpha4Q <- 1/27820000
# alpha5Q <- 1/27820000
# alpha6Q <- 1/27820000
# alpha7Q <- 1/27820000
# alpha8Q <- 1/27820000
# alpha9Q <- 1/27820000
# alpha10Q <- 1/27820000
# alpha11Q <- 1/27820000
# alpha12Q <- 1/27820000
# alpha13Q <- 1/27820000
# alpha14Q <- 1/27820000
# 
# betaI1 <- 0.000023
# beta1 <- 0.07640


# gamma0Q <- 2.88e-09
# gamma1Q <- 2.67e-09
# gamma2Q <- 2.73e-09
# gamma3Q <- 3.33e-09
# gamma4Q <- 4.73e-09
# gamma5Q <- 5.86e-09
# gamma6Q <- 7.10e-09
# 
# gamma0Q <- 1.88e-9/2
# gamma1Q <- 1.88e-9/2
# gamma2Q <- 1.00e-9/2
# gamma3Q <- 1.00e-9/2
# gamma4Q <- 1.00e-9/2
# gamma5Q <- 1.00e-9/2
# gamma6Q <- 1.00e-9/2

# gamma0Q <- 2.88e-12
# gamma1Q <- 4.67e-12
# gamma2Q <- 5.73e-12
# gamma3Q <- 7.33e-12
# gamma4Q <- 4.73e-12
# gamma5Q <- 8.86e-12
# gamma6Q <- 7.10e-12

# gamma0Q <- 1.88e-11
# gamma1Q <- 4.67e-11
# gamma2Q <- 5.73e-11
# gamma3Q <- 6.33e-11
# gamma4Q <- 4.73e-11
# gamma5Q <- 7.86e-11
# gamma6Q <- 6.10e-11

S0 <- 2782000
E0 <- 5
I0 <- 1
RE0 <- 0
RI0 <- 0
D0 <- 0
V0 <- 0

# Parameter value
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

betaI1 <- 0.23
beta1 <- 0.05340


gamma0Q <- 1.88e-6/2
gamma1Q <- 1.88e-6/2
gamma2Q <- 1.00e-5/2
gamma3Q <- 1.00e-6/2
gamma4Q <- 1.00e-6/2
gamma5Q <- 1.00e-5/2
gamma6Q <- 1.00e-5/2

# gamma0Q <- 4.71e-9
# gamma1Q <- 8.77e-8
# gamma2Q <- 6.95e-8
# gamma3Q <- 2.29e-8
# gamma4Q <- 5.14e-8
# gamma5Q <- 7.34e-9
# gamma6Q <- 9.34e-8

zeta1 <- 0.000013
rho0Q <- 0      # For when there is no vaccine
rho1Q <- 0.000391 ### assumed value for vaccine rate
#kappa0Q <- 1    # For when there is no vacccine
#kappa1Q <- -0.9 #0.9 ### assumed value


ImpI1 <- 12    # Time at which an impulse jump occurs.


#################################################################################
#
#
# Infection rate intervention points
chpt1 <- c(1, 12, 35, 48, 60, 71, 80, 87, 95, 115, 124, 136, 350, 355, 420)
alpha1 <- rep( 1/2782000,  length( chpt1 ) )

# Recovery rate changepoints
gchpt1 <- c(1, 60, 89, 130, 280, 350, 420 )
# gamma1 <- rep( 1.00e-3, length( gchpt1 ) )
gamma1 <- rep( 1.00e-6, length( gchpt1 ) )

# Vaccine intervention points.
rchpt1 <- c( 1, 426 )  # Time at which vaccine was deployed
rho1 <- c(0, 0.000391)

# When the vaccine is deployed
ImpRI1 <- 425
# ( Q1$Vaccinated[426])/( S0 - Q1$Confirmed[425] - Q1$Deaths[425] - Q1$Recovered[425] - Q1$Vaccinated[425])
#rho1I <- 0.21
rho1I <- 0.005
# Effacacy Change point
#kchpt1 <- c( 1, 426 )
#kappa1 <- c( 1, 0.99 )


gamma1Step <- rep(1.00e-9, length(gamma1) )
beta1Step <- 0.0001
zeta1Step <- 0.000001
rho1Step <- c(0,0.00001)
#kappa1Step <- c(0,0.0001)
alpha1Step <- rep( 1.00e-9, length( alpha1) ) # c(6.369427e-11, 8.000000e-11, 1.000000e-11, 5.000000e-11,5.000000e-11, 5.000000e-11, 5.000000e-11 )
betaI1Step <- 0.001
rho1IStep <- 1.00e-3


prior1am = 1
prior1b = 1
prior1g = 1
prior1z = 1
prior1r = 1

# As you add new intervention points....
#   1. Add an additional zero (initial value) to the parameter vector so that it will be recorded.
#   Example:  alpha1 <- c( alpha1, 5.902905e-08 )
#             gamma1 <- c( gamma1, 0 )
# 
#   2. Add a small value to the current step values so the model can use the 
#      the previous step values.
#   Example:  alpha1Step <- c(alpha1Step, 5.000000e-11)
#             gamma1Step <- c( gamma1Step, 0.0001 )


# Truncate things so we can get good estimates for the earlier parameters.
#alpha1 <- alpha1[1:5]
#chpt1 <- chpt1[1:6]

####### Read in the start step values from a file...
filename1 <- "StartStepValuesM1V5E3.RData"
load( filename1 )    # Everything will be named correctly since it was saved from R... below....
rho1IStep <- 1.00e-3
gamma1Step <- rep(1.00e-9, length(gamma1) )
# Comment the line above if you don't have starting start step values in a file.


tic <- Sys.time()
nMCMC1 <- 100
Q1MCMC <- MCMCSEIRD5VPred(data1, S0, E0, I0, RE0, RI0, D0, V0, 
                          alpha1, beta1, betaI1, gamma1, zeta1, rho1, rho1I,
                          n1, chpt1, rchpt1, gchpt1, ImpI1, ImpRI1,
                          nMCMC1,
                          alpha1Step,
                          beta1Step,
                          betaI1Step,
                          gamma1Step,
                          zeta1Step, 
                          rho1Step, 
                          rho1IStep,
                          prior1am, prior1as,
                          prior1b, prior1g, prior1z, prior1r)
Sys.time() - tic



#######################################################################
#
# Trace plots
#
#######################################################################
#pdf("alpha.pdf")
for( i in 1:length(alpha1) ){
  plot( Q1MCMC$alpha1[,i], type = "l" )
}
#dev.off()

#pdf("beta.pdf")
plot( Q1MCMC$beta1, type = "l" )
#dev.off()

#pdf("betaI.pdf")
plot( Q1MCMC$betaI1, type = "l" )
#dev.off()

#pdf("gamma.pdf")
for( i in 1:length( gamma1 ) ){
  plot( Q1MCMC$gamma1[,i], type = "l" )
}
#dev.off()

#pdf("gamma.pdf")
plot( Q1MCMC$zeta1, type = "l" )
#dev.off()


plot( Q1MCMC$rho1[,1], type = "l" )


#pdf("rho.pdf")
plot( Q1MCMC$rho1[,2], type = "l" )
#dev.off()

#pdf("gamma.pdf")
plot( Q1MCMC$rho1I, type = "l" )
#dev.off()

# Create files so you can read these in 
alpha1Step <- apply( Q1MCMC$alpha1, 2, sd )/4
gamma1Step <- apply( Q1MCMC$gamma1, 2, sd )/4
betaI1Step <- sd( Q1MCMC$betaI1)/4
beta1Step <- sd( Q1MCMC$beta1 )/4
zeta1Step <- sd( Q1MCMC$zeta1 )/4
rho1Step <- apply( Q1MCMC$rho1, 2, sd )/4
#kappa1Step <- apply( Q1MCMC$kappa1, 2, sd )/4
rho1IStep <- sd( Q1MCMC$rho1I)/4

# Step #So we can append other generated values to the above list 

# Create files so you can read these in
alpha1 <- Q1MCMC$alpha1[nMCMC1,]
beta1 <- Q1MCMC$beta1[nMCMC1]
betaI1 <- Q1MCMC$betaI1[nMCMC1]
gamma1 <- Q1MCMC$gamma1[nMCMC1,]
zeta1 <- Q1MCMC$zeta1[nMCMC1]
rho1 <- Q1MCMC$rho1[nMCMC1,]
#kappa1 <- Q1MCMC$kappa1[nMCMC1,]
rho1I <- Q1MCMC$rho1I[nMCMC1]


#########  Why not just use a save function?  This will create a R database and save it.  When you load it everything
#########  will be there.
save(alpha1, beta1, gamma1, zeta1, rho1, rho1I, n1Q, chpt1,
     rchpt1, gchpt1, betaI1, ImpI1, ImpRI1,
     alpha1Step, beta1Step, betaI1Step, gamma1Step, zeta1Step, rho1Step, rho1IStep, 
     file = filename1 )
########  You can add a load statement above and it will pull all these in so you don't need to start all over.

Q1MCMC20 <- Q1MCMC
save( Q1MCMC20, file = "Q1MCMC20.Rdat" )


# select a new starting point
##ind1MCMC <- sample(1:nMCMC1, 1, replace = FALSE )
#alpha1 <- Q1MCMC$alpha1[ind1MCMC,]
#beta1 <- Q1MCMC$beta1[ind1MCMC]
#betaI1 <- Q1MCMC$betaI1[ind1MCMC]
##gamma1 <- Q1MCMC$gamma1[ind1MCMC,]
#zeta1 <- Q1MCMC$eta1[ind1MCMC]
#rho1 <- Q1MCMC$rho1[ind1MCMC,]
#kappa1 <- Q1MCMC$kappa1[ind1MCMC,]
#rho1I <- Q1MCMC$rho1I[ind1MCMC]





#################################################################################

IPred1m <- apply( Q1MCMC$IPred1, 2, quantile, c(0.5,0.25,0.975) )
REPred1m <- apply( Q1MCMC$REPred1, 2, quantile, c(0.5,0.25,0.975) )
RIPred1m <- apply( Q1MCMC$RIPred1, 2, quantile, c(0.5,0.25,0.975) )
DPred1m <- apply( Q1MCMC$DPred1, 2, quantile, c(0.5,0.25,0.975) )
VPred1m <- apply( Q1MCMC$VPred1, 2, quantile, c(0.5,0.25,0.975) )


#plot( 1:n1Q, Q1$Vaccinated, col = "purple", type = "b")
#lines( 1:n1Q, Q1$AdjInfect, col = "orange", type = "b")
#lines( 1:n1Q, Q1$Recovered, col = "blue", type = "b" )
#lines( 1:n1Q, Q1$Deaths, col = "black", type = "b")
#lines( 1:n1Q, IPred1m[1,], col ="orange" )
#lines( 1:n1Q, RPred1m[1,], col ="blue")
#lines( 1:n1Q, DPred1m[1,], col ="black")
#lines( 1:n1Q, VPred1m[1,], col ="purple")
#pdf("Infected.pdf")
n1Q <- n1
plot( 1:n1Q, Q1$AdjInfect, col = "red", type = "b", lwd=2, ylab = " Actively Infected",
      xlab = "Days since 29Feb2020", main="Actively Infected")
#dev.off()
lines( 1:n1Q, IPred1m[1,], col ="red" )
lines( 1:n1Q, IPred1m[2,], col = "red", lty = 2 )
lines( 1:n1Q, IPred1m[3,], col = "red", lty = 2 )
#dev.off()
# IScale1 <- 450
# plot( 1:IScale1, Q1$AdjInfect[1:IScale1], col = "orange", type = "b")
# lines( 1:IScale1, IPred1m[1,1:IScale1], col ="orange" )
# lines( 1:IScale1, IPred1m[2,1:IScale1], col = "orange", lty = 2 )
# lines( 1:IScale1, IPred1m[3,1:IScale1], col = "orange", lty = 2 )

# plot( 1:n1Q, Q1$Recovered, col = "blue", type = "b")
# lines( 1:n1Q, REPred1m[1,], col = "blue")
# lines( 1:n1Q, REPred1m[2,], col = "blue", lty = 2)
# lines( 1:n1Q, REPred1m[3,], col = "blue", lty = 2)

#pdf("Recovered.pdf")
plot( 1:n1Q, Q1$Recovered, col = "seagreen", type = "b", lwd=2, ylab = " Recovered",
      xlab = "Days since 29Feb2020", main="Recovered")
#dev.off()
lines( 1:n1Q, RIPred1m[1,], col = "seagreen")
lines( 1:n1Q, RIPred1m[2,], col = "seagreen", lty = 2)
lines( 1:n1Q, RIPred1m[3,], col = "seagreen", lty = 2)
#dev.off()


#pdf("Vaccinated.pdf")
plot( 1:n1Q, Q1$Vaccinated, col = "darkblue", type = "b", lwd=3,
      xlab = "Days since 29April2021",ylab = "Vaccinated",main = "Vaccinated")
#dev.off()
lines( 1:n1Q, VPred1m[1,], col = "purple")
lines( 1:n1Q, VPred1m[2,], col = "purple", lty = 2)
lines( 1:n1Q, VPred1m[3,], col = "purple", lty = 2)
#dev.off()

# IScale1 <- 150
# plot( 1:IScale1, Q1$AdjInfect[1:IScale1], col = "orange", type = "b")
# lines( 1:IScale1, IPred1m[1,1:IScale1], col ="orange" )
# lines( 1:IScale1, IPred1m[2,1:IScale1], col = "orange", lty = 2 )
# lines( 1:IScale1, IPred1m[3,1:IScale1], col = "orange", lty = 2 )
# #abline( v = chpt1 )

#pdf("Deaths.pdf")
plot( 1:n1Q, Q1$Deaths, col = "black", type = "b", 
      lwd=2, xlab = "Days since 29Feb2020",ylab = "Deaths",main = "Deaths")
#dev.off()
lines( 1:n1Q, DPred1m[1,], col = "black")
lines( 1:n1Q, DPred1m[2,], col = "black", lty = 2)
lines( 1:n1Q, DPred1m[3,], col = "black", lty = 2)

lines( 1:n1Q, DPred1m[1,], col = "black")
lines( 1:n1Q, DPred1m[2,], col = "black", lty = 2)
lines( 1:n1Q, DPred1m[3,], col = "black", lty = 2)
#dev.off()


# Basic Model Fit Statistics
AIError <- sum( (Q1$AdjInfect - IPred1m[1,])^2 );AIError
RIError <- sum( (Q1$Recovered - RIPred1m[1,])^2 );RIError
DError <- sum( (Q1$Deaths - DPred1m[1,])^2 );DError
VError <- sum( (Q1$Vaccinated - VPred1m[1,])^2 );VError
TotError1 <- (n1Q-1)*(sd( Q1$AdjInfect)^2 + sd( Q1$Recovered )^2 + sd( Q1$Deaths )^2 + sd( Q1$Vaccinated )^2 )
TotError1
# Pseudo R2
PsuedoR2 <- 1 - (AIError + RIError + DError + VError)/TotError1

PsuedoR2


## The R2 looks great. Now, plotting the posterior distribution plots:


# plot( 1:n1Q, IPred1m[1,], ylim=c(min(IPred1m[2,]),max(IPred1m[3,])), type = "l", col = "darkgrey", lty = 3,
#       ylab = "Infections",
#       xlab = "Days since 29Feb2020", main="Infected-Post" )
# polygon( c(1:n1Q,n1Q:1), c( IPred1m[2,],rev( IPred1m[3,])),col = "darksalmon", border = "darksalmon")
# lines( 1:n1Q, IPred1m[1,], col = "red3", lwd = 2 )
# points( 1:n1Q, Q1$AdjInfect, col = "red3" )
# 
# 
# 
plot( 1:n1Q, IPred1m[3,], type = "l", col = "red3", lty = 3,
      ylab = "Infections",
      xlab = "Days since 29Feb2020", main = "Infected-Post")
polygon( x = c( 1:n1Q,rev(1:n1Q)), y =c( IPred1m[3,], rev( IPred1m[2,] ) ),
         col = "darksalmon", border = "darksalmon")
lines( 1:n1Q, IPred1m[1,], col = "red3", lwd = 2 )
points( 1:n1Q, Q1$AdjInfect, col = "red3" )
legend(x = "topright", legend = c("Infected", "Post"),
       lty = c(1, 1),           # Line types
       col = c("red3", "darksalmon"),           # Line colors
       lwd = 2)                 # Line width



plot( 1:n1Q, RIPred1m[3,], type = "l", col = "aquamarine4", lty = 3,
      ylab = "Recovered",
      xlab = "Days since 29Feb2020", main="Recovered-Post")
polygon( x = c( 1:n1Q,rev(1:n1Q)), y =c( RIPred1m[3,], rev( RIPred1m[2,] ) ),
         col = "aquamarine3", border = "aquamarine3") 
lines( 1:n1Q, RIPred1m[1,], col = "aquamarine4", lwd = 2 )
points( 1:n1Q, Q1$Recovered, col = "aquamarine4" )

legend(x = "topleft", legend = c("Recovered", "Post"),
       lty = c(1, 1),           # Line types
       col = c("aquamarine4", "aquamarine3"),           # Line colors
       lwd = 2)                 # Line width



#pdf("post-death.pdf")
plot( 1:n1Q, DPred1m[3,], type = "l", col = "black", lty = 3,
      ylab = "Deaths",
      xlab = "Days since 29Feb2020", main="Death-Post")
polygon( x = c( 1:n1Q,rev(1:n1Q)), y =c( DPred1m[3,], rev( DPred1m[2,] ) ),
         col = "azure4", border = "azure4") 
lines( 1:n1Q, DPred1m[1,], col = "black", lwd = 2 )
points( 1:n1Q, Q1$Deaths, col = "black" )
# legend(x = "topleft", legend = c("Deaths", "Post"),
#        lty = c(1, 1),           # Line types
#        col = c("black", "azure4"),           # Line colors
#        lwd = 2)                 # Line width
#dev.off()

#pdf("post-vaccine.pdf")
plot( 1:n1Q, VPred1m[3,], type = "l", col = "darkgoldenrod", lty = 3,
      ylab = "Vaccinated",
      xlab = "Days since 29April2021", main="Vaccinated-Post")
polygon( x = c( 1:n1Q,rev(1:n1Q)), y =c( VPred1m[3,], rev( VPred1m[2,] ) ),
         col = "chocolate", border = "chocolate") 
lines( 1:n1Q, VPred1m[1,], col = "#556B2F", lwd = 2 )
points( 1:n1Q, Q1$Vaccinated, col = "#556B2F" )
# legend(x = "topleft", legend = c("Deaths", "Post"),
#        lty = c(1, 1),           # Line types
#        col = c("#556B2F", "chocolate"),           # Line colors
#        lwd = 2)                 # Line width

#dev.off()

# plot( 1:n1Q, DPred1m[3,], type = "l", col = "darkgrey", lty = 3,
#       ylab = "Deaths",
#       xlab = "Days", main="Death-Post without vaccine" )
# polygon( x = c( 1:n1Q,rev(1:n1Q)), y =c( DPred1m[3,], rev( DPred1m[2,] ) ), col = "black", lwd = 2 )
# lines( 1:n1Q, Q1$Deaths, col = "black", lwd = 2 )


# plot( 426:631, DQuant1[1,], ylim=c(min(DQuant1[2,]),max(DQuant1[3,])), type = "l", col = "darkgrey", lty = 3,
#       ylab = "Deaths",
#       xlab = "Days", main="Death-Post without vaccine" )
# polygon( c(426:631,631:426), c(DQuant1[2,],rev(DQuant1[3,])),col = "azure4", border = "azure4")
# lines( 426:631, DQuant1[1,], col = "black", lwd = 2 )
# lines( 426:n1Q, DPred1m[1,][426:n1Q], col = "black", lwd = 2 )


# Reproduction number....

load( "Q1MCMC20.Rdat" )
ls()
n1 <- n1Q
length(Q1MCMC20 )

X1a <- MatrixBuild1( chpt1, n1 )
X1g <- MatrixBuild1( gchpt1, n1 )
X1r <- MatrixBuild1( rchpt1, n1 )
#X1k <- MatrixBuild1( kchpt1, n1 )
X1b <- rep(0, n1)
X1b[ImpI1] <- 1
X1rI <- rep(0, n1)
X1rI[ImpRI1 ] <- 1
mchpt1 <- list( X1a = X1a, X1g = X1g, X1r = X1r, X1b = X1b, X1rI= X1rI )


nMCMC1 <- nrow( Q1MCMC20$alpha1  )
Rep0 <- matrix( 0, nrow = nMCMC1, ncol = n1)
for (i in 1:nMCMC1) {
  
  SIROUT1 <- SEIRD5V( S0, E0, I0, RE0, RI0, D0, V0,
                      Q1MCMC20$alpha1[i,],Q1MCMC20$beta1[i],Q1MCMC20$betaI1[i],Q1MCMC20$gamma1[i,], 
                      Q1MCMC20$rho1[i,], Q1MCMC20$rho1I[i], Q1MCMC20$zeta1[i], n1, mchpt1)
  Rep0[i,] <- SIROUT1$Rep0
}
Rep0

Rep0m <- apply( Rep0, 2, quantile, c(0.5,0.025,0.975) )
plot( 2:n1, Rep0m[1,2:n1] , col = "palevioletred3", lwd=2, type = "l",
      xlab = "Days since 29Feb2020",
      ylab = "Effective Reproduction Number, R_e(t)",
      ylim = c(0,2.5))
# lines( 2:n1, Rep0m[2,2:n1], col = "Grey", lty = 2)
# lines( 2:n1, Rep0m[3,2:n1], col = "Grey", lty = 2)
polygon(c(rev(2:n1), 2:n1), c(rev(Rep0m[2,2:n1]), Rep0m[3,2:n1]),
        col = "aquamarine3", border = "aquamarine3")
lines( 2:n1, Rep0m[3,2:n1], col = "darkgreen", lty = 3)
abline( h =1.0, col = "black", lwd = 2, lty = 2)
abline(v=c(12,48,60,87,115, 420), col=c("steelblue", "steelblue","steelblue"), lty=c(3,3),lwd=c(3,3,3))

legend(x = "topright", legend = "R_e(t)",
       lty = c(1),           # Line types
       col = c("darkgreen"),           # Line colors
       lwd = 2)                 # Line width





#####################################
#####################################

# Rep0m <- apply( Rep0, 2, quantile, c(0.5,0.025,0.975) )
# plot( 2:n1, Rep0m[1,2:n1] , col = "palevioletred3", type = "l", lwd = 3,
#       xlab = "Days since 29Feb2020",
#       ylab = "Reproduction Number, R0(t)",
#       ylim = c(0,2.5))
# lines( 2:n1, Rep0m[2,2:n1], col = "palevioletred3", lty = 2)
# lines( 2:n1, Rep0m[3,2:n1], col = "palevioletred3", lty = 2)
# abline( h =1.0, col = "black", lwd = 2, lty = 2)
# abline(v=c(60,80,420), col=c("seagreen1", "seagreen1","seagreen1"), lty=c(3,3),lwd=c(3,3,3))


### Renaming the parameters:
alpha1 <- Q1MCMC20$alpha1
beta1 <- Q1MCMC20$beta1
betaI1 <- Q1MCMC20$betaI1
gamma1 <- Q1MCMC20$gamma1
rho1 <- Q1MCMC20$rho1
rho1I <- Q1MCMC20$rho1I
zeta1 <- Q1MCMC20$zeta1

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

mean(gamma1)
###  Mean median, sd, and quantile for rho1 #######
mean_rho1 <- apply(rho1, 2, mean);mean_rho1
sd_rho1 <- apply(rho1, 2, sd);sd_rho1
median_rho1 <- apply(rho1, 2, median);median_rho1
quantile_rho1 <- apply(rho1, 2, quantile, c(0.025,0.5,0.975));quantile_rho1


###  Mean median, sd, and quantile for rho1I ########
mean(rho1I)
sd(rho1I)
median(rho1I)
quantile(rho1I, c(0.025,0.5,0.975))

#### Now, to find the difference between the alpha's:

Conta1 <- alpha1[,2] - alpha1[,1]
Conta2 <- alpha1[,3] - alpha1[,2]
Conta3 <- alpha1[,4] - alpha1[,3]
Conta4 <- alpha1[,5] - alpha1[,4]
Conta5 <- alpha1[,6] - alpha1[,5]
Conta6 <- alpha1[,7] - alpha1[,6]
Conta7 <- alpha1[,8] - alpha1[,7]
Conta8 <- alpha1[,9] - alpha1[,8]
Conta9 <- alpha1[,10] - alpha1[,9]
Conta10 <- alpha1[,11] - alpha1[,10]
Conta11 <- alpha1[,12] - alpha1[,11]
Conta12 <- alpha1[,13] - alpha1[,12]
Conta13 <- alpha1[,14] - alpha1[,13]
Conta14 <- alpha1[,15] - alpha1[,14]

### Positive Contrast Means
mean(Conta2) ### this has a positive length
mean(Conta4) ### this has a negative length
mean(Conta7) ### this has a positive length
mean(Conta10) ### this has a positive length
mean(Conta11) ### this has a positive length
mean(Conta13) ### this has a positive length
mean(Conta14) ### this has a positive length

## Negative Contrast means........
mean(Conta1) ### this has a negative length
mean(Conta3) ### this has a negative length
mean(Conta5) ### this has a positive length
mean(Conta8) ### this has a positive length
mean(Conta6) ### this has a negative length
mean(Conta9) ### this has a negative length
mean(Conta12) ### this has a negative length

### The  Median of the contrast
median(Conta1)
median(Conta2)
median(Conta3)
median(Conta4)
median(Conta5)
median(Conta6)
median(Conta7)
median(Conta8)
median(Conta9)
median(Conta10)
median(Conta11)
median(Conta12)
median(Conta13)
median(Conta14)

### The sd of the contrast
sd(Conta1)
sd(Conta2)
sd(Conta3)
sd(Conta4)
sd(Conta5)
sd(Conta6)
sd(Conta7)
sd(Conta8)
sd(Conta9)
sd(Conta10)
sd(Conta11)
sd(Conta12)
sd(Conta13)
sd(Conta14)


quantile(Conta1,c(0.025,0.5,0.975))
quantile(Conta2,c(0.025,0.5,0.975))
quantile(Conta3,c(0.025,0.5,0.975))
quantile(Conta4,c(0.025,0.5,0.975))
quantile(Conta5,c(0.025,0.5,0.975))
quantile(Conta6,c(0.025,0.5,0.975))
quantile(Conta7,c(0.025,0.5,0.975))
quantile(Conta8,c(0.025,0.5,0.975))
quantile(Conta9,c(0.025,0.5,0.975))
quantile(Conta10,c(0.025,0.5,0.975))
quantile(Conta11,c(0.025,0.5,0.975))
quantile(Conta12,c(0.025,0.5,0.975))
quantile(Conta13,c(0.025,0.5,0.975))
quantile(Conta14,c(0.025,0.5,0.975))

##### Transmission rates across each time interval, summing over each k

# Infection rate intervention points
chpt1 <- c(1, 12, 35, 40, 60, 70, 80, 87, 95, 105, 124, 136, 350, 355, 420)

alpha_k <- alpha1[,1:15]
mean(alpha_k)
sd(alpha_k)
quantile(alpha_k,c(0.025,0.5,0.975))
##############
#alpha_k <- alpha1[,]

#for (i in 1:dim(alpha1)[2]) {
# alpha_k[,i] <- alpha1[,i]
#mean_k[i] <- mean(alpha_k[,i])
#median_k[i] <- median(alpha_k[,i])
#sd_k[i] <- sd(alpha_k[,i])
#quant_k[i] <- quantile(alpha_k[,i],c(0.025,0.5,0.975))
#}
#############

####%%%%%%%%%% Now, to find the difference between the gammas: %%%%%%%%%%%%%%%%%%%%%

Cont1 <- gamma1[,2] - gamma1[,1]
Cont2 <- gamma1[,3] - gamma1[,2]
Cont3 <- gamma1[,4] - gamma1[,3]
Cont4 <- gamma1[,5] - gamma1[,4]
Cont5 <- gamma1[,6] - gamma1[,5]
Cont6 <- gamma1[,7] - gamma1[,6]

### Positive Contrast Means
mean(Cont1) ### this has a positive length
mean(Cont2) ### this has a positive length
mean(Cont6) ### this has a positive length


## Negative Contrast means........
mean(Cont3) ### this has a positive length
mean(Cont4) ### this has a negative length
mean(Cont5) ### this has a negative length

### The  Median of the contrast
median(Cont1)
median(Cont2)
median(Cont3)
median(Cont4)
median(Cont5)
median(Cont6)

### The sd of the contrast
sd(Cont1)
sd(Cont2)
sd(Cont3)
sd(Cont4)
sd(Cont5)
sd(Cont6)

## Quantile
quantile(Cont1,c(0.025,0.5,0.975))
quantile(Cont2,c(0.025,0.5,0.975))
quantile(Cont3,c(0.025,0.5,0.975))
quantile(Cont4,c(0.025,0.5,0.975))
quantile(Cont5,c(0.025,0.5,0.975))
quantile(Cont6,c(0.025,0.5,0.975))







# nsamp1 = nMCMC1
# ###%%%%%%%%%%%%%%%%%%%%% The Posterior predictive function %%%%%%%%%%%%%%%%%%%%%%%
# Q1PostPred <- PostPred4V1(Q1MCMC, S0, E0, I0, RE0, RI0, D0, V0,
#                           alpha1, beta1, gamma1, rho1, rho1I, zeta1, n1Q, chpt1,
#                           rchpt1, gchpt1)
# 
# 
# 
# Q1PostPred <- PostPred4V1(Q1MCMC, S0, E0, I0, RE0, RI0, D0, V0, 
#                           alpha1, beta1, gamma1, zeta21, rho1, rho1I, n1Q, chpt1,
#                           rchpt1, kchpt1, gchpt1)
# 
# 
# 
# 
# 
# 
# #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# # Basic Model Fit Statistics
# AIError <- sum( (Q1$AdjInfect - Q1PostPred$I[2,] )^2 )
# RError <- sum( (Q1$Recovered - Q1PostPred$R[2,] )^2 )
# DError <- sum( (Q1$Deaths - Q1PostPred$D[2,])^2 )
# VError <- sum( (Q1$Vaccinated - Q1PostPred$V[2,])^2 )
# TotError1 <- (n1Q-1)*(sd( Q1$AdjInfect)^2 + sd( Q1$Recovered )^2 + sd( Q1$Deaths )^2 + sd( Q1$Vaccinated )^2 )
# 


# # Pseudo R2
# PsuedoR2 <- 1 - (AIError + RError + DError + VError)/TotError1






##########################################################################
#########################################################################
#######################################################################


# 
# 
# 
# 
# #prior1am = alpha1, prior1as = rep(10,length(alpha1)), prior1b = 1, prior1g = 1, prior1e = 1,
# #               prior1r =1, prior1k = 1000 
# 
# 
# 
# ##### Considering a few section
# plot( 1:100, Q1$AdjInfect[1:100], col = "red",xlab = "Days since 22Jan2020",ylab = "Active Infections",main = "Active Infections")
# lines( 1:100, Fit1Q$I[1:100], col = "red" )
# 
# #### using the whole dataset from Jan 2020 until Oct 2021
# plot(1:n1Q, Q1$AdjInfect, col = "red",xlab = "Days since 22Jan2020",ylab = "Active Infections",main = "Active Infections")
# lines( 1:n1Q, Fit1Q$I, col = "red" )
# 
# 
# #png( paste( "PlotsProp/QRecovered0-",Sys.Date(),".png" , sep = "" ) )
# plot( 1:n1Q, Q1$Recovered, col = "seagreen",xlab = "Days since 22Jan2020",ylab = "Recovered",main = "Recovered")
# lines( 1:n1Q, Fit1Q$R, col = "seagreen" )
# 
# 
# #png( paste( "PlotsProp/QDeaths0-",Sys.Date(),".png" , sep = "" ) )
# plot( 1:n1Q, Q1$Deaths, col = "black",xlab = "Days since 22Jan2020",ylab = "Deaths",main = "Deaths")
# lines( 1:n1Q, Fit1Q$D, col = "black")#dev.off()
# 
# 
# ### Plotting the vaccine
# plot( 1:n1Q, Q1$Vaccinated , col = "magenta1",xlab = "Days since 29APR2021",ylab = "Vaccinated",main = "Vaccinated")
# lines( 1:n1Q, Fit1Q$V, col = "magenta1")#dev.off()
# #dev.off()
# 
# 
# ### Adding abline to the active infected
# plot( 1:n1Q, Q1$AdjInfect, col = "blue",xlab = "Days since 22Jan2020",ylab = "Active Infections",main = "Active Infections")
# lines( 1:n1Q, Fit1Q$I, col = "blue" )
# abline( v = chpt1Q, lty = 3 )
# abline( v = chpt2Q, lty = 3 )
# abline( v = chpt3Q, lty = 3 )
# lines( 1:n1Q, qpois(0.975, Fit1Q$I), lty = 3, col = "blue" )
# lines( 1:n1Q, qpois(0.025, Fit1Q$I), lty = 3, col = "blue" )
# 
# 
# ### Adding abline to the active vaccinated
# plot( 1:n1Q, Q1$Vaccinated, col = "magenta1",xlab = "Days since 22Jan2020",ylab = "Vaccinated",main = "Vaccinated")
# lines( 1:n1Q, Fit1Q$V, col = "magenta1" )
# abline( v = chpt1Q, lty = 3 )
# abline( v = chpt2Q, lty = 3 )
# abline( v = chpt3Q, lty = 3 )
# lines( 1:n1Q, qpois(0.975, Fit1Q$V), lty = 3, col = "magenta1" )
# lines( 1:n1Q, qpois(0.025, Fit1Q$V), lty = 3, col = "magenta1" )
# 
# 
# #### Plots against time
# # plot(Cdata1$Date, Q1$AdjInfect, col = "red",xlab = "Days since 22Jan2020",ylab = "Active Infections",main = "Active Infections")
# # plot(Cdata1$Date, Q1$Recovered, col = "red",xlab = "Days since 22Jan2020",ylab = "Active Infections",main = "Active Infections")
# # plot( Cdata1$Date, Q1$Deaths, col = "black",xlab = "Days since 22Jan2020",ylab = "Deaths",main = "Deaths")
# # plot( Cdata1$Date, Q1$Vaccinated , col = "magenta1",xlab = "Days since 29APR2021",ylab = "Vaccinated",main = "Vaccinated")
# 
# ############################################################################
# #
# # MCMC Everthing...
# ############################################################################
# gamma1Step <- 0.001
# beta1Step <- 0.001
# eta1Step <- 0.001
# rho1Step <- 0.001
# kappa1Step <- 0.0001
# alpha1Step <- c(6.369427e-11, 8.000000e-11, 1.000000e-11, 5.000000e-11)
# betaI1Step <- 0.001
# 
# 
# # Parameter value 
# alpha1Q <- c( 1/15500000, -1/15510000, 1/6500000, -1/12000000)#
# beta1Q <- 1/14
# betaI1 <- 1/1.85
# gamma1Q <- 1/125
# gamma1Q <- 1/14
# eta1Q <- 1/9500
# rho1Q <- 1/100  ### assumed value
# kappa1Q <- 1 - 1/1000 #0.9 ### assumed value 
# 
# 
# 
# chpt1Q <- 50   # Time at which we believe an intervention happened.
# chpt2Q <- 60   # Time at which we believe an intervention happened.
# chpt3Q <- 75   # Time at which we believe an intervention happened.
# 
# 
# nMCMC1 <- 1000
# Q1MCMC <- MCMCSEIRD3( Q1, S0Q, E0Q, I0Q, RE0, RI0Q, D0Q, V0Q, 
#                       alpha1 = alpha1Q, 
#                       beta1Q, gamma1Q, eta1Q, rho1Q, kappa1Q, n1 = n1Q, 
#                       chpt1 = c(chpt1Q, chpt2Q, chpt3Q),
#                       betaI1, ImpI1,
#                       nsamp1 = nMCMC1,
#                       alpha1Step = alpha1Step,
#                       beta1Step = beta1Step,
#                       betaI1Step = betaI1Step,
#                       gamma1Step = gamma1Step,
#                       eta1Step = eta1Step,
#                       rho1Step = rho1Step,
#                       kappa1Step = kappa1Step,
#                       prior1am = alpha1Q, 
#                       prior1as = min(abs(alpha1Q)),
#                       prior1b = 1, prior1g = 1, prior1e = 1, prior1r = 1, prior1k = 1000)
# 
# plot( Q1MCMC$alpha1[,1], type = "l" )
# plot( Q1MCMC$alpha1[,2], type = "l" )
# plot( Q1MCMC$alpha1[,3], type = "l" )
# plot( Q1MCMC$alpha1[,4], type = "l" )
# plot( Q1MCMC$beta1, type = "l" )
# plot( Q1MCMC$betaI1, type = "l" )
# plot( Q1MCMC$gamma1, type = "l" )
# plot( Q1MCMC$eta1, type = "l" )
# plot( Q1MCMC$rho1, type = "l" )
# plot( Q1MCMC$kappa1, type = "l" )
# 
# alpha1Q <- Q1MCMC$alpha1[nMCMC1,]
# beta1Q <- Q1MCMC$beta1[nMCMC1]
# betaI1 <- Q1MCMC$beta1[nMCMC1]
# gamma1Q <- Q1MCMC$gamma1[nMCMC1]
# eta1Q <- Q1MCMC$eta1[nMCMC1]
# rho1Q <- Q1MCMC$rho1[nMCMC1]
# eta1Q <- Q1MCMC$kappa1[nMCMC1]
# 
# alpha1Step <- apply( Q1MCMC$alpha1, 2, sd )/4
# gamma1Step <- sd( Q1MCMC$gamma1 )/4
# betaI1Step <- sd( Q1MCMC$betaI1)/4
# beta1Step <- sd( Q1MCMC$beta1 )/4
# eta1Step <- sd( Q1MCMC$eta1 )/4
# rho1Step <- sd( Q1MCMC$rho1 )/4
# kappa1Step <- sd( Q1MCMC$kappa1 )/4
# 
# #alpha1 = c(alpha0Q, alpha1Q, alpha2Q, alpha3Q )
# 
# #Q2MCMC <- list( alpha1 = matrix( rep( alpha1, 100), ncol = 4, byrow = TRUE ),
# #                beta1 = rep( beta1Q, 100 ),
# #                gamma1 = rep( gamma1Q, 100),
# #                eta1 = rep( eta1Q, 100),
# #                Like1 = 1000,
# #                Post1 = 1000)
# 
# # Get the posterior predictive distribution
# Q1PostPred <- PostPred3( Q1MCMC, S0Q, E0Q, I0Q, RE0, RI0Q, D0Q, V0Q, n1Q, chpt1 = c(chpt1Q, chpt2Q, chpt3Q), ImpI1 )
# 
# # Basic Model Fit Statistics
# AIError <- sum( (Q1$AdjInfect - Q1PostPred$I[2,] )^2 )
# RError <- sum( (Q1$Recovered - Q1PostPred$R[2,] )^2 )
# DError <- sum( (Q1$Deaths - Q1PostPred$D[2,])^2 )
# VError <- sum( (Q1$Vaccinated - Q1PostPred$V[2,])^2 )
# TotError1 <- (n1Q-1)*(sd( Q1$AdjInfect)^2 + sd( Q1$Recovered )^2 + sd( Q1$Deaths )^2 + sd( Q1$Vaccinated )^2 )
# 
# # Pseudo R2
# PsuedoR2 <- 1 - (AIError + RError + DError + VError)/TotError1
# 
# 
# #png(paste( "PlotsProp/QInfectPost1-",Sys.Date(),".png", sep = "" ) )
# plot( 1:n1Q, Q1PostPred$I[3,], type = "l", col = "red", lty = 3,
#       ylab = "Infections",
#       xlab = "Days since 22Jan2020")
# polygon( x = c( 1:n1Q,rev(1:n1Q)), y =c( Q1PostPred$I[3,], rev( Q1PostPred$I[1,] ) ),
#          col = "lightpink", border = "lightpink") 
# lines( 1:n1Q, Q1PostPred$I[2,], col = "red", lwd = 2 )
# points( 1:n1Q, Q1$AdjInfect, col = "red" )
# #dev.off()
# 
# png( paste( "PlotsProp/QRecoveredPost1-",Sys.Date(),".png" , sep = "" ) )
# plot( 1:n1Q, Q1PostPred$R[3,], type = "l", col = "seagreen", lty = 3,
#       ylab = "Recovered",
#       xlab = "Days since 22Jan2020")
# polygon( x = c( 1:n1Q,rev(1:n1Q)), y =c( Q1PostPred$R[3,], rev( Q1PostPred$R[1,] ) ),
#          col = "palegreen", border = "palegreen") 
# lines( 1:n1Q, Q1PostPred$R[2,], col = "seagreen", lwd = 2 )
# points( 1:n1Q, Q1$Recovered, col = "seagreen" )
# dev.off()
# 
# 
# png( paste( "PlotsProp/QDeathsPost1-",Sys.Date(),".png" , sep = "" ) )
# plot( 1:n1Q, Q1PostPred$D[3,], type = "l", col = "black", lty = 3,
#       ylab = "Deaths",
#       xlab = "Days since 22Jan2020")
# polygon( x = c( 1:n1Q,rev(1:n1Q)), y =c( Q1PostPred$D[3,], rev( Q1PostPred$D[1,] ) ),
#          col = "grey80", border = "grey80") 
# lines( 1:n1Q, Q1PostPred$D[2,], col = "black", lwd = 2 )
# points( 1:n1Q, Q1$Deaths, col = "black" )
# dev.off()
# 
# 
# # Get the posterior predictive distribution
# n1pred <- n1Q + 7
# Q2PostPred <- PostPred3( Q1MCMC, S0Q, E0Q, I0Q, RE0, RI0Q, D0Q, V0Q, n1Q = n1pred, chpt1 = c(chpt1Q, chpt2Q, chpt3Q), ImpI1 )
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
# dev.off()
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
# dev.off()
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
# dev.off()
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
# dev.off()
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
# dev.off()


