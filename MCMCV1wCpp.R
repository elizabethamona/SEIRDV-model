######################################################
# This code estimate the model parameters using Metropolis Hasting (MH) and solves the ODE using Euler method

# Clear Everything
#rm( list= ls() )

#setwd( "G:/My Drive/RGandED/Amona/RCode")
##################################################################################################################
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
# The Model
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
V0 <- 0
Rep0 <- 1 

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

rho1I <- 0.005
#### Initial startsteps values
gamma1Step <- rep(1.00e-9, length(gamma1) )
beta1Step <- 0.0001
zeta1Step <- 0.000001
rho1Step <- c(0,0.00001)
alpha1Step <- rep( 1.00e-9, length( alpha1) ) # c(6.369427e-11, 8.000000e-11, 1.000000e-11, 5.000000e-11,5.000000e-11, 5.000000e-11, 5.000000e-11 )
betaI1Step <- 0.001
rho1IStep <- 1.00e-3


prior1am = 1
prior1b = 1
prior1g = 1
prior1z = 1
prior1r = 1


####### Read in the start step values from a file...
filename1 <- "StartStepValuesM1V5E3.RData"
load( filename1 )    # Everything will be named correctly since it was saved from R
# Comment the line above if you don't have starting start step values in a file.


tic <- Sys.time()
nMCMC1 <- 30000
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
### These are the trace plots in the summplementary material (Figure 2)
#pdf("alpha.pdf")
plot( Q1MCMC$alpha1[,7], type = "l" , ylab = "alpha")  
plot( Q1MCMC$alpha1[,10], type = "l" , ylab = "alpha")
plot( Q1MCMC$alpha1[,11], type = "l" , ylab = "alpha")
#dev.off()

#pdf("beta.pdf")
plot( Q1MCMC$beta1, type = "l" , ylab = "beta")
#dev.off()

#pdf("zeta.pdf")
plot( Q1MCMC$zeta1, type = "l" , ylab = "zeta")
#dev.off()

#pdf("rho.pdf")
plot( Q1MCMC$rho1[,2], type = "l", ylab = "rho" )
#dev.off()


# Create files so you can read these in 
alpha1Step <- apply( Q1MCMC$alpha1, 2, sd )/4
gamma1Step <- apply( Q1MCMC$gamma1, 2, sd )/4
betaI1Step <- sd( Q1MCMC$betaI1)/4
beta1Step <- sd( Q1MCMC$beta1 )/4
zeta1Step <- sd( Q1MCMC$zeta1 )/4
rho1Step <- apply( Q1MCMC$rho1, 2, sd )/4
rho1IStep <- sd( Q1MCMC$rho1I)/4

# Step
#So we can append other generated values to the above list 

# Create files so you can read these in
alpha1 <- Q1MCMC$alpha1[nMCMC1,]
beta1 <- Q1MCMC$beta1[nMCMC1]
betaI1 <- Q1MCMC$betaI1[nMCMC1]
gamma1 <- Q1MCMC$gamma1[nMCMC1,]
zeta1 <- Q1MCMC$zeta1[nMCMC1]
rho1 <- Q1MCMC$rho1[nMCMC1,]
rho1I <- Q1MCMC$rho1I[nMCMC1]


#########  Using a save function. This will create a R database and save it.  When you load it everything
#########  will be there.
save(alpha1, beta1, gamma1, zeta1, rho1, rho1I, n1Q, chpt1,
     rchpt1, gchpt1, betaI1, ImpI1, ImpRI1,
     alpha1Step, beta1Step, betaI1Step, gamma1Step, zeta1Step, rho1Step, rho1IStep, 
     file = filename1 )

########  You can add a load statement above and it will pull all these in so you don't need to start all over.

Q1MCMC20 <- Q1MCMC
save( Q1MCMC20, file = "Q1MCMC20.Rdat" )


#################################################################################
#Predicted values of the variables
IPred1m <- apply( Q1MCMC$IPred1, 2, quantile, c(0.5,0.25,0.975) )
REPred1m <- apply( Q1MCMC$REPred1, 2, quantile, c(0.5,0.25,0.975) )
RIPred1m <- apply( Q1MCMC$RIPred1, 2, quantile, c(0.5,0.25,0.975) )
DPred1m <- apply( Q1MCMC$DPred1, 2, quantile, c(0.5,0.25,0.975) )
VPred1m <- apply( Q1MCMC$VPred1, 2, quantile, c(0.5,0.25,0.975) )

####  Plots in the main document are below:
## Plotting the variables
#pdf("Infected.pdf")
n1Q <- n1
plot( 1:n1Q, Q1$AdjInfect, col = "red", type = "b", lwd=2, ylab = " Actively Infected",
      xlab = "Days since 29Feb2020", main="Actively Infected")
#dev.off()
#The above is Figure 2a in the main document

#pdf("Recovered.pdf")
plot( 1:n1Q, Q1$Recovered, col = "seagreen", type = "b", lwd=2, ylab = " Recovered",
      xlab = "Days since 29Feb2020", main="Recovered")
#dev.off()
#The above is Figure 2b in the main document

#pdf("Deaths.pdf")
plot( 1:n1Q, Q1$Deaths, col = "black", type = "b", 
      lwd=2, xlab = "Days since 29Feb2020",ylab = "Deaths",main = "Deaths")
#dev.off()
#The above is Figure 2c in the main document

#pdf("Vaccinated.pdf")
plot( 1:n1Q, Q1$Vaccinated, col = "darkblue", type = "b", lwd=3,
      xlab = "Days since 29April2021",ylab = "Vaccinated",main = "Vaccinated")
#dev.off()
#The above is Figure 2d in the main document

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


##################################################################################################
# Effective Reproduction number....
##################################################################################################
load( "Q1MCMC20.Rdat" )
ls()
n1 <- n1Q
length(Q1MCMC20 )

X1a <- MatrixBuild1( chpt1, n1 )
X1g <- MatrixBuild1( gchpt1, n1 )
X1r <- MatrixBuild1( rchpt1, n1 )
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
#Rep0

#Below is the plot of the effective reproduction number in the supplementary material (Figure 1)
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


### Renaming the parameters:
alpha1 <- Q1MCMC20$alpha1
beta1 <- Q1MCMC20$beta1
betaI1 <- Q1MCMC20$betaI1
gamma1 <- Q1MCMC20$gamma1
rho1 <- Q1MCMC20$rho1
rho1I <- Q1MCMC20$rho1I
zeta1 <- Q1MCMC20$zeta1

##########################################################################################################
### The results below are in Table 4 of the main document.......
###########################################################################################################
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
###################################################################################################################

###########################################################################################################
#### Now, to find the difference between the alpha's: (Results in Table 5 )
############################################################################################################
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
######################################################################################################

############################################################################################
####%%%%%%%%%% Now, to find the difference between the gammas (Results in Table 6)
############################################################################################
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
########################################################################################

######################
#The end
#######################










