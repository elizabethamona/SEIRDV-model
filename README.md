# SEIRDV model code

To run the main file "MCMCV1wCpp.R", you need to source the attached files called "SEIRDV6c.R" and "SEIRDV6.cpp" and also load "StartStepValuesM1V5E3.RData". Once everything is done correctly, you should be able to run this code successfully. Note that "SEIRDV6.cpp" is written in C++ in order to speed up the simulation time. 

Code Description:
This code estimates the SEIRDV model parameters using the Metropolis Hasting (MH) sampler and also solves the systems of ODE using Euler method. 

Datasets:
The datasets used for both the model analysis and model validation are attached. Save these datasets on your computer before attempting to run the code.


In addition, HTML file has been added for easy visualization of the plots in both the main paper and supplemetary material. 




