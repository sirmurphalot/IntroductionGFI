###########################################################
######### Fiducial Random Effects Simulation Study ########
############## Sigma, Beta unknown ########################
###########################################################
########## Author: Alexander Murph, UNC ###################
############## November 10th, 2020 ########################
###########################################################
library(mvtnorm)
library(geigen)
library(readr)
library(dplyr)
library(tidyr)
library(lmerTest)
library(BayesFactor)
library(rstan)

### Uncomment these if you're working with SLURM cloud computing.
# i = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
#print(paste("Running sigma value:", i))
# set.seed(i)
### Comment this if you're working with SLURM cloud computing.
i = 1
set.seed(i)

##############################################
# Parameters for this simulation:
d_list = c(6,6,3,5,6,4,6)
group1 = c(1,2,2,4,5,2,6)
group2 = c(1,2,5,4,10,2,6)
group3 = c(1,2,60,4,15,4,8)
group4 = c(1,2,0,8,20,6,8)
group5 = c(1,2,0,48,25,0,10)
group6 = c(100,100,0,0,30,0,10)
group_sizes = list(group1,group2,group3,
                   group4,group5,group6)
sigma2_a = c(.1,.5,1,.5,1,2,5,10)
sigma2_e = c(10,10,10,2,1,.5,.2,.1)

##############################################
# Simulating the Data
create_data = function(group_index, sigma_index){
  # Grab values for d, sigma a sigma e
  current_d = d_list[group_index]
  current_sig_a = sqrt(sigma2_a[sigma_index])
  current_sig_e = sqrt(sigma2_e[sigma_index])
  X = NULL
  y = NULL
  group_IDs = NULL
  current_index = 1
  
  for(index in 1:current_d){
    
    num_in_group = group_sizes[[index]][group_index]
    next_index = current_index + num_in_group - 1
    
    temp_X = matrix(1,nrow=num_in_group, ncol = 1)
    X = rbind(X,temp_X)
    y = c(y,rnorm(1,mean = 0,sd = current_sig_a) + 
      rnorm(num_in_group,mean = 0,sd = current_sig_e) )
    current_index = next_index + 1
    group_IDs = c(group_IDs, rep(index, times = num_in_group))
  }
  n = length(y)
  facID = as.factor(group_IDs)

  Sa=matrix(rep(0,n*n),ncol=n)
  for(my_fac in levels(facID)){
    temp_i=(facID==my_fac)
    Sa[temp_i,temp_i]=1
  }
  k=max(rowSums(Sa))
  return(list(n=n,d=1,k=k,X=X,Sa=Sa,y=y))
}

##############################################
# Setting up the STAN Model
fiducial_stan = stan_model(file='RandomIntercept.stan') 

##############################################
# Simulating and gathering metrics of interest
output_data = data.frame(value = NA, data = NA, metric = NA)
for(sig_num in 1:length(sigma_a)){
  for(groupings_num in 1:length(d_list)){
    stan_input=create_data(groupings_num, sig_num)
    FIDfit=sampling(fiducial_stan, data = stan_input, iter=5000L)
    samples = rstan::extract(FIDfit)
    data_info = paste("(",groupings_num,", ",sigma2_a[sig_num],", ",sigma2_e[sig_num],")", sep = "")
    
    sig2_es = samples$sigma2e
    sig2_as = samples$sigma2a
    sig2a_upper_CI = mean(sig2_as > (sigma2_a[sig_num]) )
    sig2e_upper_CI = mean(sig2_es > (sigma2_e[sig_num]))
    CI95_contains_siga = (quantile(sig2_as,probs=c(0.025,0.975))[1]<=sigma2_a[sig_num])&
      (quantile(sig_as,probs=c(0.025,0.975))[2]>=sigma2_a[sig_num])
    CI95_contains_sige = (quantile(sig2_es,probs=c(0.025,0.975))[1]<=sigma2_e[sig_num])&
      (quantile(sig_es,probs=c(0.025,0.975))[2]>=sigma2_e[sig_num])
    CI95_length_sige = quantile(sig2_es,probs=c(0.025,0.975))[2]-quantile(sig2_es,probs=c(0.025,0.975))[1] 
    CI95_length_siga = quantile(sig2_as,probs=c(0.025,0.975))[2]-quantile(sig2_as,probs=c(0.025,0.975))[1]

    new_row_1 = data.frame(value = sig2a_upper_CI, data = data_info, metric = "Sigma_a2 Upper CI")
    new_row_2 = data.frame(value = sig2e_upper_CI, data = data_info, metric = "Sigma_e2 Upper CI")
    new_row_3 = data.frame(value = CI95_contains_siga, data = data_info,
                           metric = "95% Sigma_a2 Containment")
    new_row_4 = data.frame(value = CI95_contains_sige, data = data_info, 
                           metric = "95% Sigma_e2 Containment")
    new_row_5 = data.frame(value = CI95_length_sige, data = data_info,
			   metric = "Length of 95% CI on Sigma_e2")
    new_row_6 = data.frame(value = CI95_length_siga, data = data_info,
                           metric = "Length of 95% CI on Sigma_a2")
    output_data = rbind(output_data,new_row_1,new_row_2,new_row_3,new_row_4,new_row_5,new_row_6)
  }
}
output_data = output_data[-1,]

##############################################
# Writing the data
if(i < 10){
  num = paste("000", i, sep = "")
} else if(i < 100) {
  num = paste("00", i, sep = "")
} else if(i < 1000){
  num = paste("0", i, sep = "")
}else{
  num = i
}
write.csv(output_data,file = paste("Data/simulation",num,".csv", sep = "_"))
