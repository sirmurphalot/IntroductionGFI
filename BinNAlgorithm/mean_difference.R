########################################################
######### Fiducial Binomial Simulation Coverage ########
################# n unknown, p known# ##################
########################################################
########## Author: Alexander Murph, UNC ################
############## November 10th, 2019 #####################
########################################################

source("_fid_helpers.R")
P_values = c(0.01, 0.05, 0.1, 0.2, 0.4, 0.5, 
             0.6, 0.7, 0.8, 0.9, 0.95, 0.99)
### Uncomment these if you're working with SLURM cloud computing.
# i = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
# print(paste("Running seed:", i))
# set.seed(i)
### Comment this if you're working with SLURM cloud computing.
i = 1
set.seed(i)
generate_fid_values = function(bin_data, num_of_draws, P){
  # Method that generates N values from the posterior-like fiducial distribution based
  # on the observed data.
  ## Author: Alexander Murph
  dictionaries = n_subset_probabilities(bin_data, -10,P)
  prob_list = rep(0, times = length(dictionaries[[3]]))
  n_values = dictionaries[[3]]
  fid.dist = generate_density(dictionaries, num_of_draws,P)
  return(fid.dist)
}

my_data = data.frame(P = NA, MAD = NA, parameter = NA, data_size = NA, upper_p = NA, lower_p = NA)
TrueN = 10
num_of_draws = 1000
size_of_data = c(10, 50, 100)

###### BEGIN SIMULATION
for(value in P_values){
  for(my_size in size_of_data){
    # Generate new binomial data
    bin_data = rbinom(my_size, size = TrueN, prob = value)
    # Use (deterministic) algorithm to generate GFD.
    fid_dist_for_n = generate_fid_values(bin_data, num_of_draws,value)
    # Generate Bayes posterior for comparison
    bayes_dist_for_n = generate_bayesian_posterior(bin_data, -10, value, num_of_draws)
    # Add results to our data, and repeat.
    temp_data = data.frame(P = rep(value, times = 2),
                           MAD = c( mean(abs(fid_dist_for_n - TrueN)), mean(abs(bayes_dist_for_n - TrueN))),
                           parameter = c(rep("fiducial-n", times = 1), 
                                         rep("bayesian-n", times = 1)),
                           data_size = rep(my_size, times = 2),
			   upper_p = c(mean(fid_dist_for_n < TrueN), mean(bayes_dist_for_n < TrueN)),
			   lower_p = c(mean(fid_dist_for_n > TrueN), mean(bayes_dist_for_n > TrueN)))
    my_data = rbind(my_data, temp_data)
    print("------")
    print(paste(value, ", ", my_size, " finished", sep = ""))
  }
}
###### SAVE RESULTS
write.csv(my_data, file=paste("binData/binomial_mean_difference_", i, ".csv", sep = ""))
