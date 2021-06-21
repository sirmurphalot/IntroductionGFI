#####################################################
######### Fiducial Binomial Simulation Study ########
############## n unknown, p unknown #################
#####################################################
########## Author: Alexander Murph, UNC #############
############## November 10th, 2020 ##################
#####################################################

### Uncomment these if you're working with SLURM cloud computing.
i = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
print(paste("Running seed:", i))
set.seed(i)
### Comment this if you're working with SLURM cloud computing.
# i = 1002
source("_simulation_helpers.R")
source("_belief_plaus_finder.R")

#### SET THE PARAMETERS
# Assuming your are on simulation i, generate the data.
p_values = c(0.1, 0.5, 0.9, 0.1, 0.5, 0.9)
true_ns  = c(15, 15, 15, 75, 75, 75)
MH_sigs = c(0.6,0.1,0.2,0.6,0.1,0.2)
mod_index = i%%6 + 1
MH_sigma = MH_sigs[mod_index]
p_values = c(p_values[mod_index])
TrueN = true_ns[mod_index]
data_sizes = c(100)
iterations = 16000


simulation_data = data.frame(TrueMu = NA,
                             TrueN = NA,
                             TrueP = NA,
                             DataSize = NA,
                             SimulationNumber = NA,
                             pois_belief_lower_bound_mu = NA,
                             pois_belief_upper_bound_mu = NA,
                             pois_belief_lower_bound_n = NA,
                             pois_belief_upper_bound_n = NA,
                             pois_belief_contains_truth = NA,
                             pois_plausability_lower_bound_mu = NA,
                             pois_plausability_upper_bound_mu = NA,
                             pois_plausability_lower_bound_n = NA,
                             pois_plausability_upper_bound_n = NA,
                             pois_plausability_contains_truth = NA,
                             nopois_belief_lower_bound_mu = NA,
                             nopois_belief_upper_bound_mu = NA,
                             nopois_belief_lower_bound_n = NA,
                             nopois_belief_upper_bound_n = NA,
                             nopois_belief_contains_truth = NA,
                             nopois_plausability_lower_bound_mu = NA,
                             nopois_plausability_upper_bound_mu = NA,
                             nopois_plausability_lower_bound_n = NA,
                             nopois_plausability_upper_bound_n = NA,
                             nopois_plausability_contains_truth = NA,
			     n_upper_CI = NA,
			     n_lower_CI = NA,
			     mu_upper_CI = NA,
			     mu_lower_CI = NA)

for(TrueP in p_values){
  for(sample_size in data_sizes){
    # Generate the true data
    TrueMu = TrueN * TrueP
    my_data = rbinom(sample_size, size = TrueN, prob = TrueP)
    # Calculate the GFD
    x = get_fiducial_dictionary(my_data, iterations, starting_method = "MAX", 
                                   epsilon = 1/mean(my_data), user_n = NULL, MH_sigma)
    ################################################################################################
    ## Note: If one is trying to replicate our results, you should make your own DataStore directory.
    ################################################################################################
    write.csv(x, file = paste("DataStore/fiddraw_",i,"_",TrueP*10,TrueN,".csv",sep=""))
    ################################################################################################
    ## Note: It is possible to break up this calculation into two files -- one that performs the data simulation (above)
    ##       and another that performs the belief/plausability box calculations (below).  One might want to consider this
    ##       based on computing limitations.  If you do, simply read the file in like so:
    ## x = read.csv(paste("../DataStore/fiddraw_",i,"_",TrueP*10,TrueN,".csv",sep=""))
    ##       and comment out (x = x[[1]]).
    ################################################################################################
    x = x[[1]]
    x = x[which(x$iteration_number %in% floor(iterations/2):iterations),]
    x$lower_mu = x$lower_p * x$n_value
    x$upper_mu = x$upper_p * x$n_value
    
    # Get the CIs
    n_upper_CI = mean(x$n_value > TrueN)
    n_lower_CI = mean(x$n_value < TrueN)
    mu_guesses = 0.5*x$lower_mu + 0.5*x$upper_mu
    mu_upper_CI = mean(mu_guesses > TrueMu)
    mu_lower_CI = mean(mu_guesses < TrueMu)

    # Get the bounds for the belief box
    beliefs = belief_bounds(x)
    # Bounds without screening for poisson
    pois_belief_lower_bound_mu = beliefs$mu_lower
    pois_belief_upper_bound_mu = beliefs$mu_upper
    pois_belief_lower_bound_n = beliefs$n_lower
    pois_belief_upper_bound_n = beliefs$n_upper
    pois_belief_contains_truth = (TrueMu <= pois_belief_upper_bound_mu) &
      (TrueMu >= pois_belief_lower_bound_mu) &
      (TrueN >= pois_belief_lower_bound_n) &
      (TrueN <= pois_belief_upper_bound_n)
   
    if(!pois_belief_contains_truth){
      write.csv(x,file=paste("MissedLogs/missedtruth_",TrueP*10,TrueN,i,".csv",sep=""))
    }
 
    # Get the bounds for the plausability box
    plausability = plausability_bounds(x)
    pois_plausability_lower_bound_mu = plausability$mu_lower
    pois_plausability_upper_bound_mu = plausability$mu_upper
    pois_plausability_lower_bound_n = plausability$n_lower
    pois_plausability_upper_bound_n = plausability$n_upper
    pois_plausability_contains_truth = (TrueMu <= pois_plausability_upper_bound_mu) &
      (TrueMu >= pois_plausability_lower_bound_mu) &
      (TrueN >= pois_plausability_lower_bound_n) &
      (TrueN <= pois_plausability_upper_bound_n)
    
    # Bounds with screening for poisson
    nopois_belief_lower_bound_mu = sort(x$lower_mu)[floor(length(x$lower_mu)*0.025)]
    remaining_uppers = x$upper_mu[-order(x$lower_mu)[1:floor(length(x$lower_mu)*0.025)]]
    upper_cutoff = length(remaining_uppers) - floor(length(x$upper_mu)*0.025)
    nopois_belief_upper_bound_mu = sort(remaining_uppers)[upper_cutoff]
    
    nopois_belief_lower_bound_n = quantile(x$n_value, probs = 0.025)
    nopois_belief_upper_bound_n = quantile(x$n_value, probs = 0.975)
    nopois_belief_contains_truth = (TrueMu <= pois_belief_upper_bound_mu) &
      (TrueMu >= pois_belief_lower_bound_mu) &
      (TrueN >= pois_belief_lower_bound_n) &
      (TrueN <= pois_belief_upper_bound_n)
    
    first_bound_mu = sort(x$upper_mu)[floor(length(x$upper_mu)*0.025)]
    remaining_lowers = x$lower_mu[-order(x$upper_mu)[1:floor(length(x$upper_mu)*0.025)]]
    upper_cutoff = length(remaining_lowers) - floor(length(x$lower_mu)*0.025)
    second_bound_mu = sort(remaining_lowers)[upper_cutoff]
    
    nopois_plausability_lower_bound_mu = min(first_bound_mu, second_bound_mu)
    nopois_plausability_upper_bound_mu = max(first_bound_mu, second_bound_mu)
    nopois_plausability_lower_bound_n = quantile(x$n_value, probs = 0.025)
    nopois_plausability_upper_bound_n = quantile(x$n_value, probs = 0.975)
    nopois_plausability_contains_truth = (TrueMu <= pois_plausability_upper_bound_mu) &
      (TrueMu >= pois_plausability_lower_bound_mu) &
      (TrueN >= pois_plausability_lower_bound_n) &
      (TrueN <= pois_plausability_upper_bound_n)

    # Record all results
    new_row = data.frame(TrueMu = TrueMu,
                         TrueN = TrueN,
                         TrueP = TrueP,
                         DataSize = sample_size,
                         SimulationNumber = i,
                         pois_belief_lower_bound_mu = pois_belief_lower_bound_mu,
                         pois_belief_upper_bound_mu = pois_belief_upper_bound_mu,
                         pois_belief_lower_bound_n = pois_belief_lower_bound_n,
                         pois_belief_upper_bound_n = pois_belief_upper_bound_n,
                         pois_belief_contains_truth = pois_belief_contains_truth,
                         pois_plausability_lower_bound_mu = pois_plausability_lower_bound_mu,
                         pois_plausability_upper_bound_mu = pois_plausability_upper_bound_mu,
                         pois_plausability_lower_bound_n = pois_plausability_lower_bound_n,
                         pois_plausability_upper_bound_n = pois_plausability_upper_bound_n,
                         pois_plausability_contains_truth = pois_plausability_contains_truth,
                         nopois_belief_lower_bound_mu = nopois_belief_lower_bound_mu,
                         nopois_belief_upper_bound_mu = nopois_belief_upper_bound_mu,
                         nopois_belief_lower_bound_n = nopois_belief_lower_bound_n,
                         nopois_belief_upper_bound_n = nopois_belief_upper_bound_n,
                         nopois_belief_contains_truth = nopois_belief_contains_truth,
                         nopois_plausability_lower_bound_mu = nopois_plausability_lower_bound_mu,
                         nopois_plausability_upper_bound_mu = nopois_plausability_upper_bound_mu,
                         nopois_plausability_lower_bound_n = nopois_plausability_lower_bound_n,
                         nopois_plausability_upper_bound_n = nopois_plausability_upper_bound_n,
                         nopois_plausability_contains_truth = nopois_plausability_contains_truth,
			 n_upper_CI = n_upper_CI,
                         n_lower_CI = n_lower_CI,
                         mu_upper_CI = mu_upper_CI,
                         mu_lower_CI = mu_lower_CI)
    simulation_data = rbind(simulation_data, new_row)
  }
}

if(i < 10){
 num = paste("000", i, sep = "")
} else if(i < 100) {
 num = paste("00", i, sep = "")
} else if(i < 1000){
 num = paste("0", i, sep = "")
}else{
 num = i
}

# Write each simulation's result to the BinData directory
################################################################################################
## Note: If one is trying to replicate our results, you should make your own BinData directory.
################################################################################################

write.csv(simulation_data, file = paste("BinData/binAllUnk_", num,".csv", sep = ""))


