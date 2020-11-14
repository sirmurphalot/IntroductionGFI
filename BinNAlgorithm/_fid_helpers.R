####################################################
######## Fiducial Challenge Problem Methods ########
################# Helper Functions #################
########### Author: Alexander Murph, UNC ###########
############### November 4th, 2019 #################
####################################################
third_mean_numerical_bin = function(N, s){
  # Function that returns alpha-confidence interval using the third mean of the sth and
  # sth+1 uniform order statistics, with the third mean of a bernoulli distribution.
  # Author: Alexander Murph, UNC
  if(s==0){
    c_1=0
    c_2=0
  }else{
    c_1 = log(3) - lgamma(s) - lgamma(N-s) - (2/3)*lgamma(s) - (1/3)*lgamma(N-s) + 
      lgamma(s + (2/3)) + lgamma(N-s + (1/3))
    c_2 = log(3) - lgamma(s) - lgamma(N-s) - (1/3)*lgamma(s) - (2/3)*lgamma(N-s) + 
      lgamma(s + (1/3)) + lgamma(N-s + (2/3))
    c_1 = exp(c_1)
    c_2 = exp(c_2)
  }
  dist_draw_1 = rbeta(10000, shape1 = s+1, shape2 = N-s)
  dist_draw_2 = rbeta(10000, shape1 = s, shape2 = N-s+1)
  dist_draw_3 = rbeta(floor(10000*c_1), shape1 = s+(2/3), shape2 = N-s+(1/3))
  dist_draw_4 = rbeta(floor(10000*c_2), shape1 = s+(1/3), shape2 = N-s+(2/3))
  full_dist = c(dist_draw_1,dist_draw_2,dist_draw_3,dist_draw_4)
  
  return(full_dist)
}

third_mean_numerical_pois = function(N, s){
  # Function that returns alpha-confidence interval using the third mean of the sth and
  # sth+1 uniform order statistics of a gamma distribution, using the relation between exponentials
  # and poissons.
  if(s==0){
    c_1=0
    c_2=0
  }else{
    c_1 = log(3) + (2/3)*log(s) + lgamma(s + (1/3)) - lgamma(s + 1)
    c_2 = log(3) + (1/3)*log(s) + lgamma(s + (2/3)) - lgamma(s + 1)
    c_1 = exp(c_1)
    c_2 = exp(c_2)
  }
  dist_draw_1 = rgamma(10000, shape=s, rate=1)
  dist_draw_2 = rgamma(10000, shape=s+1, rate=1)
  dist_draw_3 = rgamma(floor(10000*c_1), shape=s+(1/3), rate=1)
  dist_draw_4 = rgamma(floor(10000*c_2), shape=s+(2/3), rate=1)
  full_dist = c(dist_draw_1,dist_draw_2,dist_draw_3,dist_draw_4)/n
  
  return(full_dist)
}

continuous_CDF_nonparam = function(x, my_data = x.nonparam){
  # Linear interpolation of empirical CDF with exponential tails, assuming continuity.
  # Author: Alexander Murph
  
  ordered_data = sort(my_data)
  if(x > max(ordered_data)){
    #do exponential upper tail
    return(1 - (1/(length(my_data) + 1)) * (1 / exp(x - ordered_data[length(my_data)])))
  } else if(min(which(ordered_data >= x)) == 1){
    #do exponential lower tail
    return( (exp(x)/exp(ordered_data[1]) - (1 / exp(ordered_data[1])) ) * (1/(length(my_data) + 1) ))
  } else {
    # do linear interpolation
    index_between = min(which(ordered_data >= x))
    slope = (1/(length(my_data) + 1)) / (ordered_data[index_between] - ordered_data[index_between-1])
    intercept = index_between * (1 / (length(my_data) + 1)) - slope * (ordered_data[index_between])
    return(slope*x + intercept)
  }
}

continuous_CDF_Nonparam = function(x, my_data = x.Nonparam){
  # Linear interpolation of empirical CDF with exponential tails, assuming continuity.
  # Author: Alexander Murph
  
  ordered_data = sort(my_data)
  if(x > max(ordered_data)){
    #do exponential upper tail
    return(1 - (1/(length(my_data) + 1)) * (1 / exp(x - ordered_data[length(my_data)])))
  } else if(min(which(ordered_data >= x)) == 1){
    #do exponential lower tail
    return( (exp(x)/exp(ordered_data[1]) - (1 / exp(ordered_data[1])) ) * (1/(length(my_data) + 1) ))
  } else {
    # do linear interpolation
    index_between = min(which(ordered_data >= x))
    slope = (1/(length(my_data) + 1)) / (ordered_data[index_between] - ordered_data[index_between-1])
    intercept = index_between * (1 / (length(my_data) + 1)) - slope * (ordered_data[index_between])
    return(slope*x + intercept)
  }
}

########################################
## Binomial, p know, helper functions ##
########################################

create_dictionaries = function(M) {
  # Method that creates the index dictionary and the probability mass dictionary.
  # Fills the probability mass dictionary with appropriate number of NA values.
  # Input: Length of range of n values considered, M
  # Output: index_dict, prob_dict
  ## Author: Alexander Murph
  index_dict = list()
  prob_dict = list()
  possible_indices = 1:M
  for( i in M:1 ){
    indices_for_this_size = list()
    NAs_for_this_size = list()
    starting_index = 0
    while(TRUE) {
      if((starting_index + i) > M){
        break
      }
      indices_for_this_size[[starting_index + 1]] = possible_indices[(starting_index + 1):(starting_index + i)]
      NAs_for_this_size[[starting_index + 1]] = rep(NA, times = i)
      starting_index = starting_index + 1
    }
    index_dict[[as.character(i)]] = indices_for_this_size
    prob_dict[[as.character(i)]] = NAs_for_this_size
  }
  return(list(index_dict, prob_dict) )
}

find_prob_mass = function(index_subset, range_of_n_values, index_dict, prob_dict, my_data,P){
  # Method for determining the probability mass of a given subset.  Assumes that all subsets of strictly greater
  # size have already been calculated.  
  # Input: index_subset, the indices of the n values that make up the set whose prob mass we want;
  #        range_of_n_values, the actual values of n.  Indices in former set refer to these actual values;
  #        index_dict, the dictionary of all index subsets considered for this algorithm;
  #        prob_dict, all calculated log probability masses for the former subsets of indices (and corresponding n values);
  #        my_data, the actual observed binomial(n,P) values -- used to calculate probability mass.
  # Output: prob_mass, probability mass of the given subset.
  ## Author: Alexander Murph
  min_size_considered = length(index_subset)
  max_size_considered = as.numeric(names(index_dict)[1])
  prob_mass = 0
  for(i in min_size_considered:max_size_considered){
    for(subset in 1:length(index_dict[[as.character(i)]]) ){
      if( identical(index_subset, index_dict[[as.character(i)]][[subset]]) ){
        subset_location = subset
        min_n = range_of_n_values[min(index_dict[[as.character(i)]][[subset]])]
        max_n = range_of_n_values[max(index_dict[[as.character(i)]][[subset]])]
        
        probability_canidates = pbinom(my_data, max_n, P) - pbinom(my_data-1, min_n, P)
        probability_canidates = ifelse(probability_canidates > 0, probability_canidates, 0)
        prob_mass = prod(probability_canidates)
      }else if( all(index_subset %in% index_dict[[as.character(i)]][[subset]]) ) {
        prob_mass = prob_mass - prob_dict[[as.character(i)]][[subset]]
      }
    }
  }
  #print(prob_mass)
  prob_mass = ifelse(prob_mass < 0, 0, prob_mass)
  #prob_dict[[as.character(length(index_subset))]][[subset_location]] = prob_mass
  return(prob_mass)
}

find_range_of_feasible_n = function(my_data, epsilon, P){
  # Method to find range of singleton values that are feasible based on chosen precision value epsilon.
  # Input: my_data, vector of observations from Bin(n,P) distribution with known P but unknown n
  #        epsilon, precision value associated with the lower cutoff for feasible probability mass
  # Output: range_of_n_values, range of n values that are feasible based off of precision cutoff.
  ## Author: Alexander Murph
  n_currently_considered = max(my_data)
  current_max_prob_mass = get_prob_mass_uncorrected(my_data, n_currently_considered, n_currently_considered, P)
  feasible_n_values = c(n_currently_considered)
  while(TRUE){
    n_currently_considered = n_currently_considered + 1
    current_prob_mass = get_prob_mass_uncorrected(my_data, n_currently_considered, n_currently_considered,P) 
    if( (log10(current_prob_mass) - log10(current_max_prob_mass)) >= 0){
      feasible_n_values = c(feasible_n_values, n_currently_considered)
      current_max_prob_mass = current_prob_mass
    } else if( (log10(current_prob_mass) - log10(current_max_prob_mass)) > epsilon ) {
      feasible_n_values = c(feasible_n_values, n_currently_considered)
    } else {
      break
    }
  }
  return(list(feasible_n_values, current_max_prob_mass))
}

get_prob_mass_uncorrected = function(my_data, min_n_value, max_n_value, P){
  # Method to get fiducial prob mass for a given subset candidate for n.
  # Input: my_data, vector of observations from Bin(n,P) distribution with known P but unknown n
  #        min_n_value, minimum n value of candidate subset of n in Bin(n,P)
  #        max_n_value, maximum n value of candidate subset of n in Bin(n,P)
  # Output: prob_mass, log of fiducial probability mass for given value of n.
  ## Author: Alexander Murph
  prob_mass = prod(pbinom(my_data, max_n_value, P) - pbinom(my_data-1, min_n_value, P))
  if(prob_mass < 0){
    prob_mass = 0
  }
  return(prob_mass)
}

n_subset_probabilities = function(my_data, epsilon,P){
  # Main method to find all probability masses for every possible subset of n
  # Input: my_data, vector of observations from Bin(n,P) distribution with known P but unknown n
  #        epsilon, precision value associated with the lower cutoff for feasible probability mass
  # Output: prob_dict, dictionary of probability masses organized in tree structure.
  #         index_dict, dictionary of the subsets associated with each probability mass
  #         n_values_considered, values of n associated with each index
  ## Author: Alexander Murph
  n_values_considered_list = find_range_of_feasible_n(my_data, epsilon,P)
  n_values_considered = n_values_considered_list[[1]]
  current_max_prob_mass = n_values_considered_list[[2]]
  dictionaries = create_dictionaries(length(n_values_considered))
  index_dict = dictionaries[[1]]
  prob_dict = dictionaries[[2]]
  total_prob_mass = 0
  for( index in names(index_dict) ){
    for( subset in index_dict[[index]] ){
      prob_mass = find_prob_mass(subset, n_values_considered, index_dict, prob_dict, my_data,P)
      prob_mass = ifelse(prob_mass > 0, prob_mass, 0)
      if( (log10(prob_mass) - log10(current_max_prob_mass)) >= 0){
        prob_dict[[index]][[min(subset)]] = prob_mass
        total_prob_mass = total_prob_mass + prob_mass
        current_max_prob_mass = prob_mass
      } else if( (log10(prob_mass) - log10(current_max_prob_mass)) > epsilon ) {
        prob_dict[[index]][[min(subset)]] = prob_mass
        total_prob_mass = total_prob_mass + prob_mass
      } else {
        prob_dict[[index]][[min(subset)]] = 0
      }
    }
  }
  f = function(x) { return(x / total_prob_mass) }
  for( index in names(index_dict) ){
    prob_dict[[index]] = lapply(prob_dict[[index]], f)
  }
  return(list(index_dict, prob_dict, n_values_considered))
}

find_quantiles = function(prob_list, n_values,P){
  # Assuming you have singleton exclusive probabilities for individual values of n, calculates precise quantiles.
  # Input: prob_list, probability mass associated with each value in n_values
  #        n_values, list of feasible values for n
  # Output: the lower and upper quantiles.
  ## Author: Alexander Murph  
  lower_2_5th = 0
  middle = 0
  upper_2_5th = 0
  total_prob = 0
  updated_upper = F
  updated_lower = F
  updated_middle = F
  
  for(index in 1:length(prob_list)){
    total_prob = total_prob + prob_list[index]
    if( (total_prob >= 0.025) & (!updated_lower)){
      updated_lower = T
      difference = total_prob - 0.025
      lower_2_5th = n_values[index] - difference/prob_list[index]
    } else if( (total_prob >= 0.5) & (!updated_middle) ) {
      updated_middle = T
      difference = total_prob - 0.5
      middle = n_values[index] - difference/prob_list[index]
    }  else if((total_prob >= 0.975) & (!updated_upper)){
      updated_upper = T
      difference = total_prob - 0.975
      upper_2_5th = n_values[index] - difference/prob_list[index]
    }
  }
  return(list(lower_2_5th, middle, upper_2_5th))
}

generate_density = function(dictionaries, iterations,P){
  # Main method to simulate a distribution from the dictionary of probability masses the above creates.
  # Input: dictionaries, the filled-in dictionaries from the above methods
  # Output: list of draws from distribution for n.
  ## Author: Alexander Murph
  prob_list = rep(0, times = length(dictionaries[[3]]))
  n_values = dictionaries[[3]]
  
  for(size_index in names(dictionaries[[2]])){
    print(paste("Finding size:", size_index))
    for(index in 1:length(dictionaries[[2]][[size_index]])){
      if(size_index == '1'){
        prob_list[dictionaries[[1]][[size_index]][[index]]] = prob_list[dictionaries[[1]][[size_index]][[index]]] + 
          dictionaries[[2]][[size_index]][[index]]
      } else {
        prob_list[dictionaries[[1]][[size_index]][[index]][1]] = prob_list[dictionaries[[1]][[size_index]][[index]][1]] + 
          0.5*dictionaries[[2]][[size_index]][[index]] 
        prob_list[dictionaries[[1]][[size_index]][[index]][as.numeric(size_index)]] = 
          prob_list[dictionaries[[1]][[size_index]][[index]][as.numeric(size_index)]] + 0.5*dictionaries[[2]][[size_index]][[index]]
      }
    }
  }
  cumulative_prob = 0
  for(index in 1:length(prob_list)) {
    cumulative_prob = cumulative_prob + prob_list[index]
    prob_list[index] = cumulative_prob
  }
  draws_from_distn = c()
  for(iter in 1:iterations){
    unif_draw = runif(1)
    draws_from_distn = c(draws_from_distn, n_values[min(which(prob_list >= unif_draw))])
  }
  return(draws_from_distn)
}

### Bayesian solution helper function
generate_bayesian_posterior = function(my_data, epsilon,P, iterations){
  # Main method to find all probability masses for every possible subset of n
  # Input: my_data, vector of observations from Bin(n,P) distribution with known P but unknown n
  #        epsilon, precision value associated with the lower cutoff for feasible probability mass
  # Output: prob_dict, dictionary of probability masses organized in tree structure.
  #         index_dict, dictionary of the subsets associated with each probability mass
  #         n_values_considered, values of n associated with each index
  ## Author: Alexander Murph
  n_currently_considered = max(my_data)
  current_max_prob_mass = get_prob_mass_uncorrected(my_data, n_currently_considered, n_currently_considered, P)
  feasible_n_values = c(n_currently_considered)
  prob_list = c(current_max_prob_mass)
  total_prob_mass = current_max_prob_mass
  while(TRUE){
    n_currently_considered = n_currently_considered + 1
    current_prob_mass = get_prob_mass_uncorrected(my_data, n_currently_considered, n_currently_considered,P) 
    if( (log10(current_prob_mass) - log10(current_max_prob_mass)) >= 0){
      feasible_n_values = c(feasible_n_values, n_currently_considered)
      prob_list = c(prob_list, current_prob_mass)
      current_max_prob_mass = current_prob_mass
      total_prob_mass = total_prob_mass + current_prob_mass
    } else if( (log10(current_prob_mass) - log10(current_max_prob_mass)) > epsilon ) {
      feasible_n_values = c(feasible_n_values, n_currently_considered)
      prob_list = c(prob_list, current_prob_mass)
      total_prob_mass = total_prob_mass + current_prob_mass
    } else {
      break
    }
  }
  f = function(x) { return( log(x) - log(total_prob_mass) ) }
  prob_list = exp(unlist(lapply(prob_list,f)))
  cumulative_prob = 0
  cum_prob_list = rep(0, times = length(prob_list))
  for(index in 1:length(prob_list)) {
    cumulative_prob = cumulative_prob + prob_list[index]
    cum_prob_list[index] = cumulative_prob
  }
  draws_from_distn = c()
  for(iter in 1:iterations){
    unif_draw = runif(1)
    draws_from_distn = c(draws_from_distn, feasible_n_values[min(which(cum_prob_list >= unif_draw))])
  }
  #return(list(feasible_n_values, prob_list))
  return(draws_from_distn)
}

