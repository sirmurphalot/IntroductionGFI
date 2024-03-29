#####################################################
######### Fiducial Binomial Simulation Study ########
############## n unknown, p unknown #################
#####################################################
########## Author: Alexander Murph, UNC #############
############## November 10th, 2020 ##################
#####################################################

library(ggplot2)
library(latex2exp)
library(ggExtra)
library(gridExtra)
library(knitr)
library(mvtnorm)
library(matrixcalc)
library(gtools)

###############################
##### IMPORTANT NOTE
###############################
##### All descriptions of these methods is available in the Pseudocode
##### titled AllBinomialPseudo2020.pdf.  We would highly suggest that a
##### curious reader look here to understand the following methods.
generate_first_uniforms = function(my_data, first_n, first_mu){
  n_hat = first_n
  proposed_p = first_mu / first_n
  ordered_data = sort(my_data)
  upper_uniforms = pbeta(1 - proposed_p, shape1 = n_hat - ordered_data,
                         shape2 = ordered_data + 1)
  lower_uniforms = pbeta(1 - proposed_p, shape1 = n_hat - ordered_data + 1,
                         shape2 = ordered_data)
  generate_uniforms = (upper_uniforms + lower_uniforms) / 2
  return(generate_uniforms)
}

shift_n_uniforms = function(my_data, new_mu, new_n){
  # Given the observed data, and mu&n values you want to be true, find a set of uniforms
  # that give a solution set that includes mu&n.
  n_hat = new_n
  proposed_p = new_mu / new_n
  ordered_data = sort(my_data)
  upper_uniforms = pbeta(1 - proposed_p, shape1 = n_hat - ordered_data,
                         shape2 = ordered_data + 1)
  lower_uniforms = pbeta(1 - proposed_p, shape1 = n_hat - ordered_data + 1,
                         shape2 = ordered_data)
  return(pairlist(lowers = lower_uniforms, uppers = upper_uniforms))
}

fix_ties = function(uniforms, my_data, index_of_change){
  indices_of_ties = which(my_data == my_data[index_of_change])
  indices_of_ties = indices_of_ties[which(indices_of_ties <= index_of_change)]
  uniforms[indices_of_ties] = sort(uniforms[indices_of_ties])
  return(uniforms)
}


find_u_bounds = function(temp_run, data_value){
  upper_uniforms = c()
  lower_uniforms = c()
  for(sub_index in 1:nrow(temp_run)){
    alpha_debug = ifelse(temp_run$n_value[sub_index] - data_value < 1, 
                         0, temp_run$n_value[sub_index] - data_value)
    lowers = pbeta(1 - temp_run$upper_p[sub_index], 
                   shape1 = temp_run$n_value[sub_index] - data_value + 1,
                   shape2 = data_value )
    uppers = pbeta(1 - temp_run$lower_p[sub_index], 
                   shape1 = alpha_debug,
                   shape2 = data_value + 1)
    upper_uniforms = c(upper_uniforms, uppers)
    lower_uniforms = c(lower_uniforms, lowers)
  }
  return(pairlist(upper = max(upper_uniforms), lower = min(lower_uniforms)))
}

MH_step_n = function(my_data, ending_n_hat, epsilon_constant, temp_run, my_max, epsilon){
  ## Metropolis-Hastings n Step
  
  if(nrow(temp_run) ==0){
    return(list(my_data, ending_n_hat, epsilon_constant))
  }
  
  old_n_value = temp_run$n_value[1]
  lower_mu = temp_run$n_value[1] * temp_run$lower_p[1]
  upper_mu = temp_run$n_value[1] * temp_run$upper_p[1]
  new_mu = runif(1, min = lower_mu, max = upper_mu)
  new_n = ifelse(rbinom(1,size=1,prob=0.5), old_n_value - 1, old_n_value + 1)
  
  if(new_n < max(my_data$value)){
    # We do not create a new data because the proposed n was not valid.
    new_data = data.frame(value= character(),
                          number=character(), 
                          min_uniform=character(), 
                          max_uniform = character(),
                          stringsAsFactors = FALSE)
    
    # Otherwise find new uniforms that allow for this solution:
  } else {
    orig_data = NULL
    # Grab the original data values from the compact storage we use
    for(i in 1:nrow(my_data)){
      orig_data = c(orig_data, rep(my_data$value[i], times = my_data$number[i]))
    }
    
    # Get the theoretical upper and lower bounds based on the DGE, using the
    # data and the new n, mu.
    new_uniform_bounds = shift_n_uniforms(orig_data, new_mu, new_n)
    new_uniforms = NULL
    
    # Use these bounds to draw new uniforms.  Based on how we drew from shift_n_uniforms,
    # new n, mu should be somewhere in this solution set.
    for(i in 1:length(new_uniform_bounds$lowers)){
      new_uniforms = c(new_uniforms, runif(1, min = new_uniform_bounds$lowers[i],
                                           max = new_uniform_bounds$uppers[i]))
    }
    # Just in case, I'll make sure there are no misorderings.  Theoretically, this shouldn't happen, though.
    new_uniforms = sort(new_uniforms)
    
    # Just hold onto the max and the min for each value/number pair
    # This reflects the usual way we store uniform values for later use.
    temp_my_data = my_data
    lower_uniform_index = 1
    upper_uniform_index = 1
    for(i in 1:nrow(my_data)){
      upper_uniform_index = upper_uniform_index + my_data$number[i] - 1
      temp_my_data$min_uniform[i] = min(new_uniforms[lower_uniform_index:upper_uniform_index])
      temp_my_data$max_uniform[i] = max(new_uniforms[lower_uniform_index:upper_uniform_index])
      lower_uniform_index = upper_uniform_index + 1
      upper_uniform_index = lower_uniform_index
    }
    new_data = temp_my_data
    # Find the new solution set with these new uniforms (should exist).
    temp_temp_run = generate_p_sets(temp_my_data, 1, my_max, epsilon, final_n = ending_n_hat, 
                                    epsilon_constant)
    new_epsilon_constant = temp_temp_run$epsilon_constant
    new_n_hat = temp_temp_run$ending_n
    temp_temp_run = temp_temp_run$pset_data
    
    if(nrow(temp_temp_run)==0){
      new_data = data.frame(value= character(),
                            number=character(), 
                            min_uniform=character(), 
                            max_uniform = character(),
                            stringsAsFactors = FALSE)
    }else if(temp_temp_run$n_value[1] != new_n ){
      # We need the smallest n to be the proposed new n.  If it isn't, 
      # mark this so we can observe how likely this is.
      # print("Proposed n wasn't the smallest n in the new solution set!")
      new_data = data.frame(value= character(),
                            number=character(), 
                            min_uniform=character(), 
                            max_uniform = character(),
                            stringsAsFactors = FALSE)
    } else {
      # Now we have new uniforms and a new solution set.  Let's find the transition probability:
      # We already have everything we need for old->new
      MH_log_denominator = -log(2) - log(upper_mu - lower_mu) - 
        sum( log(new_uniform_bounds$uppers - new_uniform_bounds$lowers))
      
      # To get new->old:
      lower_mu = temp_temp_run$n_value[1] * temp_temp_run$lower_p[1]
      upper_mu = temp_temp_run$n_value[1] * temp_temp_run$upper_p[1]
      old_mu = new_mu
      # Use the DGE to find the uniform bounds for our original data:
      old_uniform_bounds = shift_n_uniforms(orig_data, old_mu, old_n_value)
      MH_log_numerator = -log(2) - log(upper_mu - lower_mu) - 
        sum(log(old_uniform_bounds$uppers - old_uniform_bounds$lowers))
    }
    # if(prod(c(temp_my_data$min_uniform==0,temp_my_data$max_uniform==0)) == 1){
    #   print("this happened")
    #   print(MH_log_numerator)
    #   print(MH_log_denominator)
    # }
  }
  
  
  if(nrow(new_data) == 0){
    did_accept = FALSE
  } else {
    accept_ratio = MH_log_numerator - MH_log_denominator
    coin_flip = log(runif(1))
    if(is.nan(accept_ratio)){
      did_accept = FALSE
    } else if(coin_flip < accept_ratio){
      # In the case that we accept, update our data (uniform values that give the new solution set)
      my_data = temp_my_data
      ending_n_hat = new_n_hat
      epsilon_constant = new_epsilon_constant
      did_accept = TRUE
    } else{
      did_accept = FALSE
    }
  }
  
  return(list(my_data, ending_n_hat, epsilon_constant))
}


MH_step_mu = function(my_data, ending_n_hat, epsilon_constant, temp_run, my_max,epsilon,sigma){
  ## Metropolis-Hastings mu Step
  if(nrow(temp_run) ==0){
    return(list(my_data, ending_n_hat, epsilon_constant))
  }
  old_n_value = temp_run$n_value[1]
  lower_mu = temp_run$n_value[1] * temp_run$lower_p[1]
  upper_mu = temp_run$n_value[1] * temp_run$upper_p[1]
  # Randomly choose a mu, then perturb it based on a N(0,sigma) distribution
  new_mu = runif(1, min = lower_mu, max = upper_mu) + rnorm(1, mean = 0, sd = sigma)
  
  # Make sure that this new mu selected is valid.
  new_p = new_mu/old_n_value
  new_n = old_n_value
  
  ## Record that we didn't accept if this is not a valid new_p.
  if((new_p <= 0)){
    # print("new p was negative!")
    return(list(my_data, ending_n_hat, epsilon_constant))
  } else if((new_p <= 0)){
    # print("new p was negative!")
    return(list(my_data, ending_n_hat, epsilon_constant))
  } else if((new_p >= 1)){
    # print("new p was above 1!")
    return(list(my_data, ending_n_hat, epsilon_constant))
  } else if((new_p >= 1)){
    # print("new p was above 1!")
    return(list(my_data, ending_n_hat, epsilon_constant))
  }
  
  orig_data = NULL
  for(i in 1:nrow(my_data)){
    orig_data = c(orig_data, rep(my_data$value[i], times = my_data$number[i]))
  }
  
  # Get the theoretical upper and lower bounds.  This is the same method used for the shift n --
  # it takes in our data, and a wanted mu&n, and spits out uniforms such that there's a soln.
  new_uniform_bounds = shift_n_uniforms(orig_data, new_mu, new_n)
  new_uniforms = NULL
  
  # Use these bounds to draw new uniforms
  for(i in 1:length(new_uniform_bounds$lowers)){
    new_uniforms = c(new_uniforms, runif(1, min = new_uniform_bounds$lowers[i],
                                         max = new_uniform_bounds$uppers[i]))
  }
  # Just in case, I'll make sure there are no mismatches here (there theoretically shouldn't be)
  new_uniforms = sort(new_uniforms)
  
  # Just hold onto the max and the min for each value/number pair.  This is how we originally stored
  # this data.
  temp_my_data = my_data
  lower_uniform_index = 1
  upper_uniform_index = 1
  for(i in 1:nrow(my_data)){
    upper_uniform_index = upper_uniform_index + my_data$number[i] - 1
    temp_my_data$min_uniform[i] = min(new_uniforms[lower_uniform_index:upper_uniform_index])
    temp_my_data$max_uniform[i] = max(new_uniforms[lower_uniform_index:upper_uniform_index])
    lower_uniform_index = upper_uniform_index + 1
    upper_uniform_index = lower_uniform_index
  }
  new_data = temp_my_data
  # Find the new solution set with these new uniforms.
  temp_temp_run = generate_p_sets(temp_my_data, 1, my_max, epsilon, final_n = ending_n_hat, 
                                  epsilon_constant)
  new_epsilon_constant = temp_temp_run$epsilon_constant
  new_n_hat = temp_temp_run$ending_n
  temp_temp_run = temp_temp_run$pset_data
  
  if(nrow(temp_temp_run)==0){
    new_data = data.frame(value= character(),
                          number=character(), 
                          min_uniform=character(), 
                          max_uniform = character(),
                          stringsAsFactors = FALSE)
  } else if(temp_temp_run$n_value[1] != new_n ){
    # We need the smallest n to be the proposed new n.  If it isn't, 
    # mark this so we can observe how likely this is.
    new_data = data.frame(value= character(),
                          number=character(), 
                          min_uniform=character(), 
                          max_uniform = character(),
                          stringsAsFactors = FALSE)
  } else {
    # Now we have new uniforms and a new solution set.  Let's find the transition probability:
    # We already have everything we need for old->new
    MH_log_denominator = - log(upper_mu - lower_mu) - 
      sum( log(new_uniform_bounds$uppers - new_uniform_bounds$lowers))
    
    # To get new->old:
    lower_mu = temp_temp_run$n_value[1] * temp_temp_run$lower_p[1]
    upper_mu = temp_temp_run$n_value[1] * temp_temp_run$upper_p[1]
    old_mu = new_mu
    old_uniform_bounds = shift_n_uniforms(orig_data, old_mu, old_n_value)
    MH_log_numerator = - log(upper_mu - lower_mu) - 
      sum(log(old_uniform_bounds$uppers - old_uniform_bounds$lowers))
    
    # NOTE: By the symmetry present in the univariate normal, the probabillity of this should
    # cancel in the ratio (so I have not included it).
  }
  
  
  if(nrow(new_data) == 0){
    did_accept = FALSE
  } else {
    accept_ratio = MH_log_numerator - MH_log_denominator
    coin_flip = log(runif(1))
    if(coin_flip < accept_ratio){
      my_data = temp_my_data
      ending_n_hat = new_n_hat
      epsilon_constant = new_epsilon_constant
      did_accept = TRUE
    } else{
      did_accept = FALSE
    }
  }
  
  return(list(my_data, ending_n_hat, epsilon_constant))
}


update_uniforms = function(my_data, first_n, epsilon, ending_n_hat, sigma, epsilon_constant){
  my_max = first_n
  
  for(index in 1:nrow(my_data)){
    fix_max = rbinom(1, size = 1, prob = 0.5)
    
    if(my_data$number[index] == 1){
      temp_my_data = my_data
      temp_my_data = my_data[-index,]
      temp_run = generate_p_sets(temp_my_data, 1, my_max, epsilon, final_n = ending_n_hat, 
                                 epsilon_constant)
      epsilon_constant = temp_run$epsilon_constant
      ending_n_hat = temp_run$ending_n
      temp_run = temp_run$pset_data
      
      uniform_bounds = find_u_bounds(temp_run, my_data$value[index])
      min_bound = uniform_bounds$lower
      max_bound = uniform_bounds$upper
      new_uniform = runif(1, min = min_bound, max = max_bound)
      if(!is.na(new_uniform)){
        my_data$min_uniform[index] = new_uniform
        my_data$max_uniform[index] = new_uniform
      }
    }else if(fix_max){
      temp_my_data = my_data
      old_min = my_data$min_uniform[index]
      temp_my_data$min_uniform[index] = my_data$max_uniform[index]
      max_bound = my_data$max_uniform[index]
      temp_run = generate_p_sets(temp_my_data, 1, my_max, epsilon, final_n = ending_n_hat, 
                                 epsilon_constant)
      epsilon_constant = temp_run$epsilon_constant
      ending_n_hat = temp_run$ending_n
      temp_run = temp_run$pset_data
      
      uniform_bounds = find_u_bounds(temp_run, my_data$value[index])
      min_bound = uniform_bounds$lower
      
      # Draw a new minimum uniform order statistic from the data that is not fixed: m - 1 observations.
      # This is the following beta.
      min_unif = rbeta(1, shape1 = 1, shape2 = my_data$number[index] - 1) * (max_bound - min_bound) + min_bound
      if(!is.na(min_unif)){
        my_data$min_uniform[index] = min_unif
      }
    } else {
      temp_my_data = my_data
      old_max = my_data$max_uniform[index]
      temp_my_data$max_uniform[index] = my_data$min_uniform[index]
      min_bound = my_data$min_uniform[index]
      temp_run = generate_p_sets(temp_my_data, 1, my_max, epsilon, final_n = ending_n_hat, 
                                 epsilon_constant)
      epsilon_constant = temp_run$epsilon_constant
      ending_n_hat = temp_run$ending_n
      temp_run = temp_run$pset_data
      
      uniform_bounds = find_u_bounds(temp_run, my_data$value[index])
      max_bound = uniform_bounds$upper
      
      # Draw a new maximum uniform order statistic from the data that is not fixed: m - 1 observations.
      # This is the following beta.
      max_unif = rbeta(1, shape1 = my_data$number[index] - 1, shape2 = 1) * (max_bound - min_bound) + min_bound
      if(!is.na(max_unif)){
        my_data$max_uniform[index] = max_unif
      }
    }
  }
  
  # Perform the MH n update step:
  update_temp = MH_step_n(my_data, ending_n_hat, epsilon_constant, temp_run, my_max, epsilon)
  my_data = update_temp[[1]]
  ending_n_hat = update_temp[[2]]
  epsilon_constant = update_temp[[3]]
  
  # Perform the MH mu update step:
  update_temp = MH_step_mu(my_data, ending_n_hat, epsilon_constant, temp_run, my_max, epsilon, sigma)
  my_data = update_temp[[1]]
  ending_n_hat = update_temp[[2]]
  epsilon_constant = update_temp[[3]]
  # 
  
  return(pairlist(the_data = my_data, ending_n = ending_n_hat, epsilon_constant = epsilon_constant))
}

find_p_bounds = function(my_data, n_hat){
  lower_ps = c()
  upper_ps = c()
  for(index in 1:nrow(my_data)){
    lower_bound_debug = ifelse(n_hat - my_data$value[index] < 1, 0, n_hat - my_data$value[index])
    new_upper = qbeta( 1 - my_data$min_uniform[index], shape1 = my_data$value[index] + 1, 
                       shape2 = lower_bound_debug)
    new_lower = qbeta(1 - my_data$min_uniform[index], shape1 = my_data$value[index], 
                      shape2 = n_hat - my_data$value[index] + 1)
    lower_ps = c(lower_ps, new_lower)
    upper_ps = c(upper_ps, new_upper)
    
    new_upper = qbeta( 1 - my_data$max_uniform[index], shape1 = my_data$value[index] + 1, 
                       shape2 = lower_bound_debug)
    new_lower = qbeta(1 - my_data$max_uniform[index], shape1 = my_data$value[index], 
                      shape2 = n_hat - my_data$value[index] + 1)
    lower_ps = c(lower_ps, new_lower)
    upper_ps = c(upper_ps, new_upper)
    
  }
  return(pairlist(lower = lower_ps, upper = upper_ps))
}

close_to_poisson = function(lower_ps, upper_ps, n_hat, my_data, epsilon){
  beta_mu_lowers = lower_ps * n_hat
  beta_mu_uppers = upper_ps * n_hat
  gamma_mu_lowers = c()
  gamma_mu_uppers = c()
  for(i in 1:nrow(my_data)){
    gamma_mu_lowers = c(gamma_mu_lowers, qgamma(1 - my_data$min_uniform[i], shape = my_data$value[i], 
                                                rate = 1))
    gamma_mu_uppers = c(gamma_mu_uppers, qgamma(1 - my_data$min_uniform[i], shape = my_data$value[i] + 1, 
                                                rate = 1))
    gamma_mu_lowers = c(gamma_mu_lowers, qgamma(1 - my_data$max_uniform[i], shape = my_data$value[i], 
                                                rate = 1))
    gamma_mu_uppers = c(gamma_mu_uppers, qgamma(1 - my_data$max_uniform[i], shape = my_data$value[i] + 1, 
                                                rate = 1))
  }
  all_values = c(gamma_mu_lowers,gamma_mu_uppers,gamma_mu_lowers,gamma_mu_uppers)
  is_infinite = sapply(all_values, is.infinite)
  if(sum(is_infinite)){
    return(TRUE)
  }
  differences_lower = abs(beta_mu_lowers - gamma_mu_lowers)
  differences_upper = abs(beta_mu_uppers - gamma_mu_uppers)
  differences = sum(c((differences_lower > epsilon), differences_upper > epsilon))
  return(differences < 1)
}

generate_p_sets = function(my_data, iter, first_n, epsilon, final_n = NULL, epsilon_constant = 1){
  temp_data = data.frame(lower_p = NA, upper_p = NA, n_value = NA, iteration_number = NA)
  n_hat = max(my_data$value)
  not_found_first = TRUE
  while(TRUE){
    if(n_hat%%1000 ==0){
      # print("Searched over 1000 n_hats")
      # print(close_to_poisson(lower_ps, upper_ps, n_hat, my_data, epsilon_constant*epsilon))
      # print(my_data)
    }
    next_run = find_p_bounds(my_data, n_hat)
    lower_ps = next_run$lower
    upper_ps = next_run$upper
    # Sometimes this strange error comes up.  My solution is to just end the generation of p_sets.
    # This may require a closer look if it means a solution that was expected does not occur.
    if(anyNA(lower_ps)|anyNA(upper_ps)){
      # print("NAs showed up in p_set generation")
      ending_n_hat = n_hat
      break
    }
    if( (n_hat > 1000)&(iter == 0) ){
      # During the first time through, we tune the epsilon constant so that we do not 
      # search an n space that is absurdly large.  This does mean that we are not allowing
      # inference on n above 1000.
      # n can go higher during later iterations, we just start epsilon at such a place that
      # this does not initially happen.
      epsilon_constant = epsilon_constant + 50
      all_info = generate_p_sets(my_data, iter, first_n, epsilon, final_n, epsilon_constant)
      temp_data = rbind(temp_data, all_info$pset_data)
      ending_n_hat = all_info$ending_n
      epsilon_constant = all_info$epsilon_constant
      break
    } else{
      # Assuming the epsilon constant needs no tuning, we proceed with exploring the n space:
      
      # If there is no solution, simply check the next n.
      if( (max(lower_ps) > min(upper_ps)) ){
        n_hat = n_hat + 1
      } else{
        # If there is a solution, add it to the solution set, and proceed to the next n.
        new_row = data.frame(lower_p = max(lower_ps), upper_p = min(upper_ps), 
                             n_value = n_hat, iteration_number = iter)
        temp_data = rbind(temp_data, new_row)
        n_hat = n_hat + 1
      }
      if( (iter == 0) ){
        if( (close_to_poisson(lower_ps, upper_ps, n_hat, my_data, epsilon_constant*epsilon) & (n_hat >= first_n))){
          ending_n_hat = n_hat
          break
        }
      } else if( n_hat > 2000) {
        ending_n_hat = 2000
        break
      } else {
        if(length(final_n)==0){
          next
        }else if( (n_hat >= final_n) &
                  (n_hat >= first_n) & 
                  (close_to_poisson(lower_ps, upper_ps, n_hat, my_data, epsilon_constant*epsilon)) ){
          ending_n_hat = n_hat
          break
        }
      }
    } 
    
    
  }
  temp_data = temp_data[-1,]
  return(list(pset_data = temp_data, ending_n = ending_n_hat, epsilon_constant = epsilon_constant))
}

get_fiducial_dictionary = function(my_data, iterations, starting_method = c("MLE", "MOM", "MAX", 
                                                                            "Carroll-Lombard", 
                                                                            "DasGupta-Rubin", 
                                                                            "USER"), 
                                   epsilon = 0.5, user_n = NULL, sigma){
  ## Make sure that we have a fresh start to recording the acceptance ratio and uniform logs:
  
  if(starting_method == "USER"){
    first_mu = user_n[[2]]
    user_n = user_n[[1]]
  } else {
    first_mu = mean(my_data)
  }
  
  m = length(my_data)
  my_data = sort(my_data)
  unique_data = unique(my_data)
  repeats = function(x){return(sum(my_data == x))}
  repeat_observations = sapply(unique_data, repeats)
  
  solution_data = data.frame(lower_p = NA, upper_p = NA, n_value = NA, iteration_number = NA)
  
  if(starting_method == "MLE"){
    f_likelihood = function(n) {
      value = -sum(lgamma(n+1) - lgamma(n-my_data+1) + lgamma(my_data + 1) + 
                     my_data*log(sum(my_data) / (m*n)) + (n - my_data)*log(1 - sum(my_data)/(m*n)))
      return(value)
    }
    n_hat = floor(max(optim(max(my_data), lower = max(my_data), upper = 1000, f_likelihood, method = "Brent")$par,
                      max(my_data)))
  } else if(starting_method == "MOM"){
    m_1 = mean(my_data)
    m_2 = mean(my_data^2)
    n_hat = max(floor( (-m_1^2) / (m_2 - m_1 - m_1^2)), max(my_data))
  } else if(starting_method == "MAX"){
    n_hat = max(my_data) + i
  } else if(starting_method == "Carroll-Lombard"){
    f_CL_likelihood = function(n) {
      value = sum(lgamma(n + 1) - lgamma(n - my_data + 1) - lgamma(my_data + 1) ) - 
        log(m*n + 2 + 1) - lgamma(m*n + 3) + lgamma(2 + sum(my_data)) +
        lgamma(m * n - sum(my_data) + 2)
      return(-value)
    }
    n_hat = floor(optim(max(my_data), lower = max(my_data), upper = 1000, f_CL_likelihood, method = "Brent")$par)
  } else if(starting_method == "DasGupta-Rubin"){
    n_1_hat = floor((max(my_data)^2 * var(my_data)) / (mean(my_data) * (max(my_data) - mean(my_data))))
    bias_quantity = 0
    for(i in 0:(n_1_hat - 2)){
      bias_quantity = bias_quantity + qbeta(1/length(my_data), shape1 = i + 1, shape2 = n_1_hat - i)
    }
    n_hat = floor(max(my_data) + bias_quantity)
  } else {
    if(is.null(user_n)){
      print("starting_method = 'USER' requires non-null user_n value.")
    }
    n_hat = user_n
  }
  
  generate_uniforms = generate_first_uniforms(my_data, n_hat, first_mu)
  mins = function(x){return(min(generate_uniforms[which(my_data == x)]))}
  maxs = function(x){return(max(generate_uniforms[which(my_data == x)]))}
  all_data = data.frame(value = unique_data, number = repeat_observations,
                        min_uniform = sapply(unique_data, mins),
                        max_uniform = sapply(unique_data, maxs))
  temp_all_data = all_data
  temp_all_data$iteration = rep(1, times = nrow(all_data))
  p_set_generation = generate_p_sets(all_data, 0, n_hat, epsilon)
  
  for(iter in 1:iterations){
    # print(paste("Starting iteration:", iter))
    epsilon_constant = p_set_generation$epsilon_constant
    temp_data = p_set_generation$pset_data
    solution_data = rbind(solution_data, temp_data)
    if(iter != iterations){
      
      generate_uniforms = update_uniforms(all_data, n_hat, epsilon, 
                                          ending_n_hat = p_set_generation$ending_n, 
                                          sigma, epsilon_constant)
      p_set_generation = generate_p_sets(all_data, iter+1, n_hat, epsilon,
                                         final_n = generate_uniforms$ending_n, 
                                         generate_uniforms$epsilon_constant)
      all_data = generate_uniforms$the_data
      
      temp_all_data = all_data
      temp_all_data$iteration = rep(iter+1, times = nrow(all_data))
    }
  }
  solution_data = solution_data[-1,]
  return(list(solution_data, n_hat, first_mu))
}