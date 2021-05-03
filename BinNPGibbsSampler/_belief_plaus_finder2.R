## Finding the Belief bounds
#### This one takes out PARTICLES following the four corners method from the first belief plaus file.
# Author: Alexander Murph

belief_mass = function(posterior_draws, lower_mu, upper_mu, lower_n, upper_n){
  iterations = unique(posterior_draws$iteration_number)
  num_included = 0
  total_number = length(iterations)
  removed_particles = c()
  for(index in iterations){
    temp_data = posterior_draws[which(posterior_draws$iteration_number == index),]
    mass_covered = (temp_data$n_value <= upper_n)&
      (temp_data$n_value >= lower_n)&
      (temp_data$upper_mu <= upper_mu)&
      (temp_data$lower_mu >= lower_mu)
    num_included = num_included + (sum(mass_covered)==nrow(temp_data))
    removed_particles = c(removed_particles, index)
  }
  return(c(num_included/total_number, removed_particles))
}


plausability_mass = function(posterior_draws, lower_mu, upper_mu, lower_n, upper_n){
  iterations = unique(posterior_draws$iteration_number)
  num_included = 0
  total_number = length(iterations)
  removed_particles = c()
  for(index in iterations){
    temp_data = posterior_draws[which(posterior_draws$iteration_number == index),]
    n_fits = (temp_data$n_value <= upper_n)&(temp_data$n_value >= lower_n)
    mu_fits = ((lower_mu<=temp_data$upper_mu)&(upper_mu>=temp_data$upper_mu))|
      ((lower_mu<=temp_data$lower_mu)&(upper_mu>=temp_data$lower_mu))|
      ((lower_mu>=temp_data$lower_mu)&(upper_mu<=temp_data$upper_mu))
    mass_covered = mean( n_fits&mu_fits )
    num_included = num_included + (sum(mass_covered)>0)
    removed_particles = c(removed_particles, index)
  }
  return(c(num_included/total_number, removed_particles))
}


belief_bounds = function(posterior_draws){
  posterior_draws$upper_mu = posterior_draws$upper_p * posterior_draws$n_value
  posterior_draws$lower_mu = posterior_draws$lower_p * posterior_draws$n_value
  all_mus = sort(unique(c(posterior_draws$upper_mu, posterior_draws$lower_mu)))
  all_ns = sort(unique(c(posterior_draws$n_value)))
  
  canidate_mu_upper = max(all_mus)
  canidate_mu_lower = min(all_mus)
  canidate_n_upper = max(all_ns)
  canidate_n_lower = min(all_ns)
  
  removed_particles = c()
  
  index = 0
  four_corners = c(0,0,0,0)
  while(TRUE){
    if(!is.null(removed_particles)){
      indices = which(posterior_draws$iteration_number == removed_particles)
      all_mus = sort(unique(c(posterior_draws$upper_mu[-indices], 
                              posterior_draws$lower_mu[-indices])))
      all_ns = sort(unique(c(posterior_draws$n_value[-indices])))
    }
    i = index%%4
    index = index + 1
    if(i == 0){
      four_corners = c(0,0,0,0)
      temp_n_upper = all_ns[which(all_ns == canidate_n_upper)-1]
      if( (length(temp_n_upper) == 0) ){
        next
      }
      mass_covered = belief_mass(posterior_draws, canidate_mu_lower, 
                                 canidate_mu_upper, canidate_n_lower, 
                                 temp_n_upper)
      removed_particles_temp = mass_covered[2:length(mass_covered)]
      mass_covered = mass_covered[1]
      
      if( (mass_covered >= 0.95)&(temp_n_upper >= canidate_n_lower)  ){
        canidate_n_upper = temp_n_upper
        four_corners[1] = 1
        removed_particles = removed_particles_temp
      } 
    } else if(i == 1) {
      temp_mu_upper = all_mus[which(all_mus == canidate_mu_upper)-1]
      
      mass_covered = belief_mass(posterior_draws, canidate_mu_lower, 
                                 temp_mu_upper, canidate_n_lower, 
                                 canidate_n_upper)
      removed_particles_temp = mass_covered[2:length(mass_covered)]
      mass_covered = mass_covered[1]
      
      if((mass_covered >= 0.95)&(temp_mu_upper >= canidate_mu_lower)){
        canidate_mu_upper = temp_mu_upper
        four_corners[2] = 1
        removed_particles = removed_particles_temp
      }
    } else if(i == 2){
      temp_n_lower = all_ns[which(all_ns == canidate_n_lower)+1]
      if(  (is.na(temp_n_lower))){
        next
      }
      
      mass_covered = belief_mass(posterior_draws, canidate_mu_lower, 
                                 canidate_mu_upper, temp_n_lower, 
                                 canidate_n_upper)
      removed_particles_temp = mass_covered[2:length(mass_covered)]
      mass_covered = mass_covered[1]
      
      if((mass_covered >= 0.95)&(temp_n_lower <= canidate_n_upper)){
        canidate_n_lower = temp_n_lower
        four_corners[3] = 1
        removed_particles = removed_particles_temp
      }
    } else {
      temp_mu_lower = all_mus[which(all_mus == canidate_mu_lower)+1]
      
      mass_covered = belief_mass(posterior_draws, temp_mu_lower, 
                                 canidate_mu_upper, canidate_n_lower, 
                                 canidate_n_upper)
      removed_particles_temp = mass_covered[2:length(mass_covered)]
      mass_covered = mass_covered[1]
      
      if((mass_covered >= 0.95)&(temp_mu_lower <= canidate_mu_upper)){
        canidate_mu_lower = temp_mu_lower
        four_corners[4] = 1
        removed_particles = removed_particles_temp
      }
      if(identical(four_corners,c(0,0,0,0))){
        break
      }
    }
    
  }
  mass_covered = belief_mass(posterior_draws, canidate_mu_lower, 
                             canidate_mu_upper, canidate_n_lower, 
                             canidate_n_upper)
  return(pairlist(mu_upper = canidate_mu_upper,
                  mu_lower = canidate_mu_lower,
                  n_upper = canidate_n_upper,
                  n_lower = canidate_n_lower,
                  mass = mass_covered))
}


plausability_bounds = function(posterior_draws){
  posterior_draws$upper_mu = posterior_draws$upper_p * posterior_draws$n_value
  posterior_draws$lower_mu = posterior_draws$lower_p * posterior_draws$n_value
  all_mus = sort(unique(c(posterior_draws$upper_mu, posterior_draws$lower_mu)))
  all_ns = sort(unique(c(posterior_draws$n_value)))
  removed_particles = c()
  
  canidate_mu_upper = max(all_mus)
  canidate_mu_lower = min(all_mus)
  canidate_n_upper = max(all_ns)
  canidate_n_lower = min(all_ns)
  index = 0
  four_corners = c(0,0,0,0)
  while(TRUE){
    if(!is.null(removed_particles)){
      indices = which(posterior_draws$iteration_number == removed_particles)
      all_mus = sort(unique(c(posterior_draws$upper_mu[-indices], 
                              posterior_draws$lower_mu[-indices])))
      all_ns = sort(unique(c(posterior_draws$n_value[-indices])))
    }
    i = index%%4
    index = index + 1
    if(i == 0){
      four_corners = c(0,0,0,0)
      temp_n_upper = all_ns[which(all_ns == canidate_n_upper)-1]
      if( (length(temp_n_upper) == 0) ){
        next
      }
      mass_covered = plausability_mass(posterior_draws, canidate_mu_lower, 
                                       canidate_mu_upper, canidate_n_lower, 
                                       temp_n_upper)
      removed_particles_temp = mass_covered[2:length(mass_covered)]
      mass_covered = mass_covered[1]
      
      if( (mass_covered >= 0.95)&(temp_n_upper >= canidate_n_lower)  ){
        canidate_n_upper = temp_n_upper
        four_corners[1] = 1
        removed_particles - removed_particles_temp
      } 
    } else if(i == 1) {
      temp_mu_upper = all_mus[which(all_mus == canidate_mu_upper)-1]
      
      mass_covered = plausability_mass(posterior_draws, canidate_mu_lower, 
                                       temp_mu_upper, canidate_n_lower, 
                                       canidate_n_upper)
      removed_particles_temp = mass_covered[2:length(mass_covered)]
      mass_covered = mass_covered[1]
      
      if((mass_covered >= 0.95)&(temp_mu_upper >= canidate_mu_lower)){
        canidate_mu_upper = temp_mu_upper
        four_corners[2] = 1
        removed_particles = removed_particles_temp
      }
    } else if(i == 2){
      temp_n_lower = all_ns[which(all_ns == canidate_n_lower)+1]
      if(  (is.na(temp_n_lower))){
        next
      }
      
      mass_covered = plausability_mass(posterior_draws, canidate_mu_lower, 
                                       canidate_mu_upper, temp_n_lower, 
                                       canidate_n_upper)
      removed_particles_temp = mass_covered[2:length(mass_covered)]
      mass_covered = mass_covered[1]
      
      if((mass_covered >= 0.95)&(temp_n_lower <= canidate_n_upper)){
        canidate_n_lower = temp_n_lower
        four_corners[3] = 1
        removed_particles = removed_particles_temp
      }
    } else {
      temp_mu_lower = all_mus[which(all_mus == canidate_mu_lower)+1]
      
      mass_covered = plausability_mass(posterior_draws, temp_mu_lower, 
                                       canidate_mu_upper, canidate_n_lower, 
                                       canidate_n_upper)
      removed_particles_temp = mass_covered[2:length(mass_covered)]
      mass_covered = mass_covered[1]
      
      if((mass_covered >= 0.95)&(temp_mu_lower <= canidate_mu_upper)){
        canidate_mu_lower = temp_mu_lower
        four_corners[4] = 1
        removed_particles = removed_particles_temp
      }
      if(identical(four_corners,c(0,0,0,0))){
        break
      }
    }
    
  }
  mass_covered = plausability_mass(posterior_draws, canidate_mu_lower, 
                                   canidate_mu_upper, canidate_n_lower, 
                                   canidate_n_upper)
  return(pairlist(mu_upper = canidate_mu_upper,
                  mu_lower = canidate_mu_lower,
                  n_upper = canidate_n_upper,
                  n_lower = canidate_n_lower,
                  mass = mass_covered))
}