## Finding the Belief bounds
#### Start at estimator and grow outward.
# Author: Alexander Murph

belief_mass = function(posterior_draws, lower_mu, upper_mu, lower_n, upper_n, total_number){
  iterations = unique(posterior_draws$iteration_number)
  num_included = 0
  remaining_data=NULL
  for(index in iterations){
    temp_data = posterior_draws[which(posterior_draws$iteration_number == index),]
    mass_covered = (temp_data$n_value <= upper_n)&
      (temp_data$n_value >= lower_n)&
      (temp_data$upper_mu <= upper_mu)&
      (temp_data$lower_mu >= lower_mu)
    num_included = num_included + (sum(mass_covered)==nrow(temp_data))
    if((!(sum(mass_covered)==nrow(temp_data)))&(is.null(remaining_data))){
      remaining_data=temp_data
    } else if ((!(sum(mass_covered)==nrow(temp_data)))&(!(is.null(remaining_data)))){
      remaining_data=rbind(remaining_data, temp_data)
    }
  }
  return(list(num_included/total_number, remaining_data))
}


plausability_mass = function(posterior_draws, lower_mu, upper_mu, lower_n, upper_n, total_number){
  iterations = unique(posterior_draws$iteration_number)
  num_included = 0
  remaining_data=NULL
  for(index in iterations){
    temp_data = posterior_draws[which(posterior_draws$iteration_number == index),]
    n_fits = (temp_data$n_value <= upper_n)&(temp_data$n_value >= lower_n)
    mu_fits = ((lower_mu<=temp_data$upper_mu)&(upper_mu>=temp_data$upper_mu))|
      ((lower_mu<=temp_data$lower_mu)&(upper_mu>=temp_data$lower_mu))|
      ((lower_mu>=temp_data$lower_mu)&(upper_mu<=temp_data$upper_mu))
    mass_covered = mean( n_fits&mu_fits )
    num_included = num_included + (sum(mass_covered)>0)
    if((!(sum(mass_covered)>0))&(is.null(remaining_data))){
      remaining_data=temp_data
    } else if ((!(sum(mass_covered)>0))&(!(is.null(remaining_data)))){
      remaining_data=rbind(remaining_data, temp_data)
    }
  }
  return(list(num_included/total_number, remaining_data))
}


belief_bounds = function(posterior_draws){
  posterior_draws$upper_mu = posterior_draws$upper_p * posterior_draws$n_value
  posterior_draws$lower_mu = posterior_draws$lower_p * posterior_draws$n_value
  all_mus = sort(unique(c(posterior_draws$upper_mu, posterior_draws$lower_mu)))
  all_ns = sort(unique(c(posterior_draws$n_value)))
  
  # I'm going to 'boost' the value to hopefully quicken this process.
  qu <- quantile(all_mus, 0.975)
  ql <- quantile(all_mus, 0.025)
  canidate_mu_upper = min(all_mus[ all_mus >= qu ])
  canidate_mu_lower = max(all_mus[ all_mus <= ql ])
  qu <- quantile(all_ns, 0.6)
  ql <- quantile(all_ns, 0.4)
  canidate_n_upper = min(all_ns[ all_ns >= qu ])
  canidate_n_lower = max(all_ns[ all_ns <= ql ])

  full_data=posterior_draws
  total_number=length(unique(posterior_draws$iteration_number))
  index = 0
  four_corners = c(0,0,0,0)
  mass_covered=0
  while(TRUE){
    print(paste("index", index, "had belief mass:",mass_covered))
    print(all_ns)
    print("HEY")
    print(paste("--->candidate n lower:",canidate_n_lower))
    print(paste("--->candidate n upper:",canidate_n_upper))
    print(paste("--->candidate mu lower:",canidate_mu_lower))
    print(paste("--->candidate mu upper:",canidate_mu_upper))
    i = index%%4
    index = index + 1
    if(i == 0){
      four_corners = c(0,0,0,0)
      temp_n_upper = all_ns[which(all_ns == canidate_n_upper)+1]
      if(canidate_n_upper == max(all_ns)){
        next
      }
      if( (length(temp_n_upper) == 0) ){
        next
      }
      mass_covered_additionally = belief_mass(posterior_draws, canidate_mu_lower, 
                                 canidate_mu_upper, canidate_n_lower, 
                                 temp_n_upper, total_number)
      mass_covered_temp=mass_covered+mass_covered_additionally[[1]]
      if( (mass_covered_temp <= 0.95)&(temp_n_upper >= canidate_n_lower)  ){
        canidate_n_upper = temp_n_upper
        four_corners[1] = 1
        mass_covered=mass_covered+mass_covered_additionally[[1]]
        posterior_draws=mass_covered_additionally[[2]]
      } else if ((mass_covered_temp >= 0.95)&(temp_n_upper >= canidate_n_lower)){
        canidate_n_upper = temp_n_upper
        four_corners[1] = 1
        mass_covered=mass_covered+mass_covered_additionally[[1]]
        posterior_draws=mass_covered_additionally[[2]]
        break
      }
    } else if(i == 1) {
      if(canidate_mu_upper == max(all_mus)){
        next
      }
      temp_mu_upper = all_mus[which(all_mus == canidate_mu_upper)+1]
      
      mass_covered_additionally = belief_mass(posterior_draws, canidate_mu_lower, 
                                              temp_mu_upper, canidate_n_lower, 
                                              canidate_n_upper, total_number)
      mass_covered_temp=mass_covered+mass_covered_additionally[[1]]
      if((mass_covered_temp <= 0.95)&(temp_mu_upper >= canidate_mu_lower)){
        canidate_mu_upper = temp_mu_upper
        four_corners[2] = 1
        mass_covered=mass_covered+mass_covered_additionally[[1]]
        posterior_draws=mass_covered_additionally[[2]]
      } else if((mass_covered_temp >= 0.95)&(temp_mu_upper >= canidate_mu_lower)){
        canidate_mu_upper = temp_mu_upper
        four_corners[2] = 1
        mass_covered=mass_covered+mass_covered_additionally[[1]]
        posterior_draws=mass_covered_additionally[[2]]
        break
      }
    } else if(i == 2){
      temp_n_lower = all_ns[which(all_ns == canidate_n_lower)-1]
      if(canidate_n_lower == min(all_ns)){
        next
      }
      if(  (is.na(temp_n_lower))){
        next
      }
      
      mass_covered_additionally = belief_mass(posterior_draws, canidate_mu_lower, 
                                              canidate_mu_upper, temp_n_lower, 
                                              canidate_n_upper, total_number)
      mass_covered_temp=mass_covered+mass_covered_additionally[[1]]
      
      if((mass_covered_temp <= 0.95)&(temp_n_lower <= canidate_n_upper)){
        canidate_n_lower = temp_n_lower
        four_corners[3] = 1
        mass_covered=mass_covered+mass_covered_additionally[[1]]
        posterior_draws=mass_covered_additionally[[2]]
      } else if((mass_covered_temp >= 0.95)&(temp_n_lower <= canidate_n_upper)){
        canidate_n_lower = temp_n_lower
        four_corners[3] = 1
        mass_covered=mass_covered+mass_covered_additionally[[1]]
        posterior_draws=mass_covered_additionally[[2]]
        break
      }
    } else {
      if(canidate_mu_lower == min(all_mus)){
        next
      }
      temp_mu_lower = all_mus[which(all_mus == canidate_mu_lower)-1]
      
      mass_covered_additionally = belief_mass(posterior_draws, temp_mu_lower, 
                                              canidate_mu_upper, canidate_n_lower, 
                                              canidate_n_upper, total_number)
      mass_covered_temp=mass_covered+mass_covered_additionally[[1]]
      
      if((mass_covered_temp <= 0.95)&(temp_mu_lower <= canidate_mu_upper)){
        canidate_mu_lower = temp_mu_lower
        four_corners[4] = 1
        mass_covered=mass_covered+mass_covered_additionally[[1]]
        posterior_draws=mass_covered_additionally[[2]]
      } else if((mass_covered_temp >= 0.95)&(temp_mu_lower <= canidate_mu_upper)){
        canidate_mu_lower = temp_mu_lower
        four_corners[4] = 1
        mass_covered=mass_covered+mass_covered_additionally[[1]]
        posterior_draws=mass_covered_additionally[[2]]
        break
      }
      if(identical(four_corners,c(0,0,0,0))){
        break
      }
    }
    
  }
  # mass_covered = belief_mass(full_data, canidate_mu_lower,
  #                            canidate_mu_upper, canidate_n_lower,
  #                            canidate_n_upper,total_number)
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
  
  qu <- quantile(all_mus, 0.975)
  ql <- quantile(all_mus, 0.025)
  canidate_mu_upper = min(all_mus[ all_mus >= qu ])
  canidate_mu_lower = max(all_mus[ all_mus <= ql ])
  qu <- quantile(all_ns, 0.6)
  ql <- quantile(all_ns, 0.4)
  canidate_n_upper = min(all_ns[ all_ns >= qu ])
  canidate_n_lower = max(all_ns[ all_ns <= ql ])
  
  full_data=posterior_draws
  total_number=length(unique(posterior_draws$iteration_number))
  index = 0
  four_corners = c(0,0,0,0)
  mass_covered=0
  while(TRUE){
    print(paste("index", index, "had plaus mass:",mass_covered))
    print(paste("--->candidate n lower:",canidate_n_lower))
    print(paste("--->candidate n upper:",canidate_n_upper))
    print(paste("--->candidate mu lower:",canidate_mu_lower))
    print(paste("--->candidate mu upper:",canidate_mu_upper))
    i = index%%4
    index = index + 1
    if(i == 0){
      four_corners = c(0,0,0,0)
      temp_n_upper = all_ns[which(all_ns == canidate_n_upper)+1]
      if( (length(temp_n_upper) == 0) ){
        next
      }
      if(canidate_n_upper == max(all_ns)){
        next
      }
      mass_covered_additionally = plausability_mass(posterior_draws, canidate_mu_lower, 
                                       canidate_mu_upper, canidate_n_lower, 
                                       temp_n_upper, total_number)
      mass_covered_temp=mass_covered+mass_covered_additionally[[1]]
      
      if( (mass_covered_temp <= 0.95)&(temp_n_upper >= canidate_n_lower)  ){
        canidate_n_upper = temp_n_upper
        four_corners[1] = 1
        mass_covered=mass_covered+mass_covered_additionally[[1]]
        posterior_draws = mass_covered_additionally[[2]]
      } else if ((mass_covered_temp >= 0.95)&(temp_n_upper >= canidate_n_lower)){
        canidate_n_upper = temp_n_upper
        four_corners[1] = 1
        mass_covered=mass_covered+mass_covered_additionally[[1]]
        posterior_draws = mass_covered_additionally[[2]]
        break
      }
    } else if(i == 1) {
      if(canidate_mu_upper == max(all_mus)){
        next
      }
      temp_mu_upper = all_mus[which(all_mus == canidate_mu_upper)+1]
      
      mass_covered_additionally = plausability_mass(posterior_draws, canidate_mu_lower, 
                                                    temp_mu_upper, canidate_n_lower, 
                                                    canidate_n_upper, total_number)
      mass_covered_temp=mass_covered+mass_covered_additionally[[1]]
      
      if((mass_covered_temp <= 0.95)&(temp_mu_upper >= canidate_mu_lower)){
        canidate_mu_upper = temp_mu_upper
        four_corners[2] = 1
        mass_covered=mass_covered+mass_covered_additionally[[1]]
        posterior_draws = mass_covered_additionally[[2]]
      } else if ((mass_covered_temp >= 0.95)&(temp_mu_upper >= canidate_mu_lower)){
        canidate_mu_upper = temp_mu_upper
        four_corners[2] = 1
        mass_covered=mass_covered+mass_covered_additionally[[1]]
        posterior_draws = mass_covered_additionally[[2]]
        break
      }
    } else if(i == 2){
      if(canidate_n_lower == min(all_ns)){
        next
      }
      temp_n_lower = all_ns[which(all_ns == canidate_n_lower)-1]
      if(  (is.na(temp_n_lower))){
        next
      }
      
      mass_covered_additionally = plausability_mass(posterior_draws, canidate_mu_lower, 
                                                    canidate_mu_upper, temp_n_lower, 
                                                    canidate_n_upper, total_number)
      mass_covered_temp=mass_covered+mass_covered_additionally[[1]]
      
      if((mass_covered_temp <= 0.95)&(temp_n_lower <= canidate_n_upper)){
        canidate_n_lower = temp_n_lower
        four_corners[3] = 1
        mass_covered=mass_covered+mass_covered_additionally[[1]]
        posterior_draws=mass_covered_additionally[[2]]
      } else if((mass_covered_temp >= 0.95)&(temp_n_lower <= canidate_n_upper)) {
        canidate_n_lower = temp_n_lower
        four_corners[3] = 1
        mass_covered=mass_covered+mass_covered_additionally[[1]]
        posterior_draws=mass_covered_additionally[[2]]
        break
      }
    } else {
      if(canidate_mu_lower == min(all_mus)){
        next
      }
      temp_mu_lower = all_mus[which(all_mus == canidate_mu_lower)-1]
      
      mass_covered_additionally = plausability_mass(posterior_draws, temp_mu_lower, 
                                                    canidate_mu_upper, canidate_n_lower, 
                                                    canidate_n_upper, total_number)
      mass_covered_temp=mass_covered+mass_covered_additionally[[1]]
      
      if((mass_covered_temp <= 0.95)&(temp_mu_lower <= canidate_mu_upper)){
        canidate_mu_lower = temp_mu_lower
        four_corners[4] = 1
        mass_covered=mass_covered+mass_covered_additionally[[1]]
        posterior_draws=mass_covered_additionally[[2]]
      } else if((mass_covered_temp >= 0.95)&(temp_mu_lower <= canidate_mu_upper)){
        canidate_mu_lower = temp_mu_lower
        four_corners[4] = 1
        mass_covered=mass_covered+mass_covered_additionally[[1]]
        posterior_draws=mass_covered_additionally[[2]]
        break
      }
      if(identical(four_corners,c(0,0,0,0))){
        break
      }
    }
    
  }
  # mass_covered = plausability_mass(full_data, canidate_mu_lower, 
  #                                  canidate_mu_upper, canidate_n_lower, 
  #                                  canidate_n_upper, total_number)
  return(pairlist(mu_upper = canidate_mu_upper,
                  mu_lower = canidate_mu_lower,
                  n_upper = canidate_n_upper,
                  n_lower = canidate_n_lower,
                  mass = mass_covered))
}
